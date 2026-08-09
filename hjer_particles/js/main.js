/**
 * HJER Particle Physics Visualizer - Main Application
 *
 * Interactive visualization of particle physics on the 600-cell:
 * - 5 cosets (inscribed 24-cells) with proper 5-coloring
 * - Gauge edge classification: U(1) / SU(2) / SU(3)
 * - Soliton creation with spectral splitting
 * - Kibble-Zurek cooling simulation
 *
 * Based on results from exp164-181 of the 600-cell theory.
 */

import { PolytopeGraph } from './polytope.js';
import { Rotation4D } from './rotation4d.js';
import { projectAll, stereographic } from './projection.js';
import { Renderer } from './renderer.js';
import { computeCosets, verifyCosets } from './cosets.js';
import { classifyEdges } from './gauge.js';
import { SolitonManager } from './soliton.js';
import { PottsModel } from './potts.js';
import { drawSpectrum } from './spectrum.js';
import { UI } from './ui.js';
import { cosetColor, COSET_COLORS, VERTEX_SIZES, GAUGE_COLORS } from './colormap.js';
import { TOUR_STEPS } from './tour.js';

// Vertex type colors for "gauge" mode
const TYPE_COLORS = {
    A: [1.0, 1.0, 1.0],     // white (gauge vertices, 8 total)
    B: [0.6, 0.6, 0.7],     // gray (Higgs, 16 total)
    C: [0.25, 0.25, 0.35],  // dim (fermion, 96 total)
};

// Even dimmer type colors
const DIM_TYPE_COLORS = {
    A: [0.5, 0.5, 0.55],
    B: [0.35, 0.35, 0.4],
    C: [0.15, 0.15, 0.2],
};

class ParticlePhysicsApp {
    constructor() {
        this.polytope = new PolytopeGraph();
        this.rotation4d = new Rotation4D();
        this.renderer = new Renderer(document.getElementById('viewport'));
        this.ui = null;

        // Physics state
        this.coset = null;        // per-vertex coset (0-4)
        this.types = null;        // per-vertex type ('A', 'B', 'C')
        this.edgeType = null;     // per-edge gauge type
        this.edgeSubType = null;
        this.gaugeCounts = null;

        this.soliton = null;      // SolitonManager
        this.potts = null;        // PottsModel

        // Working buffers
        this.rotated4D = null;
        this.projected3D = null;
        this.visible = null;
        this.vertexColors = null;
        this.vertexScales = null;

        // Display state
        this.mode = 'tour';        // 'tour' | 'coset' | 'gauge' | 'soliton' | 'cooling'
        this.autoRotate = false;
        this.autoRotateSpeed = 0.005;
        this.cosetVisible = [true, true, true, true, true];
        this.gaugeVisible = { U1: true, SU2: true, SU3: true };
        this.flipMode = false;
        this.pairMode = false;
        this._pairFirst = -1;

        // Cooling animation
        this._coolingActive = false;
        this._coolSchedule = null;
        this._coolStep = 0;
        this._coolInterval = null;

        // Soliton animation
        this._glowPhase = 0;
        this._lastFlipResult = null;
        this._splitAnimProgress = 0;

        // Tour state
        this._tourStep = 0;
        this._tourFeaturedVertex = -1;
        this._tourNeighborhood = null; // Set of vertex indices
        this._tourEdgeCounts = null;   // { U1, SU2, SU3 } for featured vertex

        // Spectral canvas
        this.spectrumCanvas = null;
        this.spectrumCtx = null;

        // FPS
        this._frameCount = 0;
        this._lastFpsTime = performance.now();
        this._fps = 0;

        this.running = false;
    }

    async init() {
        this.ui = new UI(this);

        // Welcome overlay
        this._initWelcome();

        // Create spectral canvas
        this._createSpectrumCanvas();

        // Legend
        this.legendEl = document.getElementById('legend');
        this._updateLegend();

        // Event handling
        this.renderer.container.addEventListener('click', (e) => this._onClick(e));

        await this._loadPolytope();
        this.running = true;
        this._animate();
    }

    _initWelcome() {
        const overlay = document.getElementById('welcome-overlay');
        const closeBtn = document.getElementById('welcome-close');
        if (overlay && closeBtn) {
            closeBtn.addEventListener('click', () => {
                overlay.classList.add('hidden');
            });
        }
    }

    async _loadPolytope() {
        await this.polytope.load('data/600cell.json');

        const n = this.polytope.numVertices;
        this.rotated4D = new Float32Array(n * 4);
        this.projected3D = new Float32Array(n * 3);
        this.visible = new Uint8Array(n);
        this.vertexColors = new Array(n);
        this.vertexScales = new Float32Array(n);

        // Compute cosets
        const result = computeCosets(this.polytope.vertices4D, this.polytope.adjacency);
        this.coset = result.coset;
        this.types = result.types;

        // Verify
        const verify = verifyCosets(this.coset, this.types, this.polytope.edges);
        console.log('Coset verification:', verify);

        // Classify edges
        const edgeResult = classifyEdges(
            this.polytope.edges,
            this.polytope.adjacency,
            this.types
        );
        this.edgeType = edgeResult.edgeType;
        this.edgeSubType = edgeResult.edgeSubType;
        this.gaugeCounts = edgeResult.counts;

        console.log('Gauge counts:', this.gaugeCounts);
        console.log('Sub-type counts:', edgeResult.subCounts);

        // Init renderer
        this.renderer.initPolytope(n, this.polytope.numEdges, this.edgeType);

        // Init soliton manager
        this.soliton = new SolitonManager(this.coset, this.polytope.adjacency, this.polytope.edges);

        // Init Potts model
        this.potts = new PottsModel(n, this.polytope.adjacency, this.polytope.edges);
        this.potts.initFromCoset(this.coset);

        // Compute vertex scales from types
        for (let i = 0; i < n; i++) {
            this.vertexScales[i] = VERTEX_SIZES[this.types[i]] || 1.0;
        }

        // Init tour
        this._initTour();

        // Set initial colors and show tour step 0
        this._updateColors();
        this._updateEdgeVisuals();
        this.setTourStep(0);
    }

    /**
     * Initialize tour: pick featured vertex, compute its neighborhood.
     */
    _initTour() {
        // Pick a C-type vertex in coset 1 (blue) for demonstrations
        for (let i = 0; i < this.polytope.numVertices; i++) {
            if (this.types[i] === 'C' && this.coset[i] === 1) {
                this._tourFeaturedVertex = i;
                break;
            }
        }

        // Compute neighborhood (featured + 12 neighbors)
        const nbrs = this.polytope.adjacency[this._tourFeaturedVertex];
        this._tourNeighborhood = new Set([this._tourFeaturedVertex, ...nbrs]);

        // Count edges by gauge type for the featured vertex
        this._tourEdgeCounts = { U1: 0, SU2: 0, SU3: 0 };
        for (let e = 0; e < this.polytope.edges.length; e++) {
            const [a, b] = this.polytope.edges[e];
            if (a === this._tourFeaturedVertex || b === this._tourFeaturedVertex) {
                this._tourEdgeCounts[this.edgeType[e]]++;
            }
        }
        // Expected: U1=1, SU2=3, SU3=8
        console.log('Tour featured vertex:', this._tourFeaturedVertex,
                    'Edge counts:', this._tourEdgeCounts);
    }

    _createSpectrumCanvas() {
        const specPanel = document.getElementById('spectrum-panel');
        if (!specPanel) return;

        this.spectrumCanvas = document.createElement('canvas');
        const dpr = window.devicePixelRatio || 1;
        const w = specPanel.clientWidth || 280;
        const h = specPanel.clientHeight || 400;
        this.spectrumCanvas.width = w * dpr;
        this.spectrumCanvas.height = h * dpr;
        this.spectrumCanvas.style.width = '100%';
        this.spectrumCanvas.style.height = '100%';
        specPanel.appendChild(this.spectrumCanvas);
        this.spectrumCtx = this.spectrumCanvas.getContext('2d');
        this.spectrumCtx.scale(dpr, dpr);
        this._specW = w;
        this._specH = h;
    }

    // --- Mode management ---

    setMode(mode) {
        const prevMode = this.mode;
        this.mode = mode;
        this.flipMode = false;
        this.pairMode = false;
        this._pairFirst = -1;

        // Auto-randomize when entering cooling mode (so it looks different)
        if (mode === 'cooling' && prevMode !== 'cooling') {
            this.potts.randomize();
        }

        // Reset soliton if leaving tour (soliton step may have created one)
        if (prevMode === 'tour' && this.soliton) {
            this.soliton.reset();
        }

        // Show tour step if entering tour mode
        if (mode === 'tour') {
            this.setTourStep(this._tourStep);
        }

        this._updateColors();
        this._updateEdgeVisuals();
        this._updateLegend();
    }

    // --- Tour navigation ---

    setTourStep(step) {
        const prevStep = this._tourStep;
        this._tourStep = Math.max(0, Math.min(step, TOUR_STEPS.length - 1));
        const tourStep = TOUR_STEPS[this._tourStep];

        // Handle soliton step (index 9): create/reset soliton
        if (this._tourStep === 9 && prevStep !== 9 && this.soliton) {
            this.soliton.reset();
            const target = this.soliton.pickTargetCoset(this._tourFeaturedVertex);
            this.soliton.flipVertex(this._tourFeaturedVertex, target);
            this._splitAnimProgress = 0;
        } else if (prevStep === 9 && this._tourStep !== 9 && this.soliton) {
            this.soliton.reset();
        }

        this._updateColors();
        this._updateEdgeVisuals();
        this._updateLegend();

        // Update UI panel
        if (this.ui) {
            this.ui.updateTourStep(
                this._tourStep, TOUR_STEPS.length,
                tourStep.title, tourStep.text
            );
        }
    }

    nextTourStep() {
        this.setTourStep(this._tourStep + 1);
    }

    prevTourStep() {
        this.setTourStep(this._tourStep - 1);
    }

    /**
     * Update the viewport legend overlay based on current mode.
     */
    _updateLegend() {
        if (!this.legendEl) return;

        const legends = {
            coset: `
                <div class="legend-title">Cosets (5 x 24-cell)</div>
                <div class="legend-row"><span class="legend-dot" style="background:#e64033"></span><span class="legend-label">Coset 0 - Red (A+B = 24-cell = 2T)</span></div>
                <div class="legend-row"><span class="legend-dot" style="background:#3380f2"></span><span class="legend-label">Coset 1 - Blue</span></div>
                <div class="legend-row"><span class="legend-dot" style="background:#33cc4d"></span><span class="legend-label">Coset 2 - Green</span></div>
                <div class="legend-row"><span class="legend-dot" style="background:#f2c019"></span><span class="legend-label">Coset 3 - Gold</span></div>
                <div class="legend-row"><span class="legend-dot" style="background:#b34de6"></span><span class="legend-label">Coset 4 - Purple</span></div>
                <div class="legend-row"><span class="legend-line" style="background:#ffe633"></span><span class="legend-label">U(1) 48 edges</span></div>
                <div class="legend-row"><span class="legend-line" style="background:#33ccff"></span><span class="legend-label">SU(2) 288 edges</span></div>
                <div class="legend-row"><span class="legend-line" style="background:#ff4d33"></span><span class="legend-label">SU(3) 384 edges</span></div>
                <div class="legend-hint">All 720 edges inter-coset (proper 5-coloring)</div>
            `,
            gauge: `
                <div class="legend-title">Gauge Channels</div>
                <div class="legend-row"><span class="legend-dot" style="background:#fff"></span><span class="legend-label">A-type: 8 gauge vertices (large)</span></div>
                <div class="legend-row"><span class="legend-dot" style="background:#999ab3"></span><span class="legend-label">B-type: 16 Higgs vertices (medium)</span></div>
                <div class="legend-row"><span class="legend-dot" style="background:#404059"></span><span class="legend-label">C-type: 96 fermion vertices (small)</span></div>
                <div class="legend-row"><span class="legend-line" style="background:#ffe633"></span><span class="legend-label">U(1) - 48 CC edges, perfect matching</span></div>
                <div class="legend-row"><span class="legend-line" style="background:#33ccff"></span><span class="legend-label">SU(2) - 288 edges (AC + BC)</span></div>
                <div class="legend-row"><span class="legend-line" style="background:#ff4d33"></span><span class="legend-label">SU(3) - 384 CC edges (with plaquettes)</span></div>
                <div class="legend-hint">1 + 3 + 8 = 12 = degree = dim(SM)</div>
            `,
            soliton: `
                <div class="legend-title">Soliton = Topological Defect</div>
                <div class="legend-row"><span class="legend-dot" style="background:#fff; box-shadow:0 0 6px #ff8"></span><span class="legend-label">Soliton (pulsing) = color was changed</span></div>
                <div class="legend-row"><span class="legend-line" style="background:#fff"></span><span class="legend-label">Broken edge = two neighbors now same color</span></div>
                <div class="legend-hint">Normal state: no two neighbors share color (E=0).</div>
                <div class="legend-hint">One flip always breaks exactly 3 edges (E=3).</div>
                <div class="legend-hint">A pair of adjacent flips: E=4 (not 6!) = binding.</div>
            `,
            cooling: `
                <div class="legend-title">Kibble-Zurek Cooling</div>
                <div class="legend-row"><span class="legend-dot" style="background:#e64033"></span><span class="legend-dot" style="background:#3380f2"></span><span class="legend-dot" style="background:#33cc4d"></span><span class="legend-dot" style="background:#f2c019"></span><span class="legend-dot" style="background:#b34de6"></span><span class="legend-label">5 Potts states (random)</span></div>
                <div class="legend-hint">Randomized at T=inf. Cool to see defects freeze.</div>
                <div class="legend-hint">Expect ~88 frozen defects (topological relics).</div>
                <div class="legend-hint">Matter = relict of cosmological phase transition.</div>
            `,
        };

        if (this.mode === 'tour') {
            this.legendEl.innerHTML = this._getTourLegend();
        } else {
            this.legendEl.innerHTML = legends[this.mode] || '';
        }
    }

    _getTourLegend() {
        const step = TOUR_STEPS[this._tourStep];
        let html = `<div class="legend-title">Tour: ${step.title}</div>`;

        if (step.highlight === 'neighborhood') {
            html += `
                <div class="legend-row"><span class="legend-dot" style="background:#ff8; box-shadow:0 0 6px #ff8"></span><span class="legend-label">Featured fermion (pulsing)</span></div>
                <div class="legend-row"><span class="legend-dot" style="background:#fff"></span><span class="legend-label">A-type neighbor (1)</span></div>
                <div class="legend-row"><span class="legend-dot" style="background:#999"></span><span class="legend-label">B-type neighbors (2)</span></div>
                <div class="legend-row"><span class="legend-dot" style="background:#404059"></span><span class="legend-label">C-type neighbors (9)</span></div>
            `;
            if (step.edges === 'neighborhood_gauge') {
                html += `
                    <div class="legend-row"><span class="legend-line" style="background:#ffe633"></span><span class="legend-label">U(1): 1 edge</span></div>
                    <div class="legend-row"><span class="legend-line" style="background:#33ccff"></span><span class="legend-label">SU(2): 3 edges</span></div>
                    <div class="legend-row"><span class="legend-line" style="background:#ff4d33"></span><span class="legend-label">SU(3): 8 edges</span></div>
                `;
            }
        } else if (step.edges === 'U1') {
            html += `
                <div class="legend-row"><span class="legend-line" style="background:#ffe633"></span><span class="legend-label">U(1): 48 edges (perfect matching)</span></div>
                <div class="legend-hint">Every fermion has exactly 1 yellow edge.</div>
            `;
        } else if (step.edges === 'SU3') {
            html += `
                <div class="legend-row"><span class="legend-line" style="background:#ff4d33"></span><span class="legend-label">SU(3): 384 edges (288 plaquettes)</span></div>
                <div class="legend-hint">The ONLY gauge type with plaquettes.</div>
            `;
        } else if (step.highlight === 'soliton') {
            html += `
                <div class="legend-row"><span class="legend-dot" style="background:#ff8; box-shadow:0 0 6px #ff8"></span><span class="legend-label">Soliton (flipped vertex)</span></div>
                <div class="legend-hint">Spectrum splits into 3 sub-levels.</div>
            `;
        } else {
            html += `<div class="legend-hint">Step ${this._tourStep + 1} of ${TOUR_STEPS.length}</div>`;
        }

        return html;
    }

    /**
     * Update edge layer opacity/visibility based on mode.
     */
    _updateEdgeVisuals() {
        if (!this.renderer.edgeLayers) return;

        if (this.mode === 'tour') {
            const step = TOUR_STEPS[this._tourStep];
            if (step.edges === 'U1') {
                // Only U(1) edges visible, bright
                this.renderer.edgeLayers['U1'].visible = true;
                this.renderer.edgeLayers['U1'].material.opacity = 0.95;
                this.renderer.edgeLayers['SU2'].visible = false;
                this.renderer.edgeLayers['SU3'].visible = false;
            } else if (step.edges === 'SU3') {
                // Only SU(3) edges visible
                this.renderer.edgeLayers['U1'].visible = false;
                this.renderer.edgeLayers['SU2'].visible = false;
                this.renderer.edgeLayers['SU3'].visible = true;
                this.renderer.edgeLayers['SU3'].material.opacity = 0.7;
            } else if (step.edges === 'dim') {
                // All edges very dim
                for (const type of ['U1', 'SU2', 'SU3']) {
                    this.renderer.edgeLayers[type].visible = true;
                    this.renderer.edgeLayers[type].material.opacity = 0.12;
                }
            } else if (step.edges === 'neighborhood_dim' || step.edges === 'neighborhood_gauge') {
                // Neighborhood: gauge colors, bright
                for (const type of ['U1', 'SU2', 'SU3']) {
                    this.renderer.edgeLayers[type].visible = true;
                    this.renderer.edgeLayers[type].material.opacity =
                        type === 'U1' ? 0.95 : (type === 'SU2' ? 0.8 : 0.65);
                }
            } else {
                // Default: gauge colors, moderate
                for (const type of ['U1', 'SU2', 'SU3']) {
                    this.renderer.edgeLayers[type].visible = true;
                    this.renderer.edgeLayers[type].material.opacity =
                        type === 'U1' ? 0.7 : (type === 'SU2' ? 0.45 : 0.35);
                }
            }
        } else if (this.mode === 'gauge') {
            // Gauge mode: edges bright and prominent
            for (const type of ['U1', 'SU2', 'SU3']) {
                const layer = this.renderer.edgeLayers[type];
                if (layer) {
                    layer.visible = this.gaugeVisible[type];
                    layer.material.opacity = type === 'U1' ? 0.9 : (type === 'SU2' ? 0.7 : 0.55);
                }
            }
        } else if (this.mode === 'cooling') {
            // Cooling: all edges dim uniform gray
            for (const type of ['U1', 'SU2', 'SU3']) {
                const layer = this.renderer.edgeLayers[type];
                if (layer) {
                    layer.visible = true;
                    layer.material.opacity = 0.15;
                }
            }
        } else {
            // Coset / soliton: normal gauge colors with moderate opacity
            for (const type of ['U1', 'SU2', 'SU3']) {
                const layer = this.renderer.edgeLayers[type];
                if (layer) {
                    layer.visible = this.gaugeVisible[type];
                    layer.material.opacity = type === 'U1' ? 0.7 : (type === 'SU2' ? 0.45 : 0.35);
                }
            }
        }
    }

    // --- Toggle handlers ---

    toggleCoset(cosetIdx, visible) {
        this.cosetVisible[cosetIdx] = visible;
    }

    toggleGauge(gaugeType, visible) {
        this.gaugeVisible[gaugeType] = visible;
        this.renderer.setEdgeLayerVisible(gaugeType, visible);
    }

    toggleFlipMode() {
        this.flipMode = !this.flipMode;
        if (this.flipMode) this.pairMode = false;
        this._pairFirst = -1;
    }

    togglePairMode() {
        this.pairMode = !this.pairMode;
        if (this.pairMode) this.flipMode = false;
        this._pairFirst = -1;
    }

    // --- Soliton operations ---

    /**
     * Auto-create a single soliton on a random C-type vertex (for demo).
     * Resets first so we start from clean ground state.
     */
    createDemoSoliton() {
        this.soliton.reset();
        // Pick a random C-type vertex (fermion) for demo
        const cVerts = [];
        for (let i = 0; i < this.polytope.numVertices; i++) {
            if (this.types[i] === 'C') cVerts.push(i);
        }
        const idx = cVerts[Math.floor(Math.random() * cVerts.length)];
        const target = this.soliton.pickTargetCoset(idx);
        const result = this.soliton.flipVertex(idx, target);
        if (result) {
            this._lastFlipResult = result;
            this._splitAnimProgress = 0;
        }
        this._updateColors();
    }

    /**
     * Auto-create a soliton pair on two adjacent C-type vertices (for demo).
     */
    createDemoPair() {
        this.soliton.reset();
        // Find a pair of adjacent C-type vertices
        for (const [i, j] of this.polytope.edges) {
            if (this.types[i] === 'C' && this.types[j] === 'C') {
                const result = this.soliton.createPair(i, j);
                if (result) {
                    this._lastFlipResult = { energy: result.energy, vertex: i };
                    this._splitAnimProgress = 0;
                }
                this._updateColors();
                return;
            }
        }
    }

    resetSolitons() {
        this.soliton.reset();
        this._lastFlipResult = null;
        this._splitAnimProgress = 0;
        this._updateColors();
    }

    // --- Potts operations ---

    setTemperature(T) {
        this.potts.temperature = T;
    }

    randomizePotts() {
        this.potts.randomize();
        this._stopCooling();
        this._updateColors();
    }

    setGroundState() {
        this.potts.initFromCoset(this.coset);
        this._stopCooling();
        this._updateColors();
    }

    startCooling() {
        if (this._coolingActive) {
            this._stopCooling();
            this.ui.coolBtn.textContent = 'Start Cooling';
            return;
        }

        this._coolingActive = true;
        this.ui.coolBtn.textContent = 'Stop Cooling';
        this.ui.coolBtn.classList.add('active');

        const speed = this.ui.getCoolSpeed();
        this._coolSchedule = PottsModel.coolingSchedule(this.potts.temperature, 0.01, 200);
        this._coolStep = 0;

        this._coolInterval = setInterval(() => {
            if (this._coolStep >= this._coolSchedule.length) {
                this._stopCooling();
                return;
            }

            this.potts.temperature = this._coolSchedule[this._coolStep];
            for (let i = 0; i < speed; i++) {
                this.potts.sweep();
            }
            this._coolStep++;
            this._updateColors();
        }, 50);
    }

    _stopCooling() {
        this._coolingActive = false;
        if (this._coolInterval) {
            clearInterval(this._coolInterval);
            this._coolInterval = null;
        }
        if (this.ui.coolBtn) {
            this.ui.coolBtn.textContent = 'Start Cooling';
            this.ui.coolBtn.classList.remove('active');
        }
    }

    // --- Click handling ---

    _onClick(e) {
        if (this.mode === 'soliton' && (this.flipMode || this.pairMode)) {
            const idx = this.renderer.pickVertex(e, this.projected3D, this.polytope.numVertices);
            if (idx < 0) return;

            if (this.flipMode) {
                const target = this.soliton.pickTargetCoset(idx);
                const result = this.soliton.flipVertex(idx, target);
                if (result) {
                    this._lastFlipResult = result;
                    this._splitAnimProgress = 0;
                    console.log(`Soliton at v${idx}: E=${result.energy}, coset ${result.oldCoset}->${result.newCoset}`);
                }
            } else if (this.pairMode) {
                if (this._pairFirst < 0) {
                    this._pairFirst = idx;
                    console.log(`Pair: first vertex = ${idx}. Click adjacent vertex...`);
                } else {
                    const result = this.soliton.createPair(this._pairFirst, idx);
                    if (result) {
                        this._lastFlipResult = { energy: result.energy, vertex: this._pairFirst };
                        this._splitAnimProgress = 0;
                        console.log(`Pair created: v${this._pairFirst}-v${idx}, E=${result.energy}`);
                    } else {
                        console.log(`Vertices ${this._pairFirst} and ${idx} are not adjacent`);
                    }
                    this._pairFirst = -1;
                }
            }

            this._updateColors();
        }
    }

    // --- Color computation (mode-dependent) ---

    _updateColors() {
        const n = this.polytope.numVertices;
        if (!n) return;

        switch (this.mode) {
            case 'coset':
                // Vertices bright by coset color
                for (let i = 0; i < n; i++) {
                    this.vertexColors[i] = cosetColor(this.coset[i]);
                }
                break;

            case 'gauge':
                // Vertices colored by TYPE: A=white, B=gray, C=dim
                // This makes the gauge edges the visual focus
                for (let i = 0; i < n; i++) {
                    this.vertexColors[i] = TYPE_COLORS[this.types[i]] || [0.3, 0.3, 0.3];
                }
                break;

            case 'soliton':
                // Vertices by current (possibly flipped) coset
                for (let i = 0; i < n; i++) {
                    this.vertexColors[i] = cosetColor(this.soliton.currentCoset[i]);
                }
                break;

            case 'cooling':
                // Vertices by Potts state (randomized / cooling)
                const state = this.potts.getState();
                for (let i = 0; i < n; i++) {
                    this.vertexColors[i] = cosetColor(state[i]);
                }
                break;

            case 'tour':
                this._updateTourColors(n);
                break;
        }
    }

    _updateTourColors(n) {
        const step = TOUR_STEPS[this._tourStep];

        if (step.vertices === 'coset') {
            for (let i = 0; i < n; i++) {
                this.vertexColors[i] = cosetColor(this.coset[i]);
            }
        } else if (step.vertices === 'type') {
            for (let i = 0; i < n; i++) {
                this.vertexColors[i] = TYPE_COLORS[this.types[i]];
            }
        } else if (step.vertices === 'dim_type') {
            for (let i = 0; i < n; i++) {
                this.vertexColors[i] = DIM_TYPE_COLORS[this.types[i]];
            }
        } else if (step.vertices === 'neighborhood') {
            // Non-neighborhood very dim, neighborhood bright by type
            for (let i = 0; i < n; i++) {
                if (this._tourNeighborhood && this._tourNeighborhood.has(i)) {
                    if (i === this._tourFeaturedVertex) {
                        this.vertexColors[i] = [1.0, 1.0, 0.5]; // bright yellow
                    } else {
                        this.vertexColors[i] = TYPE_COLORS[this.types[i]];
                    }
                } else {
                    this.vertexColors[i] = [0.08, 0.08, 0.12];
                }
            }
        } else if (step.vertices === 'soliton') {
            for (let i = 0; i < n; i++) {
                this.vertexColors[i] = cosetColor(this.soliton.currentCoset[i]);
            }
        } else {
            for (let i = 0; i < n; i++) {
                this.vertexColors[i] = cosetColor(this.coset[i]);
            }
        }
    }

    // --- Animation loop ---

    _animate() {
        if (!this.running) return;
        requestAnimationFrame(() => this._animate());

        const n = this.polytope.numVertices;
        if (n === 0) { this.renderer.render(); return; }

        // --- View rotation ---
        if (this.autoRotate) {
            this.rotation4d.addAngle('xw', this.autoRotateSpeed);
            this.rotation4d.addAngle('yz', this.autoRotateSpeed * 0.7);
            this.ui.syncSliders(this.rotation4d.angles);
        }

        // --- Project vertices ---
        this.rotation4d.transformAll(this.polytope.vertices4D, this.rotated4D);
        this.visible = projectAll(this.rotated4D, this.projected3D, n);

        // --- Apply coset visibility filter ---
        const effectiveVisible = new Uint8Array(this.visible);
        if (this.mode === 'coset' || this.mode === 'soliton') {
            for (let i = 0; i < n; i++) {
                const c = this.coset[i];
                if (!this.cosetVisible[c]) {
                    effectiveVisible[i] = 0;
                }
            }
        }
        // In cooling mode, filter by current Potts state
        if (this.mode === 'cooling') {
            const state = this.potts.getState();
            for (let i = 0; i < n; i++) {
                if (!this.cosetVisible[state[i]]) {
                    effectiveVisible[i] = 0;
                }
            }
        }

        // --- Glow animation for solitons ---
        this._glowPhase += 0.05;
        let glowSet = null;
        if (this.mode === 'soliton' && this.soliton) {
            const flipped = this.soliton.getFlippedVertices();
            if (flipped.length > 0) {
                glowSet = new Set(flipped);
            }
        }

        // Tour: glow featured vertex on neighborhood/soliton steps
        if (this.mode === 'tour' && this._tourFeaturedVertex >= 0) {
            const tourStep = TOUR_STEPS[this._tourStep];
            if (tourStep.highlight === 'neighborhood' || tourStep.highlight === 'soliton') {
                glowSet = new Set([this._tourFeaturedVertex]);
                // Also glow soliton-flipped vertices on soliton step
                if (tourStep.highlight === 'soliton' && this.soliton) {
                    for (const v of this.soliton.getFlippedVertices()) {
                        glowSet.add(v);
                    }
                }
            }
        }

        // --- Splitting animation progress ---
        if (this._splitAnimProgress < 1) {
            this._splitAnimProgress = Math.min(1, this._splitAnimProgress + 0.03);
        }

        // --- Vertex scales ---
        let scales = this.vertexScales;
        if (this.mode === 'gauge') {
            // Exaggerate type differences in gauge mode
            const gaugeScales = new Float32Array(n);
            for (let i = 0; i < n; i++) {
                if (this.types[i] === 'A') gaugeScales[i] = 2.5;
                else if (this.types[i] === 'B') gaugeScales[i] = 1.8;
                else gaugeScales[i] = 0.8;
            }
            scales = gaugeScales;
        } else if (this.mode === 'tour') {
            const tourStep = TOUR_STEPS[this._tourStep];
            if (tourStep.highlight === 'neighborhood') {
                // Neighborhood: bright+large. Others: tiny.
                const tourScales = new Float32Array(n);
                for (let i = 0; i < n; i++) {
                    if (this._tourNeighborhood && this._tourNeighborhood.has(i)) {
                        if (i === this._tourFeaturedVertex) {
                            tourScales[i] = 2.2;
                        } else {
                            tourScales[i] = (VERTEX_SIZES[this.types[i]] || 1.0) * 1.5;
                        }
                    } else {
                        tourScales[i] = 0.15;
                    }
                }
                scales = tourScales;
            } else if (tourStep.vertices === 'type' || tourStep.vertices === 'dim_type') {
                // Type-based sizes
                const typeScales = new Float32Array(n);
                for (let i = 0; i < n; i++) {
                    if (this.types[i] === 'A') typeScales[i] = 2.2;
                    else if (this.types[i] === 'B') typeScales[i] = 1.6;
                    else typeScales[i] = 0.9;
                }
                scales = typeScales;
            }
        }

        // --- Update renderer ---
        this.renderer.updateVertices(
            this.projected3D, effectiveVisible, n,
            this.vertexColors, scales,
            glowSet, this._glowPhase
        );

        this.renderer.updateEdges(
            this.polytope.edges, this.edgeType,
            this.projected3D, effectiveVisible
        );

        // Broken edges (soliton mode only)
        if (this.mode === 'soliton' && this.soliton) {
            const broken = this.soliton.getBrokenEdges();
            this.renderer.updateBrokenEdges(
                broken, this.polytope.edges,
                this.projected3D, effectiveVisible
            );
        } else {
            this.renderer.updateBrokenEdges(null, null, null, null);
        }

        // --- Spectral panel ---
        if (this.spectrumCtx) {
            const solitonActive =
                (this.mode === 'soliton' && this.soliton && this.soliton.solitonCount > 0) ||
                (this.mode === 'tour' && this._tourStep === 9 && this.soliton && this.soliton.solitonCount > 0);
            drawSpectrum(
                this.spectrumCtx,
                this._specW || 280,
                this._specH || 400,
                solitonActive,
                this.soliton ? this.soliton.solitonCount : 0,
                this._splitAnimProgress
            );
        }

        // --- FPS ---
        this._frameCount++;
        const now = performance.now();
        if (now - this._lastFpsTime >= 1000) {
            this._fps = this._frameCount;
            this._frameCount = 0;
            this._lastFpsTime = now;
        }

        // --- Update UI info ---
        const typeCounts = { A: 0, B: 0, C: 0 };
        for (const t of this.types) typeCounts[t]++;

        this.ui.updateInfo({
            vertices: n,
            edges: this.polytope.numEdges,
            fps: this._fps,
            mode: this.mode,
            cosetValid: true,
            gaugeU1: this.gaugeCounts.U1,
            gaugeSU2: this.gaugeCounts.SU2,
            gaugeSU3: this.gaugeCounts.SU3,
            typeA: typeCounts.A,
            typeB: typeCounts.B,
            typeC: typeCounts.C,
        });

        // Mode-specific info
        if (this.mode === 'soliton' && this.soliton) {
            this.ui.updateSolitonInfo({
                count: this.soliton.solitonCount,
                totalEnergy: this.soliton.computeTotalEnergy(),
                brokenEdges: this.soliton.getBrokenEdges().length,
                lastFlip: this._lastFlipResult,
            });
        }

        if (this.mode === 'cooling') {
            this.ui.updateCoolingInfo({
                temperature: this.potts.temperature,
                energy: this.potts.energy,
                defects: this.potts.defectCount,
                sweeps: this.potts.sweepCount,
            });
        }

        this.renderer.render();
    }
}

// --- Boot ---
const app = new ParticlePhysicsApp();
app.init().catch(err => console.error('ParticlePhysics init failed:', err));
