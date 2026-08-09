/**
 * HJER Engine - Main bootstrap and animation loop.
 *
 * ARCHITECTURE:
 * - Physics runs on S^3 (the 3-sphere where 600-cell vertices live)
 * - View rotation is ONLY for visualization (does not affect physics)
 * - Particle position & trail are in 4D, projected to 3D only for rendering
 * - No free parameters in physics: alpha_G = M (Planck units)
 */

import { PolytopeGraph } from './polytope.js';
import { Rotation4D } from './rotation4d.js';
import { projectAll, stereographic } from './projection.js';
import { Renderer } from './renderer.js';
import { STG } from './stg.js';
import { S3Particle } from './verlet.js';
import { TMO } from './tmo.js';
import { scalarToRGB, defaultColor } from './colormap.js';
import { UI } from './ui.js';

class HJERApp {
    constructor() {
        this.polytope = new PolytopeGraph();
        this.rotation4d = new Rotation4D();
        this.renderer = new Renderer(document.getElementById('viewport'));
        this.stg = new STG(this.polytope);
        this.particle = new S3Particle();
        this.tmo = new TMO();
        this.ui = null;

        // Working buffers
        this.rotated4D = null;
        this.projected3D = null;
        this.visible = null;
        this.vertexColors = null;

        // State
        this.autoRotate = true;
        this.autoRotateSpeed = 0.005;
        this.placeMassMode = false;
        this.running = false;

        // FPS
        this._frameCount = 0;
        this._lastFpsTime = performance.now();
        this._fps = 0;
    }

    async init() {
        this.ui = new UI(this);
        this.renderer.container.addEventListener('click', (e) => this._onClick(e));
        await this.loadPolytope('600cell');
        this.running = true;
        this._animate();
    }

    async loadPolytope(name) {
        const url = `data/${name}.json`;
        try {
            await this.polytope.load(url);
        } catch (err) {
            console.error(`Failed to load ${url}:`, err);
            return;
        }

        const n = this.polytope.numVertices;
        this.rotated4D = new Float32Array(n * 4);
        this.projected3D = new Float32Array(n * 3);
        this.visible = new Uint8Array(n);
        this.vertexColors = new Array(n);
        for (let i = 0; i < n; i++) this.vertexColors[i] = defaultColor();

        this.renderer.initPolytope(n, this.polytope.numEdges);
        this.stg = new STG(this.polytope);
        this.particle.active = false;

        const rcs = this.tmo.checkPolytope(this.polytope);
        console.log(`TMO RCS: ${rcs.detail} => ${rcs.valid ? 'OK' : 'FAIL'}`);
    }

    _onClick(e) {
        if (!this.placeMassMode) return;

        const idx = this.renderer.pickVertex(e, this.projected3D, this.polytope.numVertices);
        if (idx >= 0) {
            const massVal = this.ui.getMassValue();
            this.stg.placeMass(idx, massVal);
            console.log(`Mass ${massVal} M_P placed at vertex ${idx}`);
            // Stay in placement mode - user clicks button again to exit
        }
    }

    /**
     * Launch particle in orbit around first mass.
     * All physics on S^3 -- no projected-space quantities.
     */
    launchParticle() {
        if (this.stg.masses.length === 0) {
            console.log('Place a mass first!');
            return;
        }

        const massVert = this.stg.masses[0].vertex;
        const massVal = this.stg.masses[0].mass;
        const dist = this.polytope.bfs(massVert);

        // Find a vertex at graph distance 3 from mass
        let startVert = -1;
        for (let i = 0; i < dist.length; i++) {
            if (dist[i] === 3) { startVert = i; break; }
        }
        if (startVert === -1) {
            for (let i = 0; i < dist.length; i++) {
                if (dist[i] === 2) { startVert = i; break; }
            }
        }
        if (startVert === -1) return;

        // Particle starts at this vertex's S^3 position
        const pos = [...this.polytope.vertices4D[startVert]];
        const massPos = this.polytope.vertices4D[massVert];

        // Compute tangent direction on S^3 from particle toward mass
        const dot_pm = pos[0] * massPos[0] + pos[1] * massPos[1] +
                       pos[2] * massPos[2] + pos[3] * massPos[3];
        const radial = [
            massPos[0] - dot_pm * pos[0],
            massPos[1] - dot_pm * pos[1],
            massPos[2] - dot_pm * pos[2],
            massPos[3] - dot_pm * pos[3],
        ];
        const radLen = Math.sqrt(radial[0] ** 2 + radial[1] ** 2 +
                                  radial[2] ** 2 + radial[3] ** 2);

        if (radLen < 1e-8) return;

        // Normalize radial direction
        for (let i = 0; i < 4; i++) radial[i] /= radLen;

        // Find a tangent vector perpendicular to radial (for tangential velocity)
        // Use Gram-Schmidt: pick an arbitrary vector, orthogonalize to pos and radial
        let arb = [0, 1, 0, 0];
        // Remove pos component
        let d1 = arb[0] * pos[0] + arb[1] * pos[1] + arb[2] * pos[2] + arb[3] * pos[3];
        for (let i = 0; i < 4; i++) arb[i] -= d1 * pos[i];
        // Remove radial component
        let d2 = arb[0] * radial[0] + arb[1] * radial[1] + arb[2] * radial[2] + arb[3] * radial[3];
        for (let i = 0; i < 4; i++) arb[i] -= d2 * radial[i];
        let arbLen = Math.sqrt(arb[0] ** 2 + arb[1] ** 2 + arb[2] ** 2 + arb[3] ** 2);

        if (arbLen < 1e-8) {
            // Try another vector
            arb = [0, 0, 1, 0];
            d1 = arb[0] * pos[0] + arb[1] * pos[1] + arb[2] * pos[2] + arb[3] * pos[3];
            for (let i = 0; i < 4; i++) arb[i] -= d1 * pos[i];
            d2 = arb[0] * radial[0] + arb[1] * radial[1] + arb[2] * radial[2] + arb[3] * radial[3];
            for (let i = 0; i < 4; i++) arb[i] -= d2 * radial[i];
            arbLen = Math.sqrt(arb[0] ** 2 + arb[1] ** 2 + arb[2] ** 2 + arb[3] ** 2);
        }

        for (let i = 0; i < 4; i++) arb[i] /= arbLen;

        // Orbital speed from STG: v_orbit = sqrt(forceMag * r_angular)
        // forceMag = M / (d^2 * C^3), r_angular = angular distance on S^3
        const graphDist = 3;
        const C = 1 + massVal / graphDist;
        const forceMag = massVal / (graphDist * graphDist * C * C * C);
        const r_angular = radLen;  // this IS the angular distance (tangent length on unit S^3)
        const speed = Math.sqrt(forceMag * r_angular);

        const vel = [
            arb[0] * speed,
            arb[1] * speed,
            arb[2] * speed,
            arb[3] * speed,
        ];

        this.particle.init(pos, vel);
        console.log(`Particle on S^3: d=${graphDist}, C=${C.toFixed(3)}, F=${forceMag.toFixed(5)}, v=${speed.toFixed(5)}`);
    }

    /**
     * Apply view rotation and project a 4D point to 3D.
     * View rotation does NOT affect physics.
     */
    _projectPoint4D(p4d) {
        const M = this.rotation4d.getMatrix();
        const r = [
            M[0] * p4d[0] + M[1] * p4d[1] + M[2] * p4d[2] + M[3] * p4d[3],
            M[4] * p4d[0] + M[5] * p4d[1] + M[6] * p4d[2] + M[7] * p4d[3],
            M[8] * p4d[0] + M[9] * p4d[1] + M[10] * p4d[2] + M[11] * p4d[3],
            M[12] * p4d[0] + M[13] * p4d[1] + M[14] * p4d[2] + M[15] * p4d[3],
        ];
        return stereographic(r[0], r[1], r[2], r[3]);
    }

    _animate() {
        if (!this.running) return;
        requestAnimationFrame(() => this._animate());

        const n = this.polytope.numVertices;
        if (n === 0) { this.renderer.render(); return; }

        // --- View rotation (visualization only) ---
        if (this.autoRotate) {
            this.rotation4d.addAngle('xw', this.autoRotateSpeed);
            this.rotation4d.addAngle('yz', this.autoRotateSpeed * 0.7);
            this.ui.syncSliders(this.rotation4d.angles);
        }

        // --- Project vertices for rendering ---
        this.rotation4d.transformAll(this.polytope.vertices4D, this.rotated4D);
        this.visible = projectAll(this.rotated4D, this.projected3D, n);

        // --- Vertex colors from C field ---
        if (this.stg.compression) {
            const { min, max } = this.stg.getRange();
            for (let i = 0; i < n; i++) {
                this.vertexColors[i] = scalarToRGB(this.stg.compression[i], min, max);
            }
        } else {
            for (let i = 0; i < n; i++) {
                this.vertexColors[i] = defaultColor();
            }
        }

        this.renderer.updateVertices(this.projected3D, this.visible, n, this.vertexColors);
        this.renderer.updateEdges(this.polytope.edges, this.projected3D, this.visible);

        // --- Physics: particle on S^3 (independent of view rotation) ---
        if (this.particle.active && this.stg.compression) {
            const stg = this.stg;
            for (let i = 0; i < 5; i++) {
                this.particle.step((pos) => stg.accelerationS3(pos));
            }

            // Project particle position through view rotation + stereographic
            const p3d = this._projectPoint4D(this.particle.pos);
            if (p3d) {
                this.renderer.particleMesh.position.set(p3d[0], p3d[1], p3d[2]);
                this.renderer.particleMesh.visible = true;
            } else {
                this.renderer.particleMesh.visible = false;
            }

            // Project trail
            const trail4D = this.particle.getTrail4D();
            const tAttr = this.renderer.trailLine.geometry.getAttribute('position');
            let trailCount = 0;
            for (let i = 0; i < trail4D.length; i++) {
                const tp = this._projectPoint4D(trail4D[i]);
                if (tp) {
                    tAttr.array[trailCount * 3] = tp[0];
                    tAttr.array[trailCount * 3 + 1] = tp[1];
                    tAttr.array[trailCount * 3 + 2] = tp[2];
                    trailCount++;
                }
            }
            tAttr.needsUpdate = true;
            this.renderer.trailLine.geometry.setDrawRange(0, trailCount);
            this.renderer.trailLine.visible = trailCount > 1;
        }

        // --- TMO ---
        const tdi = this.tmo.updateDrift(this.stg.compression);

        // --- FPS ---
        this._frameCount++;
        const now = performance.now();
        if (now - this._lastFpsTime >= 1000) {
            this._fps = this._frameCount;
            this._frameCount = 0;
            this._lastFpsTime = now;
        }

        // --- Info panel ---
        const rcs = this.tmo.checkPolytope(this.polytope);
        const { min: cMin, max: cMax } = this.stg.getRange();
        this.ui.updateInfo({
            name: this.polytope.name,
            vertices: this.polytope.numVertices,
            edges: this.polytope.numEdges,
            degree: this.polytope.degree,
            diameter: this.polytope.graphDiameter,
            fps: this._fps,
            chi: rcs.chi,
            chiValid: rcs.valid,
            massCount: this.stg.masses.length,
            cMin: cMin.toFixed(3),
            cMax: cMax.toFixed(3),
            drift: tdi.drift.toFixed(4),
            driftThreshold: tdi.threshold.toFixed(4),
            driftAlert: tdi.alert,
            particleActive: this.particle.active,
            particleSpeed: this.particle.active ? this.particle.getSpeed().toFixed(5) : '-',
            particleNearest: this.particle.active ?
                this.stg.nearestVertex(this.particle.pos) : '-',
        });

        this.renderer.render();
    }
}

const app = new HJERApp();
app.init().catch(err => console.error('HJER init failed:', err));
