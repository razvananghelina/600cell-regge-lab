/**
 * Three.js renderer for Particle Physics Visualizer.
 *
 * Extended from HJER renderer with:
 * - Per-edge coloring (gauge types)
 * - Per-vertex sizing by type (A/B/C)
 * - Soliton glow effect
 * - Multiple edge layers (one per gauge type for toggle visibility)
 */

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const BASE_VERTEX_RADIUS = 0.04;
const VERTEX_SEGMENTS = 8;

export class Renderer {
    constructor(container) {
        this.container = container;

        // Scene
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x0a0a1a);

        // Camera
        const aspect = container.clientWidth / container.clientHeight;
        this.camera = new THREE.PerspectiveCamera(60, aspect, 0.1, 100);
        this.camera.position.set(3, 2, 4);

        // WebGL Renderer
        this.webgl = new THREE.WebGLRenderer({ antialias: true });
        this.webgl.setSize(container.clientWidth, container.clientHeight);
        this.webgl.setPixelRatio(window.devicePixelRatio);
        container.appendChild(this.webgl.domElement);

        // Controls
        this.controls = new OrbitControls(this.camera, this.webgl.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.05;

        // Lights
        const ambient = new THREE.AmbientLight(0x404060, 1.0);
        this.scene.add(ambient);
        const dir = new THREE.DirectionalLight(0xffffff, 1.0);
        dir.position.set(5, 5, 5);
        this.scene.add(dir);

        // Geometry objects
        this.vertexMesh = null;       // InstancedMesh for vertices
        this.edgeLayers = {};         // { 'U1': LineSegments, 'SU2': ..., 'SU3': ... }
        this.brokenEdgeLines = null;  // LineSegments for broken (same-coset) edges

        // Working
        this._dummy = new THREE.Object3D();
        this._color = new THREE.Color();

        // Raycaster
        this.raycaster = new THREE.Raycaster();
        this._mouse = new THREE.Vector2();

        // Resize handler
        this._onResize = () => {
            const w = container.clientWidth;
            const h = container.clientHeight;
            this.camera.aspect = w / h;
            this.camera.updateProjectionMatrix();
            this.webgl.setSize(w, h);
        };
        window.addEventListener('resize', this._onResize);
    }

    /**
     * Initialize vertex instances and edge layers.
     * @param {number} numVertices
     * @param {number} numEdges
     * @param {string[]} edgeTypes - per-edge gauge type
     */
    initPolytope(numVertices, numEdges, edgeTypes) {
        // Clean up old
        if (this.vertexMesh) this.scene.remove(this.vertexMesh);
        for (const key in this.edgeLayers) {
            this.scene.remove(this.edgeLayers[key]);
        }
        if (this.brokenEdgeLines) this.scene.remove(this.brokenEdgeLines);

        // Vertex InstancedMesh
        const sphereGeo = new THREE.SphereGeometry(BASE_VERTEX_RADIUS, VERTEX_SEGMENTS, VERTEX_SEGMENTS);
        const sphereMat = new THREE.MeshPhongMaterial({ color: 0xffffff });
        this.vertexMesh = new THREE.InstancedMesh(sphereGeo, sphereMat, numVertices);
        this.vertexMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
        this.scene.add(this.vertexMesh);

        // Count edges per gauge type
        const typeCounts = { U1: 0, SU2: 0, SU3: 0 };
        if (edgeTypes) {
            for (const t of edgeTypes) typeCounts[t]++;
        }

        // Create separate LineSegments per gauge type for toggling
        const gaugeColors = {
            U1:  0xffe633,
            SU2: 0x33ccff,
            SU3: 0xff4d33,
        };

        this.edgeLayers = {};
        for (const type of ['U1', 'SU2', 'SU3']) {
            const count = typeCounts[type] || numEdges;
            const geo = new THREE.BufferGeometry();
            const positions = new Float32Array(count * 2 * 3);
            geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
            geo.setDrawRange(0, 0);

            const mat = new THREE.LineBasicMaterial({
                color: gaugeColors[type],
                transparent: true,
                opacity: type === 'U1' ? 0.7 : (type === 'SU2' ? 0.45 : 0.35),
            });

            this.edgeLayers[type] = new THREE.LineSegments(geo, mat);
            this.edgeLayers[type].userData = { gaugeType: type };
            this.scene.add(this.edgeLayers[type]);
        }

        // Broken edges layer (white, shown when solitons exist)
        const brokenGeo = new THREE.BufferGeometry();
        const brokenPos = new Float32Array(100 * 2 * 3); // max ~100 broken edges
        brokenGeo.setAttribute('position', new THREE.BufferAttribute(brokenPos, 3));
        brokenGeo.setDrawRange(0, 0);
        const brokenMat = new THREE.LineBasicMaterial({
            color: 0xffffff,
            transparent: true,
            opacity: 0.9,
            linewidth: 2,
        });
        this.brokenEdgeLines = new THREE.LineSegments(brokenGeo, brokenMat);
        this.brokenEdgeLines.visible = false;
        this.scene.add(this.brokenEdgeLines);
    }

    /**
     * Update vertex positions, colors, and sizes.
     * @param {Float32Array} positions3D
     * @param {Uint8Array} visible
     * @param {number} count
     * @param {Array<[number,number,number]>} colors - per-vertex RGB
     * @param {number[]|null} scales - per-vertex scale factors (from vertex type)
     * @param {Set<number>|null} glowSet - set of vertex indices that should glow
     * @param {number} glowPhase - animation phase for glow
     */
    updateVertices(positions3D, visible, count, colors, scales, glowSet, glowPhase) {
        if (!this.vertexMesh) return;

        for (let i = 0; i < count; i++) {
            const off = i * 3;
            if (visible[i]) {
                this._dummy.position.set(positions3D[off], positions3D[off + 1], positions3D[off + 2]);
                const s = scales ? scales[i] : 1.0;
                // Glow: make soliton vertices larger and pulsing
                let finalScale = s;
                if (glowSet && glowSet.has(i)) {
                    finalScale = s * (1.5 + 0.5 * Math.sin(glowPhase));
                }
                this._dummy.scale.set(finalScale, finalScale, finalScale);
            } else {
                this._dummy.position.set(0, 0, 0);
                this._dummy.scale.set(0, 0, 0);
            }
            this._dummy.updateMatrix();
            this.vertexMesh.setMatrixAt(i, this._dummy.matrix);

            if (colors) {
                this._color.setRGB(colors[i][0], colors[i][1], colors[i][2]);
            } else {
                this._color.setRGB(1.0, 0.84, 0.0);
            }
            this.vertexMesh.setColorAt(i, this._color);
        }

        this.vertexMesh.instanceMatrix.needsUpdate = true;
        if (this.vertexMesh.instanceColor) this.vertexMesh.instanceColor.needsUpdate = true;
    }

    /**
     * Update edges per gauge type.
     * @param {Array<number[]>} edges
     * @param {string[]} edgeTypes
     * @param {Float32Array} positions3D
     * @param {Uint8Array} visible
     */
    updateEdges(edges, edgeTypes, positions3D, visible) {
        // Build per-type edge lists
        const byType = { U1: [], SU2: [], SU3: [] };

        for (let e = 0; e < edges.length; e++) {
            const [i, j] = edges[e];
            if (!visible[i] || !visible[j]) continue;
            const type = edgeTypes[e];
            byType[type].push(e);
        }

        for (const type of ['U1', 'SU2', 'SU3']) {
            const layer = this.edgeLayers[type];
            if (!layer) continue;

            const posAttr = layer.geometry.getAttribute('position');
            const arr = posAttr.array;
            let idx = 0;
            let drawCount = 0;

            for (const e of byType[type]) {
                const [i, j] = edges[e];
                const oi = i * 3, oj = j * 3;

                if (idx + 6 > arr.length) break;

                arr[idx++] = positions3D[oi];
                arr[idx++] = positions3D[oi + 1];
                arr[idx++] = positions3D[oi + 2];
                arr[idx++] = positions3D[oj];
                arr[idx++] = positions3D[oj + 1];
                arr[idx++] = positions3D[oj + 2];
                drawCount += 2;
            }

            posAttr.needsUpdate = true;
            layer.geometry.setDrawRange(0, drawCount);
        }
    }

    /**
     * Update broken edge display.
     * @param {number[]} brokenEdgeIndices - indices into edges array
     * @param {Array<number[]>} edges
     * @param {Float32Array} positions3D
     * @param {Uint8Array} visible
     */
    updateBrokenEdges(brokenEdgeIndices, edges, positions3D, visible) {
        if (!this.brokenEdgeLines) return;

        if (!brokenEdgeIndices || brokenEdgeIndices.length === 0) {
            this.brokenEdgeLines.visible = false;
            return;
        }

        const posAttr = this.brokenEdgeLines.geometry.getAttribute('position');
        const arr = posAttr.array;
        let idx = 0;
        let drawCount = 0;

        for (const e of brokenEdgeIndices) {
            const [i, j] = edges[e];
            if (!visible[i] || !visible[j]) continue;
            const oi = i * 3, oj = j * 3;

            if (idx + 6 > arr.length) break;

            arr[idx++] = positions3D[oi];
            arr[idx++] = positions3D[oi + 1];
            arr[idx++] = positions3D[oi + 2];
            arr[idx++] = positions3D[oj];
            arr[idx++] = positions3D[oj + 1];
            arr[idx++] = positions3D[oj + 2];
            drawCount += 2;
        }

        posAttr.needsUpdate = true;
        this.brokenEdgeLines.geometry.setDrawRange(0, drawCount);
        this.brokenEdgeLines.visible = drawCount > 0;
    }

    /**
     * Toggle visibility of a gauge edge layer.
     * @param {string} gaugeType - 'U1', 'SU2', 'SU3'
     * @param {boolean} visible
     */
    setEdgeLayerVisible(gaugeType, visible) {
        if (this.edgeLayers[gaugeType]) {
            this.edgeLayers[gaugeType].visible = visible;
        }
    }

    /**
     * Raycast click to find nearest vertex.
     */
    pickVertex(event, positions3D, count) {
        const rect = this.container.getBoundingClientRect();
        this._mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
        this._mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;

        this.raycaster.setFromCamera(this._mouse, this.camera);
        const results = this.raycaster.intersectObject(this.vertexMesh);

        if (results.length > 0) {
            return results[0].instanceId;
        }
        return -1;
    }

    render() {
        this.controls.update();
        this.webgl.render(this.scene, this.camera);
    }

    dispose() {
        window.removeEventListener('resize', this._onResize);
        this.webgl.dispose();
    }
}
