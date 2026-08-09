/**
 * Three.js renderer: InstancedMesh for vertices, LineSegments for edges.
 * Uses Three.js r160+ from CDN (ES module).
 */

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const VERTEX_RADIUS = 0.04;
const VERTEX_SEGMENTS = 8;
const EDGE_COLOR_DEFAULT = 0x4488cc;
const PARTICLE_COLOR = 0x00ff88;
const TRAIL_COLOR = 0x00ff88;

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

        // Renderer
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

        // Geometry objects (initialized later)
        this.vertexMesh = null;      // InstancedMesh
        this.edgeLines = null;       // LineSegments
        this.particleMesh = null;    // Mesh (sphere)
        this.trailLine = null;       // Line

        // Working arrays
        this._dummy = new THREE.Object3D();
        this._color = new THREE.Color();

        // Raycaster for click picking
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
     * Initialize vertex instances and edge lines for a polytope.
     * @param {number} numVertices
     * @param {number} numEdges
     */
    initPolytope(numVertices, numEdges) {
        // Clean up old
        if (this.vertexMesh) this.scene.remove(this.vertexMesh);
        if (this.edgeLines) this.scene.remove(this.edgeLines);

        // Vertex InstancedMesh
        const sphereGeo = new THREE.SphereGeometry(VERTEX_RADIUS, VERTEX_SEGMENTS, VERTEX_SEGMENTS);
        const sphereMat = new THREE.MeshPhongMaterial({ color: 0xffcc00 });
        this.vertexMesh = new THREE.InstancedMesh(sphereGeo, sphereMat, numVertices);
        this.vertexMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
        this.scene.add(this.vertexMesh);

        // Edge LineSegments (allocate max buffer)
        const lineGeo = new THREE.BufferGeometry();
        const linePositions = new Float32Array(numEdges * 2 * 3);
        lineGeo.setAttribute('position', new THREE.BufferAttribute(linePositions, 3));
        lineGeo.setDrawRange(0, numEdges * 2);
        const lineMat = new THREE.LineBasicMaterial({
            color: EDGE_COLOR_DEFAULT,
            transparent: true,
            opacity: 0.35,
        });
        this.edgeLines = new THREE.LineSegments(lineGeo, lineMat);
        this.scene.add(this.edgeLines);

        // Particle (test mass)
        if (!this.particleMesh) {
            const pGeo = new THREE.SphereGeometry(0.08, 12, 12);
            const pMat = new THREE.MeshPhongMaterial({
                color: PARTICLE_COLOR,
                emissive: PARTICLE_COLOR,
                emissiveIntensity: 0.5,
            });
            this.particleMesh = new THREE.Mesh(pGeo, pMat);
            this.particleMesh.visible = false;
            this.scene.add(this.particleMesh);
        }

        // Trail line
        if (!this.trailLine) {
            const tGeo = new THREE.BufferGeometry();
            const tPos = new Float32Array(800 * 3);
            tGeo.setAttribute('position', new THREE.BufferAttribute(tPos, 3));
            tGeo.setDrawRange(0, 0);
            const tMat = new THREE.LineBasicMaterial({
                color: TRAIL_COLOR,
                transparent: true,
                opacity: 0.6,
            });
            this.trailLine = new THREE.Line(tGeo, tMat);
            this.trailLine.visible = false;
            this.scene.add(this.trailLine);
        }
    }

    /**
     * Update vertex positions and colors.
     * @param {Float32Array} positions3D - flat [x0,y0,z0, ...]
     * @param {Uint8Array} visible - visibility mask
     * @param {number} count
     * @param {Array<[number,number,number]>|null} colors - per-vertex RGB or null for default
     */
    updateVertices(positions3D, visible, count, colors) {
        if (!this.vertexMesh) return;

        for (let i = 0; i < count; i++) {
            const off = i * 3;
            if (visible[i]) {
                this._dummy.position.set(positions3D[off], positions3D[off + 1], positions3D[off + 2]);
                this._dummy.scale.set(1, 1, 1);
            } else {
                this._dummy.position.set(0, 0, 0);
                this._dummy.scale.set(0, 0, 0);
            }
            this._dummy.updateMatrix();
            this.vertexMesh.setMatrixAt(i, this._dummy.matrix);

            if (colors) {
                this._color.setRGB(colors[i][0], colors[i][1], colors[i][2]);
            } else {
                this._color.setRGB(1.0, 0.84, 0.0); // gold
            }
            this.vertexMesh.setColorAt(i, this._color);
        }

        this.vertexMesh.instanceMatrix.needsUpdate = true;
        if (this.vertexMesh.instanceColor) this.vertexMesh.instanceColor.needsUpdate = true;
    }

    /**
     * Update edge line positions.
     * @param {Array<number[]>} edges - [[i,j], ...]
     * @param {Float32Array} positions3D - flat vertex positions
     * @param {Uint8Array} visible
     */
    updateEdges(edges, positions3D, visible) {
        if (!this.edgeLines) return;

        const posAttr = this.edgeLines.geometry.getAttribute('position');
        const arr = posAttr.array;
        let idx = 0;
        let drawCount = 0;

        for (const [i, j] of edges) {
            if (!visible[i] || !visible[j]) continue;

            const oi = i * 3, oj = j * 3;
            arr[idx++] = positions3D[oi];
            arr[idx++] = positions3D[oi + 1];
            arr[idx++] = positions3D[oi + 2];
            arr[idx++] = positions3D[oj];
            arr[idx++] = positions3D[oj + 1];
            arr[idx++] = positions3D[oj + 2];
            drawCount += 2;
        }

        posAttr.needsUpdate = true;
        this.edgeLines.geometry.setDrawRange(0, drawCount);
    }

    /**
     * Update test particle position and trail.
     * @param {VerletParticle} particle
     */
    updateParticle(particle) {
        if (!particle || !particle.active) {
            if (this.particleMesh) this.particleMesh.visible = false;
            if (this.trailLine) this.trailLine.visible = false;
            return;
        }

        // Particle sphere
        this.particleMesh.position.set(particle.pos[0], particle.pos[1], particle.pos[2]);
        this.particleMesh.visible = true;

        // Trail
        const trailPos = particle.getTrailPositions();
        const tAttr = this.trailLine.geometry.getAttribute('position');
        tAttr.array.set(trailPos);
        tAttr.needsUpdate = true;
        this.trailLine.geometry.setDrawRange(0, particle.trailLength);
        this.trailLine.visible = particle.trailLength > 1;
    }

    /**
     * Raycast click to find nearest vertex.
     * @param {MouseEvent} event
     * @param {Float32Array} positions3D
     * @param {number} count
     * @returns {number} vertex index or -1
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

    /**
     * Render one frame.
     */
    render() {
        this.controls.update();
        this.webgl.render(this.scene, this.camera);
    }

    /**
     * Dispose resources.
     */
    dispose() {
        window.removeEventListener('resize', this._onResize);
        this.webgl.dispose();
    }
}
