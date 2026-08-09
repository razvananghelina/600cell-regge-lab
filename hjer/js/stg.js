/**
 * Surface-Tension Gravity (STG) module.
 *
 * PHYSICS (DERIVED):
 * - C(d) = 1 + M/d  where d = graph distance, M = mass in Planck units
 * - On the graph, edge_length = l_Planck = 1 (natural units)
 * - So alpha_G = M (not a free parameter!)
 * - STG geodesic acceleration: a = nabla(C) / C^3
 *
 * IMPLEMENTATION:
 * - Force magnitude: M / (d^2 * C^3) -- from graph distance (TOPOLOGICAL)
 * - Force direction: tangent vector on S^3 toward mass (GEOMETRIC)
 * - This is the correct hybrid: discrete topology + continuous geometry
 */

export class STG {
    constructor(polytope) {
        this.polytope = polytope;
        this.masses = [];           // [{vertex: i, mass: m}, ...]
        this.compression = null;    // Float32Array: C[i] for each vertex
        this.massDists = [];        // BFS distance arrays, one per mass
    }

    /**
     * Place a mass at a vertex.
     * Mass is in Planck units (M=1 = one Planck mass).
     * @param {number} vertexIndex
     * @param {number} mass - mass in Planck units (default 1.0)
     */
    placeMass(vertexIndex, mass = 1.0) {
        const existing = this.masses.find(m => m.vertex === vertexIndex);
        if (existing) {
            existing.mass = mass;
        } else {
            this.masses.push({ vertex: vertexIndex, mass });
        }
        this._recompute();
    }

    removeMass(vertexIndex) {
        this.masses = this.masses.filter(m => m.vertex !== vertexIndex);
        this._recompute();
    }

    clearMasses() {
        this.masses = [];
        this.compression = null;
        this.massDists = [];
    }

    /**
     * Recompute compression field.
     * C[i] = 1 + sum_k( M_k / dist(i, k) )
     * No free parameters -- alpha_G = M in Planck units.
     */
    _recompute() {
        const n = this.polytope.numVertices;
        this.compression = new Float32Array(n).fill(1.0);
        this.massDists = [];

        for (const { vertex, mass } of this.masses) {
            const dist = this.polytope.bfs(vertex);
            this.massDists.push(dist);
            for (let i = 0; i < n; i++) {
                if (dist[i] > 0) {
                    this.compression[i] += mass / dist[i];
                } else if (i === vertex) {
                    // At source: use d=0.5 (half-edge, natural cutoff)
                    this.compression[i] += mass / 0.5;
                }
            }
        }
    }

    /**
     * Get C range for colormap.
     */
    getRange() {
        if (!this.compression) return { min: 1, max: 1 };
        let min = Infinity, max = -Infinity;
        for (let i = 0; i < this.compression.length; i++) {
            if (this.compression[i] < min) min = this.compression[i];
            if (this.compression[i] > max) max = this.compression[i];
        }
        return { min, max };
    }

    /**
     * Find the nearest graph vertex to a point on S^3.
     * Uses dot product (cos of angular distance on S^3).
     * @param {number[]} pos4D - unit 4D vector
     * @returns {number} vertex index
     */
    nearestVertex(pos4D) {
        let maxDot = -Infinity;
        let nearest = 0;
        const verts = this.polytope.vertices4D;

        for (let i = 0; i < verts.length; i++) {
            const v = verts[i];
            const dot = pos4D[0] * v[0] + pos4D[1] * v[1] +
                        pos4D[2] * v[2] + pos4D[3] * v[3];
            if (dot > maxDot) {
                maxDot = dot;
                nearest = i;
            }
        }
        return nearest;
    }

    /**
     * Compute STG acceleration on S^3 at a given position.
     *
     * DERIVATION:
     * STG metric: ds^2 = -c^2 dt^2 / C^2 + C^2 * d_sigma^2
     * Geodesic for slow particle: a = nabla(C) / C^3
     *
     * On the graph:
     * - |nabla(C)| at distance d from mass M: dC/dd = -M/d^2
     * - C at distance d: 1 + M/d
     * - So |a| = (M/d^2) / (1+M/d)^3 = M / (d^2 * C^3)
     * - Direction: tangent vector on S^3 from pos toward mass vertex
     *
     * @param {number[]} pos4D - particle position (unit 4D vector on S^3)
     * @returns {number[]} [ax, ay, az, aw] tangent acceleration on S^3
     */
    accelerationS3(pos4D) {
        if (this.masses.length === 0 || !this.compression) return [0, 0, 0, 0];

        // Find nearest graph vertex for graph distance lookup
        const nearest = this.nearestVertex(pos4D);

        let ax = 0, ay = 0, az = 0, aw = 0;

        for (let m = 0; m < this.masses.length; m++) {
            const { vertex: massVert, mass } = this.masses[m];
            const graphDist = this.massDists[m][nearest];
            if (graphDist <= 0) continue;  // at the mass itself

            // C at particle's nearest vertex
            const C = this.compression[nearest];

            // Force magnitude: M / (d^2 * C^3)  [DERIVED from STG geodesic]
            const C3 = C * C * C;
            const forceMag = mass / (graphDist * graphDist * C3);

            // Direction on S^3: tangent vector from pos toward mass vertex
            // For two points p, q on S^3: tangent at p toward q is:
            //   t = q - (q . p) * p,  then normalize
            const mv = this.polytope.vertices4D[massVert];
            const dot = pos4D[0] * mv[0] + pos4D[1] * mv[1] +
                        pos4D[2] * mv[2] + pos4D[3] * mv[3];

            let tx = mv[0] - dot * pos4D[0];
            let ty = mv[1] - dot * pos4D[1];
            let tz = mv[2] - dot * pos4D[2];
            let tw = mv[3] - dot * pos4D[3];

            const tLen = Math.sqrt(tx * tx + ty * ty + tz * tz + tw * tw);
            if (tLen < 1e-10) continue;  // antipodal or coincident

            // Normalize and apply force
            const invT = forceMag / tLen;
            ax += tx * invT;
            ay += ty * invT;
            az += tz * invT;
            aw += tw * invT;
        }

        return [ax, ay, az, aw];
    }
}
