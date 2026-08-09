/**
 * Soliton physics on the 600-cell.
 *
 * A soliton = vertex flipped from its ground-state coset to a different coset.
 * Key results from exp164-167, exp177-179:
 *   E_flip = 3 (topological, uniform for all vertices)
 *   E_pair = 4 (two adjacent flips)
 *   E_bag  = 6 = 2 * E_flip (spectral energy of bag)
 *   Binding = 33% (pair: 6 - 4 = 2 saved)
 *   Each eigenvalue group splits into (m-3, 2, 1) sub-levels
 *
 * Ground state: proper 5-coloring (all edges inter-coset, E = 0).
 * Flip vertex v from coset c0 to coset c1:
 *   - 12 neighbors of v: each was in some coset != c0
 *   - After flip: 3 neighbors that were in coset c1 now have same-coset edge (broken)
 *   - E_flip = 3 (always, because each vertex has exactly 3 neighbors per other coset)
 */

export class SolitonManager {
    /**
     * @param {number[]} groundCoset - original coset assignment (0-4)
     * @param {Array<number[]>} adjacency
     * @param {Array<number[]>} edges
     */
    constructor(groundCoset, adjacency, edges) {
        this.groundCoset = [...groundCoset];
        this.currentCoset = [...groundCoset];
        this.adjacency = adjacency;
        this.edges = edges;
        this.n = groundCoset.length;

        // Track which vertices are flipped
        this.flipped = new Uint8Array(this.n);
        this.solitonCount = 0;
    }

    /**
     * Flip a vertex to a new coset.
     * @param {number} vertex - vertex index
     * @param {number} targetCoset - coset to flip to (0-4)
     * @returns {{ energy: number, broken: number[], vertex: number }}
     */
    flipVertex(vertex, targetCoset) {
        const oldCoset = this.currentCoset[vertex];
        if (targetCoset === oldCoset) return null;

        this.currentCoset[vertex] = targetCoset;
        this.flipped[vertex] = 1;
        this.solitonCount++;

        // Count broken edges (same-coset edges = defects)
        const broken = [];
        for (const nb of this.adjacency[vertex]) {
            if (this.currentCoset[nb] === targetCoset) {
                broken.push(nb);
            }
        }

        return {
            energy: broken.length,
            broken,
            vertex,
            oldCoset,
            newCoset: targetCoset,
        };
    }

    /**
     * Create a soliton-antisoliton pair on adjacent vertices.
     * @param {number} v1 - first vertex
     * @param {number} v2 - second vertex (must be neighbor of v1)
     * @returns {{ energy: number, details: object[] }|null}
     */
    createPair(v1, v2) {
        // Check adjacency
        if (!this.adjacency[v1].includes(v2)) return null;

        const c1 = this.currentCoset[v1];
        const c2 = this.currentCoset[v2];

        // Flip v1 to v2's coset, v2 to v1's coset (swap)
        const r1 = this.flipVertex(v1, c2);
        const r2 = this.flipVertex(v2, c1);

        return {
            energy: r1.energy + r2.energy,
            details: [r1, r2],
        };
    }

    /**
     * Compute total energy (number of same-coset edges).
     * In ground state this is 0 (proper 5-coloring).
     * @returns {number}
     */
    computeTotalEnergy() {
        let energy = 0;
        for (const [i, j] of this.edges) {
            if (this.currentCoset[i] === this.currentCoset[j]) {
                energy++;
            }
        }
        return energy;
    }

    /**
     * Get list of broken edges (same-coset).
     * @returns {number[]} indices into edges array
     */
    getBrokenEdges() {
        const broken = [];
        for (let e = 0; e < this.edges.length; e++) {
            const [i, j] = this.edges[e];
            if (this.currentCoset[i] === this.currentCoset[j]) {
                broken.push(e);
            }
        }
        return broken;
    }

    /**
     * Get list of flipped vertex indices.
     * @returns {number[]}
     */
    getFlippedVertices() {
        const result = [];
        for (let i = 0; i < this.n; i++) {
            if (this.flipped[i]) result.push(i);
        }
        return result;
    }

    /**
     * Reset to ground state.
     */
    reset() {
        this.currentCoset = [...this.groundCoset];
        this.flipped = new Uint8Array(this.n);
        this.solitonCount = 0;
    }

    /**
     * Pick the best target coset for a flip (maximally different from neighbors).
     * For ground state: any coset != current will give E=3.
     * For perturbed state: picks coset that minimizes new broken edges.
     * @param {number} vertex
     * @returns {number} best target coset
     */
    pickTargetCoset(vertex) {
        const current = this.currentCoset[vertex];
        let bestCoset = -1;
        let bestBroken = Infinity;

        for (let c = 0; c < 5; c++) {
            if (c === current) continue;
            let broken = 0;
            for (const nb of this.adjacency[vertex]) {
                if (this.currentCoset[nb] === c) broken++;
            }
            if (broken < bestBroken) {
                bestBroken = broken;
                bestCoset = c;
            }
        }
        return bestCoset;
    }

    /**
     * Compute E_bag for a single soliton (13-vertex bag: defect + 12 neighbors).
     * E_bag = 6.0 EXACT, UNIFORM.
     * @param {number} vertex - the soliton vertex
     * @returns {number}
     */
    computeBagEnergy(vertex) {
        // Bag = vertex + its 12 neighbors
        const bag = new Set([vertex, ...this.adjacency[vertex]]);
        let energy = 0;

        for (const [i, j] of this.edges) {
            if (bag.has(i) && bag.has(j)) {
                if (this.currentCoset[i] === this.currentCoset[j]) {
                    energy++;
                }
            }
        }
        return energy;
    }
}
