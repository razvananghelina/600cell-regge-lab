/**
 * 5-state Antiferromagnetic Potts Model on 600-cell.
 *
 * Monte Carlo Metropolis simulation for Kibble-Zurek cooling.
 * Ground state = proper 5-coloring (all edges inter-coset, E = 0).
 * T = infinity: random coloring. T = 0: ground state.
 * Cooling reveals frozen topological defects (~88 in slow cool).
 *
 * Key results (exp179):
 * - T_c ~ 1 (NOT 3). Entropic effect dominates.
 * - Slow cool: ~88 frozen defects. Quench: ~95.
 * - Matter = relict of cosmological phase transition.
 * - Annihilation emits predominantly into lambda=14 mode (mult=36).
 */

export class PottsModel {
    /**
     * @param {number} nVertices
     * @param {Array<number[]>} adjacency
     * @param {Array<number[]>} edges
     */
    constructor(nVertices, adjacency, edges) {
        this.n = nVertices;
        this.adjacency = adjacency;
        this.edges = edges;
        this.q = 5; // 5-state Potts

        // State: coset assignment for each vertex
        this.state = new Uint8Array(this.n);
        this.temperature = 5.0;

        // Statistics
        this.energy = 0;
        this.defectCount = 0;
        this.sweepCount = 0;
    }

    /**
     * Initialize with random state (T = infinity).
     */
    randomize() {
        for (let i = 0; i < this.n; i++) {
            this.state[i] = Math.floor(Math.random() * this.q);
        }
        this._updateStats();
        this.sweepCount = 0;
    }

    /**
     * Initialize from a given coset assignment (ground state).
     * @param {number[]} coset
     */
    initFromCoset(coset) {
        for (let i = 0; i < this.n; i++) {
            this.state[i] = coset[i];
        }
        this._updateStats();
        this.sweepCount = 0;
    }

    /**
     * Get current state as array.
     * @returns {Uint8Array}
     */
    getState() {
        return this.state;
    }

    /**
     * One Metropolis sweep (visit all vertices once in random order).
     */
    sweep() {
        // Fisher-Yates shuffle for visit order
        const order = Array.from({ length: this.n }, (_, i) => i);
        for (let i = this.n - 1; i > 0; i--) {
            const j = Math.floor(Math.random() * (i + 1));
            [order[i], order[j]] = [order[j], order[i]];
        }

        for (const v of order) {
            this._metropolisStep(v);
        }

        this._updateStats();
        this.sweepCount++;
    }

    /**
     * Metropolis step for a single vertex.
     * AF Potts: E = -J * sum_{<ij>} (1 - delta(s_i, s_j))
     * But we use the convention where E counts SAME-coset edges (defects).
     * Lower E = better (fewer defects).
     */
    _metropolisStep(v) {
        const oldState = this.state[v];

        // Propose random new state (different from current)
        let newState;
        do {
            newState = Math.floor(Math.random() * this.q);
        } while (newState === oldState);

        // Count same-coset neighbors for old and new states
        let oldSame = 0, newSame = 0;
        for (const nb of this.adjacency[v]) {
            if (this.state[nb] === oldState) oldSame++;
            if (this.state[nb] === newState) newSame++;
        }

        // Delta E = newSame - oldSame (want to minimize same-coset edges)
        const dE = newSame - oldSame;

        // Accept if dE <= 0, or with probability exp(-dE/T)
        if (dE <= 0 || (this.temperature > 1e-10 && Math.random() < Math.exp(-dE / this.temperature))) {
            this.state[v] = newState;
        }
    }

    /**
     * Update energy and defect count.
     */
    _updateStats() {
        let energy = 0;
        for (const [i, j] of this.edges) {
            if (this.state[i] === this.state[j]) {
                energy++;
            }
        }
        this.energy = energy;

        // Defect = vertex with at least one same-coset neighbor
        let defects = 0;
        for (let i = 0; i < this.n; i++) {
            let hasSame = false;
            for (const nb of this.adjacency[i]) {
                if (this.state[nb] === this.state[i]) {
                    hasSame = true;
                    break;
                }
            }
            if (hasSame) defects++;
        }
        this.defectCount = defects;
    }

    /**
     * Run multiple sweeps.
     * @param {number} numSweeps
     */
    runSweeps(numSweeps) {
        for (let i = 0; i < numSweeps; i++) {
            this.sweep();
        }
    }

    /**
     * Get cooling schedule temperatures.
     * Linear from T_start to T_end over numSteps.
     * @param {number} tStart
     * @param {number} tEnd
     * @param {number} numSteps
     * @returns {number[]}
     */
    static coolingSchedule(tStart = 5.0, tEnd = 0.01, numSteps = 100) {
        const temps = [];
        for (let i = 0; i < numSteps; i++) {
            const t = tStart * (1 - i / (numSteps - 1)) + tEnd * (i / (numSteps - 1));
            temps.push(Math.max(t, tEnd));
        }
        return temps;
    }
}
