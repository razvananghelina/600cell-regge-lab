/**
 * TMO: Topological Monitoring Observables
 * - RCS: Retea Curbura Scalara (chi = V - E + F - C = 0 check)
 * - TDI: Topological Drift Indicator (phase drift monitoring)
 */

import { PHI } from './constants.js';

export class TMO {
    constructor() {
        this.chi = null;           // Euler characteristic
        this.phaseDrift = 0;       // accumulated phase drift
        this.driftThreshold = Math.log(PHI);  // ln(phi) ~ 0.4812
        this.driftHistory = [];
        this.rcsValid = false;
    }

    /**
     * RCS: Compute Euler characteristic chi = V - E + F - C
     * For 600-cell: V=120, E=720, F=1200, C=600 => chi = 0
     * For 120-cell: V=600, E=1200, F=720, C=120 => chi = 0
     * (These are the cells of the polytope on S^3, so chi of S^3 = 0)
     *
     * @param {object} params - {V, E, F, C} counts
     * @returns {{chi: number, valid: boolean}}
     */
    computeRCS(params) {
        const { V, E, F, C } = params;
        this.chi = V - E + F - C;
        this.rcsValid = (this.chi === 0);
        return { chi: this.chi, valid: this.rcsValid };
    }

    /**
     * RCS check for loaded polytope data.
     * Uses known face/cell counts based on polytope type.
     * @param {PolytopeGraph} polytope
     * @returns {{chi: number, valid: boolean, detail: string}}
     */
    checkPolytope(polytope) {
        let F, C;
        if (polytope.numVertices === 120) {
            // 600-cell
            F = 1200;
            C = 600;
        } else if (polytope.numVertices === 600) {
            // 120-cell
            F = 720;
            C = 120;
        } else {
            return { chi: null, valid: false, detail: 'Unknown polytope' };
        }

        const result = this.computeRCS({
            V: polytope.numVertices,
            E: polytope.numEdges,
            F,
            C,
        });

        const detail = `V=${polytope.numVertices} - E=${polytope.numEdges} + F=${F} - C=${C} = ${this.chi}`;
        return { ...result, detail };
    }

    /**
     * TDI: Update phase drift from STG compression field.
     * Monitors if compression causes topological instability.
     * Phase drift = max deviation of C from 1.
     *
     * @param {Float32Array} compression - C[i] array
     * @returns {{drift: number, threshold: number, alert: boolean}}
     */
    updateDrift(compression) {
        if (!compression) {
            return { drift: 0, threshold: this.driftThreshold, alert: false };
        }

        let maxDev = 0;
        for (let i = 0; i < compression.length; i++) {
            const dev = Math.abs(compression[i] - 1);
            if (dev > maxDev) maxDev = dev;
        }

        this.phaseDrift = maxDev;
        this.driftHistory.push(maxDev);
        if (this.driftHistory.length > 200) {
            this.driftHistory.shift();
        }

        const alert = this.phaseDrift > this.driftThreshold;
        return { drift: this.phaseDrift, threshold: this.driftThreshold, alert };
    }

    /**
     * Reset drift monitoring.
     */
    resetDrift() {
        this.phaseDrift = 0;
        this.driftHistory = [];
    }
}
