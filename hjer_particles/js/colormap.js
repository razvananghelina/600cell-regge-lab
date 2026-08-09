/**
 * Color maps for the Particle Physics Visualizer.
 *
 * Coset colors (5 distinct, high contrast):
 *   0: Red, 1: Blue, 2: Green, 3: Gold, 4: Purple
 *
 * Gauge channel colors:
 *   U(1): Yellow, SU(2): Cyan, SU(3): Red
 *
 * Temperature colormap for Potts cooling.
 */

/**
 * 5 coset colors (RGB in [0,1]).
 */
export const COSET_COLORS = [
    [0.9, 0.25, 0.2],   // 0: Red
    [0.2, 0.5, 0.95],   // 1: Blue
    [0.2, 0.8, 0.3],    // 2: Green
    [0.95, 0.75, 0.1],  // 3: Gold
    [0.7, 0.3, 0.9],    // 4: Purple
];

/**
 * Coset color names for legend.
 */
export const COSET_NAMES = ['Red', 'Blue', 'Green', 'Gold', 'Purple'];

/**
 * Gauge type colors (RGB in [0,1]).
 */
export const GAUGE_COLORS = {
    U1:  [1.0, 0.9, 0.2],   // yellow
    SU2: [0.2, 0.8, 1.0],   // cyan
    SU3: [1.0, 0.3, 0.2],   // red-orange
};

/**
 * Vertex type sizes (relative scale factor).
 * A = large (gauge/24-cell), B = medium, C = small (fermion).
 */
export const VERTEX_SIZES = {
    A: 1.8,
    B: 1.3,
    C: 1.0,
};

/**
 * Get coset color for a vertex.
 * @param {number} cosetIdx - 0 to 4
 * @returns {[number, number, number]}
 */
export function cosetColor(cosetIdx) {
    return COSET_COLORS[cosetIdx] || [0.5, 0.5, 0.5];
}

/**
 * Get gauge edge color.
 * @param {string} gaugeType - 'U1', 'SU2', 'SU3'
 * @returns {[number, number, number]}
 */
export function gaugeEdgeColor(gaugeType) {
    return GAUGE_COLORS[gaugeType] || [0.4, 0.4, 0.4];
}

/**
 * Temperature to color (for Potts display).
 * Cold (T=0): blue. Hot (T=5): red-white.
 * @param {number} T
 * @param {number} maxT
 * @returns {[number, number, number]}
 */
export function temperatureColor(T, maxT = 5.0) {
    let t = Math.max(0, Math.min(1, T / maxT));
    if (t < 0.5) {
        const s = t * 2;
        return [s, s * 0.5, 1.0];
    } else {
        const s = (t - 0.5) * 2;
        return [1.0, 0.5 * (1 - s), 1 - s];
    }
}

/**
 * Soliton highlight color (bright pulsing white-yellow).
 * @param {number} phase - animation phase (0 to 1)
 * @returns {[number, number, number]}
 */
export function solitonGlowColor(phase) {
    const p = 0.6 + 0.4 * Math.sin(phase * Math.PI * 2);
    return [1.0, 0.95 * p, 0.6 * p];
}

/**
 * Broken edge color (white).
 * @returns {[number, number, number]}
 */
export function brokenEdgeColor() {
    return [1.0, 1.0, 1.0];
}

/**
 * Default vertex color when no mode is active.
 * @returns {[number, number, number]}
 */
export function defaultColor() {
    return [1.0, 0.84, 0.0]; // gold
}
