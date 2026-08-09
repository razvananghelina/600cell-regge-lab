/**
 * Scalar to RGB colormap for STG compression visualization.
 * Maps C(r) values to a blue (uncompressed) -> red (compressed) gradient.
 */

/**
 * Map a scalar value to an RGB color.
 * @param {number} value - the scalar (e.g. C(r))
 * @param {number} min - minimum expected value
 * @param {number} max - maximum expected value
 * @returns {[number, number, number]} RGB in [0,1]
 */
export function scalarToRGB(value, min, max) {
    // Normalize to [0, 1]
    let t = (max > min) ? (value - min) / (max - min) : 0;
    t = Math.max(0, Math.min(1, t));

    // Blue (cold, uncompressed) -> White (mid) -> Red (hot, compressed)
    let r, g, b;
    if (t < 0.5) {
        // Blue to white
        const s = t * 2;
        r = s;
        g = s;
        b = 1;
    } else {
        // White to red
        const s = (t - 0.5) * 2;
        r = 1;
        g = 1 - s;
        b = 1 - s;
    }

    return [r, g, b];
}

/**
 * Map graph distance to RGB.
 * Closer = warmer (red), farther = cooler (blue).
 * @param {number} dist - graph distance
 * @param {number} maxDist - maximum distance (graph diameter)
 * @returns {[number, number, number]} RGB in [0,1]
 */
export function distanceToRGB(dist, maxDist) {
    if (dist < 0) return [0.3, 0.3, 0.3]; // unreachable = gray
    return scalarToRGB(maxDist - dist, 0, maxDist);
}

/**
 * Default vertex color (golden).
 * @returns {[number, number, number]}
 */
export function defaultColor() {
    return [1.0, 0.84, 0.0]; // gold
}
