/**
 * Stereographic projection 4D -> 3D with anamorphic clamp.
 * Projects from the north pole (0,0,0,1).
 */

const MAX_RADIUS = 10;

/**
 * Project a single rotated 4D point to 3D.
 * @param {number} x
 * @param {number} y
 * @param {number} z
 * @param {number} w
 * @returns {[number, number, number]|null} 3D point or null if at pole
 */
export function stereographic(x, y, z, w) {
    if (Math.abs(w - 1) < 1e-6) return null;
    const scale = 1 / (1 - w);
    const px = x * scale;
    const py = y * scale;
    const pz = z * scale;

    // Clamp to max radius
    const r = Math.sqrt(px * px + py * py + pz * pz);
    if (r > MAX_RADIUS) {
        const f = MAX_RADIUS / r;
        return [px * f, py * f, pz * f];
    }
    return [px, py, pz];
}

/**
 * Project all rotated 4D vertices to 3D positions.
 * @param {Float32Array} rotated4D - flat [x0,y0,z0,w0, x1,...] from Rotation4D.transformAll
 * @param {Float32Array} out3D - flat [x0,y0,z0, x1,y1,z1, ...]
 * @param {number} count - number of vertices
 * @returns {Uint8Array} visibility mask (1=visible, 0=at pole)
 */
export function projectAll(rotated4D, out3D, count) {
    const visible = new Uint8Array(count);

    for (let i = 0; i < count; i++) {
        const off4 = i * 4;
        const off3 = i * 3;
        const w = rotated4D[off4 + 3];

        if (Math.abs(w - 1) < 1e-6) {
            out3D[off3] = 0;
            out3D[off3 + 1] = 0;
            out3D[off3 + 2] = 0;
            visible[i] = 0;
            continue;
        }

        const scale = 1 / (1 - w);
        let px = rotated4D[off4] * scale;
        let py = rotated4D[off4 + 1] * scale;
        let pz = rotated4D[off4 + 2] * scale;

        const r = Math.sqrt(px * px + py * py + pz * pz);
        if (r > MAX_RADIUS) {
            const f = MAX_RADIUS / r;
            px *= f;
            py *= f;
            pz *= f;
        }

        out3D[off3] = px;
        out3D[off3 + 1] = py;
        out3D[off3 + 2] = pz;
        visible[i] = 1;
    }

    return visible;
}
