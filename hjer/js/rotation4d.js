/**
 * 4D rotation matrices for the 6 rotation planes.
 * Planes: XY, XZ, XW, YZ, YW, ZW
 */

const PLANES = ['xy', 'xz', 'xw', 'yz', 'yw', 'zw'];

/**
 * Create a 4x4 rotation matrix for a given plane and angle.
 * Returns a flat Float64Array(16) in row-major order.
 */
function makeRotation4D(angle, plane) {
    // Start with identity
    const R = new Float64Array(16);
    R[0] = 1; R[5] = 1; R[10] = 1; R[15] = 1;

    const c = Math.cos(angle);
    const s = Math.sin(angle);

    // Indices: R[row*4 + col]
    switch (plane) {
        case 'xy':  // rotate in XY plane (indices 0,1)
            R[0] = c;  R[1] = -s;
            R[4] = s;  R[5] = c;
            break;
        case 'xz':  // rotate in XZ plane (indices 0,2)
            R[0] = c;  R[2] = -s;
            R[8] = s;  R[10] = c;
            break;
        case 'xw':  // rotate in XW plane (indices 0,3)
            R[0] = c;  R[3] = -s;
            R[12] = s; R[15] = c;
            break;
        case 'yz':  // rotate in YZ plane (indices 1,2)
            R[5] = c;  R[6] = -s;
            R[9] = s;  R[10] = c;
            break;
        case 'yw':  // rotate in YW plane (indices 1,3)
            R[5] = c;  R[7] = -s;
            R[13] = s; R[15] = c;
            break;
        case 'zw':  // rotate in ZW plane (indices 2,3)
            R[10] = c; R[11] = -s;
            R[14] = s; R[15] = c;
            break;
    }
    return R;
}

/**
 * Multiply two 4x4 matrices (row-major).
 */
function mul4x4(A, B) {
    const C = new Float64Array(16);
    for (let i = 0; i < 4; i++) {
        for (let j = 0; j < 4; j++) {
            let sum = 0;
            for (let k = 0; k < 4; k++) {
                sum += A[i * 4 + k] * B[k * 4 + j];
            }
            C[i * 4 + j] = sum;
        }
    }
    return C;
}

/**
 * Apply a 4x4 matrix to a 4D vector [x,y,z,w].
 */
function transform4D(M, v) {
    return [
        M[0] * v[0] + M[1] * v[1] + M[2] * v[2] + M[3] * v[3],
        M[4] * v[0] + M[5] * v[1] + M[6] * v[2] + M[7] * v[3],
        M[8] * v[0] + M[9] * v[1] + M[10] * v[2] + M[11] * v[3],
        M[12] * v[0] + M[13] * v[1] + M[14] * v[2] + M[15] * v[3],
    ];
}

export class Rotation4D {
    constructor() {
        // 6 rotation angles (radians)
        this.angles = { xy: 0, xz: 0, xw: 0, yz: 0, yw: 0, zw: 0 };
        this._matrix = null;
        this._dirty = true;
    }

    /**
     * Set angle for a specific plane.
     */
    setAngle(plane, radians) {
        this.angles[plane] = radians;
        this._dirty = true;
    }

    /**
     * Add delta to a specific plane angle.
     */
    addAngle(plane, delta) {
        this.angles[plane] += delta;
        this._dirty = true;
    }

    /**
     * Get the composed 4x4 rotation matrix.
     */
    getMatrix() {
        if (this._dirty) {
            // Compose all 6 rotations: XY * XZ * XW * YZ * YW * ZW
            let M = makeRotation4D(this.angles.xy, 'xy');
            M = mul4x4(M, makeRotation4D(this.angles.xz, 'xz'));
            M = mul4x4(M, makeRotation4D(this.angles.xw, 'xw'));
            M = mul4x4(M, makeRotation4D(this.angles.yz, 'yz'));
            M = mul4x4(M, makeRotation4D(this.angles.yw, 'yw'));
            M = mul4x4(M, makeRotation4D(this.angles.zw, 'zw'));
            this._matrix = M;
            this._dirty = false;
        }
        return this._matrix;
    }

    /**
     * Transform a single 4D vertex.
     */
    transformVertex(v) {
        return transform4D(this.getMatrix(), v);
    }

    /**
     * Transform all vertices in-place into output array.
     * @param {Array<number[]>} src - source 4D vertices
     * @param {Float32Array} dst - destination flat array [x0,y0,z0,w0, x1,...]
     */
    transformAll(src, dst) {
        const M = this.getMatrix();
        for (let i = 0; i < src.length; i++) {
            const v = src[i];
            const off = i * 4;
            dst[off] = M[0] * v[0] + M[1] * v[1] + M[2] * v[2] + M[3] * v[3];
            dst[off + 1] = M[4] * v[0] + M[5] * v[1] + M[6] * v[2] + M[7] * v[3];
            dst[off + 2] = M[8] * v[0] + M[9] * v[1] + M[10] * v[2] + M[11] * v[3];
            dst[off + 3] = M[12] * v[0] + M[13] * v[1] + M[14] * v[2] + M[15] * v[3];
        }
    }
}

export { PLANES };
