/**
 * Velocity Verlet integrator for test particle on S^3.
 *
 * PHYSICS:
 * - Particle lives on the unit 3-sphere (4D unit vector)
 * - Velocity is tangent to S^3 (perpendicular to position)
 * - Acceleration from STG geodesic: a = nabla(C) / C^3  (projected to tangent plane)
 * - Force magnitude uses GRAPH distance (topological, discrete)
 * - Force direction uses S^3 geometry (continuous, smooth)
 *
 * This is the correct hybrid for a discrete-topology theory:
 * topology is discrete (600-cell graph), geometry is continuous (S^3).
 */

const MAX_TRAIL = 800;

export class S3Particle {
    constructor() {
        this.pos = [1, 0, 0, 0];   // 4D unit vector on S^3
        this.vel = [0, 0, 0, 0];   // 4D tangent vector (vel . pos = 0)
        this.acc = [0, 0, 0, 0];   // 4D tangent acceleration
        this.dt = 0.005;
        this.active = false;
        this.nearestVertex = -1;

        // Trail: circular buffer of 4D positions
        this.trail = new Float32Array(MAX_TRAIL * 4);
        this.trailHead = 0;
        this.trailLength = 0;
    }

    /**
     * Initialize particle on S^3.
     * @param {number[]} pos4D - unit 4D vector (point on S^3)
     * @param {number[]} vel4D - tangent 4D vector (vel . pos = 0)
     */
    init(pos4D, vel4D) {
        this.pos = [...pos4D];
        this._normalize();
        this.vel = this._projectToTangent(vel4D);
        this.acc = [0, 0, 0, 0];
        this.active = true;
        this.trailHead = 0;
        this.trailLength = 0;
    }

    /**
     * Normalize position to unit sphere.
     */
    _normalize() {
        const n = Math.sqrt(
            this.pos[0] ** 2 + this.pos[1] ** 2 +
            this.pos[2] ** 2 + this.pos[3] ** 2
        );
        if (n > 1e-12) {
            this.pos[0] /= n; this.pos[1] /= n;
            this.pos[2] /= n; this.pos[3] /= n;
        }
    }

    /**
     * Project a 4D vector onto the tangent plane at current position.
     * T(v) = v - (v . pos) * pos
     */
    _projectToTangent(v) {
        const dot = v[0] * this.pos[0] + v[1] * this.pos[1] +
                    v[2] * this.pos[2] + v[3] * this.pos[3];
        return [
            v[0] - dot * this.pos[0],
            v[1] - dot * this.pos[1],
            v[2] - dot * this.pos[2],
            v[3] - dot * this.pos[3],
        ];
    }

    /**
     * Velocity Verlet step on S^3.
     *
     * 1. pos += vel*dt + 0.5*acc*dt^2
     * 2. Renormalize to S^3
     * 3. Compute new acceleration
     * 4. vel += 0.5*(acc_old + acc_new)*dt
     * 5. Re-orthogonalize velocity
     *
     * @param {function} getAcceleration - (pos4D) => [ax, ay, az, aw] tangent acceleration
     */
    step(getAcceleration) {
        if (!this.active) return;

        const dt = this.dt;

        // 1. Position update (Euler step in embedding space)
        for (let i = 0; i < 4; i++) {
            this.pos[i] += this.vel[i] * dt + 0.5 * this.acc[i] * dt * dt;
        }

        // 2. Project back to S^3
        this._normalize();

        // 3. Store old acceleration, compute new
        const oldAcc = [...this.acc];
        this.acc = getAcceleration(this.pos);

        // 4. Velocity update
        for (let i = 0; i < 4; i++) {
            this.vel[i] += 0.5 * (oldAcc[i] + this.acc[i]) * dt;
        }

        // 5. Re-orthogonalize velocity to tangent plane
        this.vel = this._projectToTangent(this.vel);

        // Record trail (4D position)
        const toff = this.trailHead * 4;
        this.trail[toff] = this.pos[0];
        this.trail[toff + 1] = this.pos[1];
        this.trail[toff + 2] = this.pos[2];
        this.trail[toff + 3] = this.pos[3];
        this.trailHead = (this.trailHead + 1) % MAX_TRAIL;
        if (this.trailLength < MAX_TRAIL) this.trailLength++;
    }

    /**
     * Get speed (magnitude of velocity tangent vector).
     */
    getSpeed() {
        return Math.sqrt(
            this.vel[0] ** 2 + this.vel[1] ** 2 +
            this.vel[2] ** 2 + this.vel[3] ** 2
        );
    }

    /**
     * Get trail as array of 4D positions (oldest first).
     * @returns {Array<number[]>}
     */
    getTrail4D() {
        const out = [];
        for (let i = 0; i < this.trailLength; i++) {
            const idx = ((this.trailHead - this.trailLength + i + MAX_TRAIL) % MAX_TRAIL) * 4;
            out.push([
                this.trail[idx], this.trail[idx + 1],
                this.trail[idx + 2], this.trail[idx + 3],
            ]);
        }
        return out;
    }
}
