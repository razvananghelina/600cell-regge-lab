// Physical and geometric constants for 600-cell theory

export const PHI = (1 + Math.sqrt(5)) / 2;           // 1.6180339887...
export const ALPHA = 7.2973525693e-3;                  // fine structure constant
export const ALPHA_INV = 1 / ALPHA;                    // ~137.036
export const ALPHA_S = 0.1179;                         // strong coupling at M_Z
export const SIN2_THETA_W = 0.23122;                   // Weinberg angle
export const PI = Math.PI;

// Derived geometric constants
export const PHI_SQ = PHI * PHI;                       // phi^2 = phi + 1
export const PHI_INV = 1 / PHI;                        // 1/phi = phi - 1
export const ICOSA_ALPHA = 20 * PHI ** 4;              // ~137.508 (600-cell formula)
