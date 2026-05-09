/*
 * schneider3d.h — 3D port of Philip J. Schneider's "Algorithm for
 * Automatically Fitting Digitized Curves" (Graphics Gems I, 1990).
 *
 * Original: Point2/Vector2 ops, public-domain code.  This file ports
 * the algorithm to 3D (Point3/Vector3) by extending the inner-product,
 * squared-length and distance ops; the fitting flow is unchanged
 * (chord-length param + Newton-Raphson reparam + recursive bisection).
 *
 * This is a vendored, modified copy used only for cross-comparison
 * benchmark inside the liscio test harness.  Schneider's algorithm
 * is licensed under EULA terms (essentially public domain for
 * non-commercial; see Graphics Gems repo erich666/GraphicsGems).
 */
#ifndef SCHNEIDER3D_H
#define SCHNEIDER3D_H

#ifdef __cplusplus
extern "C" {
#endif

/* Output cubic Bezier piece (P0..P3 are 3D control points). */
typedef struct {
    double P0[3], P1[3], P2[3], P3[3];
} schneider3d_bezier_t;

/* Caller-provided sink: invoked once per emitted cubic Bezier. */
typedef void (*schneider3d_emit_t)(const schneider3d_bezier_t *b, void *user);

/* Fit cubic Bezier curve(s) to N input 3D points (XYZ flat array,
 * length 3*N).  `error` is squared distance tolerance (units²),
 * matching the original Schneider API.  The algorithm recursively
 * bisects until each piece's max squared deviation < error.
 * Each piece is delivered via the emit callback. */
void schneider3d_fit(const double *xyz_flat, int n_points,
                     double error_sq,
                     schneider3d_emit_t emit, void *user);

#ifdef __cplusplus
}
#endif
#endif /* SCHNEIDER3D_H */
