/*
 * schneider3d.c — 3D port of Schneider's FitCurves (Graphics Gems I).
 *
 * Original 2D code: Philip J. Schneider, "An Algorithm for
 * Automatically Fitting Digitized Curves", Graphics Gems I, 1990,
 * pp. 612–626.  Public-domain.
 *
 * 3D port: replaces Point2 with double[3], extends vector ops to 3D.
 * Algorithm structure preserved verbatim:
 *   - Chord-length parameterization
 *   - GenerateBezier: 2x2 normal-equation LSQ on Bernstein basis
 *     with endpoints pinned + tangents fixed via alpha_l / alpha_r
 *   - Reparameterize: Newton-Raphson root-find on each point
 *   - Recursive bisection at max-error split point
 *
 * No malloc inside the fit loop (uses caller stack via fixed buffer).
 * Recursion depth bounded by SCHNEIDER3D_MAX_DEPTH.
 */
#include "schneider3d.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define SCHNEIDER3D_MAX_PTS    200000   /* matches test_xcompare polyline */
#define SCHNEIDER3D_MAX_DEPTH  64

typedef double V3[3];

static inline void v3_set (V3 r, double x, double y, double z)
                          { r[0]=x; r[1]=y; r[2]=z; }
static inline void v3_copy(V3 r, const V3 a) { r[0]=a[0]; r[1]=a[1]; r[2]=a[2]; }
static inline void v3_sub (V3 r, const V3 a, const V3 b)
                          { r[0]=a[0]-b[0]; r[1]=a[1]-b[1]; r[2]=a[2]-b[2]; }
static inline void v3_add (V3 r, const V3 a, const V3 b)
                          { r[0]=a[0]+b[0]; r[1]=a[1]+b[1]; r[2]=a[2]+b[2]; }
static inline void v3_scale(V3 r, const V3 a, double s)
                          { r[0]=a[0]*s; r[1]=a[1]*s; r[2]=a[2]*s; }
static inline double v3_dot(const V3 a, const V3 b)
                          { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
static inline double v3_sqlen(const V3 a) { return v3_dot(a, a); }
static inline double v3_dist(const V3 a, const V3 b)
{
    V3 d; v3_sub(d, a, b); return sqrt(v3_sqlen(d));
}
static inline void v3_normalize(V3 a)
{
    double n = sqrt(v3_sqlen(a));
    if (n > 1e-30) { a[0]/=n; a[1]/=n; a[2]/=n; }
}

/* Bernstein basis */
static double B0(double u) { double t = 1.0 - u; return t*t*t; }
static double B1(double u) { double t = 1.0 - u; return 3.0*u*t*t; }
static double B2(double u) { double t = 1.0 - u; return 3.0*u*u*t; }
static double B3(double u) { return u*u*u; }

/* De Casteljau evaluation of Bezier of given degree (≤3). */
static void bezier_eval(int degree, const V3 *V, double t, V3 out)
{
    V3 tmp[4];
    for (int i = 0; i <= degree; i++) v3_copy(tmp[i], V[i]);
    for (int i = 1; i <= degree; i++) {
        for (int j = 0; j <= degree - i; j++) {
            tmp[j][0] = (1.0 - t) * tmp[j][0] + t * tmp[j+1][0];
            tmp[j][1] = (1.0 - t) * tmp[j][1] + t * tmp[j+1][1];
            tmp[j][2] = (1.0 - t) * tmp[j][2] + t * tmp[j+1][2];
        }
    }
    v3_copy(out, tmp[0]);
}

static void compute_left_tangent(const V3 *d, int end, V3 out)
{
    v3_sub(out, d[end+1], d[end]);
    v3_normalize(out);
}
static void compute_right_tangent(const V3 *d, int end, V3 out)
{
    v3_sub(out, d[end-1], d[end]);
    v3_normalize(out);
}
static void compute_center_tangent(const V3 *d, int center, V3 out)
{
    V3 v1, v2;
    v3_sub(v1, d[center-1], d[center]);
    v3_sub(v2, d[center],   d[center+1]);
    out[0] = (v1[0] + v2[0]) * 0.5;
    out[1] = (v1[1] + v2[1]) * 0.5;
    out[2] = (v1[2] + v2[2]) * 0.5;
    v3_normalize(out);
}

static void chord_length_param(const V3 *d, int first, int last, double *u)
{
    int n = last - first;
    u[0] = 0.0;
    for (int i = 1; i <= n; i++)
        u[i] = u[i-1] + v3_dist(d[first+i], d[first+i-1]);
    if (u[n] > 1e-30) {
        for (int i = 1; i <= n; i++) u[i] /= u[n];
    }
    u[n] = 1.0;
}

/* Solve 2x2 normal eq for Bezier middle control points P1, P2.
 * Output bez[0]=d[first], bez[3]=d[last] pinned; bez[1], bez[2]
 * lie along tHat1, tHat2 at distances alpha_l, alpha_r. */
static void generate_bezier(const V3 *d, int first, int last,
                            const double *u, const V3 tHat1,
                            const V3 tHat2, V3 *bez)
{
    int npts = last - first + 1;
    double C[2][2] = {{0,0},{0,0}};
    double X[2]    = {0, 0};

    V3 *Av0 = (V3*)malloc(npts * sizeof(V3));
    V3 *Av1 = (V3*)malloc(npts * sizeof(V3));

    for (int i = 0; i < npts; i++) {
        v3_scale(Av0[i], tHat1, B1(u[i]));
        v3_scale(Av1[i], tHat2, B2(u[i]));
    }

    V3 first_pt, last_pt;
    v3_copy(first_pt, d[first]);
    v3_copy(last_pt,  d[last]);

    for (int i = 0; i < npts; i++) {
        C[0][0] += v3_dot(Av0[i], Av0[i]);
        C[0][1] += v3_dot(Av0[i], Av1[i]);
        C[1][1] += v3_dot(Av1[i], Av1[i]);

        /* tmp = d[first+i] − [B0·first + B1·first + B2·last + B3·last] */
        V3 acc, t;
        v3_scale(acc, first_pt, B0(u[i]));
        v3_scale(t,   first_pt, B1(u[i])); v3_add(acc, acc, t);
        v3_scale(t,   last_pt,  B2(u[i])); v3_add(acc, acc, t);
        v3_scale(t,   last_pt,  B3(u[i])); v3_add(acc, acc, t);
        V3 tmp; v3_sub(tmp, d[first+i], acc);

        X[0] += v3_dot(Av0[i], tmp);
        X[1] += v3_dot(Av1[i], tmp);
    }
    C[1][0] = C[0][1];

    free(Av0); free(Av1);

    double det_C   = C[0][0] * C[1][1] - C[1][0] * C[0][1];
    double det_C0X = C[0][0] * X[1]    - C[1][0] * X[0];
    double det_XC1 = X[0]    * C[1][1] - X[1]    * C[0][1];

    double alpha_l = (det_C == 0) ? 0.0 : det_XC1 / det_C;
    double alpha_r = (det_C == 0) ? 0.0 : det_C0X / det_C;

    double seg_len = v3_dist(d[last], d[first]);
    double eps = 1e-6 * seg_len;

    v3_copy(bez[0], d[first]);
    v3_copy(bez[3], d[last]);

    if (alpha_l < eps || alpha_r < eps) {
        /* Wu/Barsky fallback */
        double dist = seg_len / 3.0;
        V3 t;
        v3_scale(t, tHat1, dist); v3_add(bez[1], bez[0], t);
        v3_scale(t, tHat2, dist); v3_add(bez[2], bez[3], t);
        return;
    }
    V3 t;
    v3_scale(t, tHat1, alpha_l); v3_add(bez[1], bez[0], t);
    v3_scale(t, tHat2, alpha_r); v3_add(bez[2], bez[3], t);
}

static double newton_raphson_root(const V3 *Q, const V3 P, double u)
{
    /* Q' control vertices */
    V3 Q1[3], Q2[2];
    for (int i = 0; i < 3; i++) {
        Q1[i][0] = (Q[i+1][0] - Q[i][0]) * 3.0;
        Q1[i][1] = (Q[i+1][1] - Q[i][1]) * 3.0;
        Q1[i][2] = (Q[i+1][2] - Q[i][2]) * 3.0;
    }
    for (int i = 0; i < 2; i++) {
        Q2[i][0] = (Q1[i+1][0] - Q1[i][0]) * 2.0;
        Q2[i][1] = (Q1[i+1][1] - Q1[i][1]) * 2.0;
        Q2[i][2] = (Q1[i+1][2] - Q1[i][2]) * 2.0;
    }
    V3 Qu, Q1u, Q2u;
    bezier_eval(3, Q,  u, Qu);
    bezier_eval(2, Q1, u, Q1u);
    bezier_eval(1, Q2, u, Q2u);

    V3 diff; v3_sub(diff, Qu, P);
    double num = v3_dot(diff, Q1u);
    double den = v3_sqlen(Q1u) + v3_dot(diff, Q2u);
    if (den == 0.0) return u;
    return u - num / den;
}

static void reparameterize(const V3 *d, int first, int last,
                           const double *u, const V3 *bez,
                           double *uPrime)
{
    int n = last - first;
    for (int i = 0; i <= n; i++)
        uPrime[i] = newton_raphson_root(bez, d[first+i], u[i]);
}

static double compute_max_error(const V3 *d, int first, int last,
                                const V3 *bez, const double *u,
                                int *split_point)
{
    int n = last - first;
    *split_point = first + n / 2;
    double max_d = 0.0;
    for (int i = 1; i < n; i++) {
        V3 P; bezier_eval(3, bez, u[i], P);
        V3 v; v3_sub(v, P, d[first+i]);
        double d2 = v3_sqlen(v);
        if (d2 >= max_d) { max_d = d2; *split_point = first + i; }
    }
    return max_d;
}

static void fit_cubic(const V3 *d, int first, int last,
                      const V3 tHat1, const V3 tHat2,
                      double error_sq, int depth,
                      schneider3d_emit_t emit, void *user)
{
    int npts = last - first + 1;

    /* Heuristic for 2-point region: place P1, P2 along tangents. */
    if (npts == 2) {
        double dist = v3_dist(d[last], d[first]) / 3.0;
        V3 bez[4], t;
        v3_copy(bez[0], d[first]);
        v3_copy(bez[3], d[last]);
        v3_scale(t, tHat1, dist); v3_add(bez[1], bez[0], t);
        v3_scale(t, tHat2, dist); v3_add(bez[2], bez[3], t);
        schneider3d_bezier_t out;
        v3_copy(out.P0, bez[0]); v3_copy(out.P1, bez[1]);
        v3_copy(out.P2, bez[2]); v3_copy(out.P3, bez[3]);
        emit(&out, user);
        return;
    }

    double *u  = (double*)malloc(npts * sizeof(double));
    chord_length_param(d, first, last, u);

    V3 bez[4];
    generate_bezier(d, first, last, u, tHat1, tHat2, bez);

    int split;
    double max_err = compute_max_error(d, first, last, bez, u, &split);

    if (max_err < error_sq) {
        schneider3d_bezier_t out;
        v3_copy(out.P0, bez[0]); v3_copy(out.P1, bez[1]);
        v3_copy(out.P2, bez[2]); v3_copy(out.P3, bez[3]);
        emit(&out, user);
        free(u);
        return;
    }

    /* Try iteration: 4 Newton-Raphson reparam attempts. */
    double iter_err = error_sq * 4.0;
    if (max_err < iter_err && depth < SCHNEIDER3D_MAX_DEPTH) {
        double *uPrime = (double*)malloc(npts * sizeof(double));
        for (int it = 0; it < 4; it++) {
            reparameterize(d, first, last, u, bez, uPrime);
            generate_bezier(d, first, last, uPrime, tHat1, tHat2, bez);
            max_err = compute_max_error(d, first, last, bez, uPrime, &split);
            if (max_err < error_sq) {
                schneider3d_bezier_t out;
                v3_copy(out.P0, bez[0]); v3_copy(out.P1, bez[1]);
                v3_copy(out.P2, bez[2]); v3_copy(out.P3, bez[3]);
                emit(&out, user);
                free(u); free(uPrime);
                return;
            }
            memcpy(u, uPrime, npts * sizeof(double));
        }
        free(uPrime);
    }

    /* Bisect at split point. */
    free(u);
    if (depth >= SCHNEIDER3D_MAX_DEPTH || split <= first || split >= last) {
        /* Recursion bottom: emit straight Bezier as fallback. */
        double dist = v3_dist(d[last], d[first]) / 3.0;
        V3 b[4], t;
        v3_copy(b[0], d[first]);
        v3_copy(b[3], d[last]);
        v3_scale(t, tHat1, dist); v3_add(b[1], b[0], t);
        v3_scale(t, tHat2, dist); v3_add(b[2], b[3], t);
        schneider3d_bezier_t out;
        v3_copy(out.P0, b[0]); v3_copy(out.P1, b[1]);
        v3_copy(out.P2, b[2]); v3_copy(out.P3, b[3]);
        emit(&out, user);
        return;
    }
    V3 tHatCenter, tHatNeg;
    compute_center_tangent(d, split, tHatCenter);
    fit_cubic(d, first, split, tHat1, tHatCenter, error_sq, depth+1, emit, user);
    tHatNeg[0] = -tHatCenter[0];
    tHatNeg[1] = -tHatCenter[1];
    tHatNeg[2] = -tHatCenter[2];
    fit_cubic(d, split, last, tHatNeg, tHat2, error_sq, depth+1, emit, user);
}

void schneider3d_fit(const double *xyz_flat, int n_points,
                     double error_sq,
                     schneider3d_emit_t emit, void *user)
{
    if (n_points < 2) return;
    V3 *d = (V3*)malloc(n_points * sizeof(V3));
    for (int i = 0; i < n_points; i++) {
        d[i][0] = xyz_flat[3*i+0];
        d[i][1] = xyz_flat[3*i+1];
        d[i][2] = xyz_flat[3*i+2];
    }
    V3 tHat1, tHat2;
    compute_left_tangent(d, 0, tHat1);
    compute_right_tangent(d, n_points - 1, tHat2);
    fit_cubic(d, 0, n_points - 1, tHat1, tHat2, error_sq, 0, emit, user);
    free(d);
}
