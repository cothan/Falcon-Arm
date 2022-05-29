/*
 * Addition of two complex numbers (d = a + b).
 */
#define FPC_ADD(d_re, d_im, a_re, a_im, b_re, b_im) \
    d_re = a_re + b_re;                       \
    d_im = a_im + b_im;

/*
 * Addition of two complex numbers (d = a + jb).
 */
#define FPC_ADDJ(d_re, d_im, a_re, a_im, b_re, b_im) \
    d_re = a_re - b_im;                        \
    d_im = a_im + b_re;
/*
 * Subtraction of two complex numbers (d = a - b).
 */
#define FPC_SUB(d_re, d_im, a_re, a_im, b_re, b_im) \
    d_re = a_re - b_re;                       \
    d_im = a_im - b_im;

/*
 * Subtraction of two complex numbers (d = a - jb).
 */
#define FPC_SUBJ(d_re, d_im, a_re, a_im, b_re, b_im) \
    d_re = a_re + b_im;                        \
    d_im = a_im - b_re;

/*
 * Multplication of two complex numbers (d = a * b).
 */
#define FPC_MUL(d_re, d_im, a_re, a_im, b_re, b_im) \
    d_re = a_re * b_re - a_im * b_im;     \
    d_im = a_re * b_im + a_im * b_re;

/*
 * Multplication of two complex numbers (d = a * conj(b)).
 * a is swapped from: a_re|a_im to a_im|a_re
 * b is swapped from: b_re|b_im to b_im|b_re
 */
#define FPC_MUL_CONJ(d_re, d_im, a_im, a_re, b_im, b_re) \
    d_re = b_re * a_re + a_im * b_im;          \
    d_im = b_im * a_re - a_im * b_re;

/*
 * Multplication of two complex numbers (d = a * conj(jb)).
 */
#define FPC_MUL_CONJ_J(d_re, d_im, a_re, a_im, b_re, b_im) \
    d_re = a_im * b_re - b_im * a_re;            \
    d_im = -(a_im * b_im + b_re * a_re);

/*
 * Multplication of two complex numbers (d = a * - conj(jb)).
 * b is swapped from: b_re|b_im to b_im|b_re
 */
#define FPC_MUL_CONJ_J_m(d_re, d_im, a_re, a_im, b_im, b_re) \
    d_re = a_re * b_re - a_im * b_im;              \
    d_im = a_im * b_re + a_re * b_im;

