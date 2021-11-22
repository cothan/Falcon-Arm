#include <arm_neon.h>
#include "macrous.h"
#include "inner.h"
#include "ntt.h"
#include "ntt_consts.h"
/* ===================================================================== */
/*
 * Compute NTT on a ring element.
 */
void mq_NTT(uint16_t *a, unsigned logn)
{
    size_t n, t, m;

    n = (size_t)1 << logn;
    t = n;
    for (m = 1; m < n; m <<= 1)
    {
        size_t ht, i, j1;

        ht = t >> 1;
        for (i = 0, j1 = 0; i < m; i++, j1 += t)
        {
            size_t j, j2;
            uint32_t s;

            s = GMb[m + i];
            j2 = j1 + ht;
            for (j = j1; j < j2; j++)
            {
                uint32_t u, v;

                u = a[j];
                v = mq_montymul(a[j + ht], s);
                a[j] = (uint16_t)mq_add(u, v);
                a[j + ht] = (uint16_t)mq_sub(u, v);
            }
        }
        t = ht;
    }
}

void neon_invNTT(int16_t a[FALCON_N])
{
    // Total SIMD registers: 28 = 16 + 8 + 4
    uint16x8x4_t v0, v1, v2, v3; // 16
    uint16x8x4_t zl, zh; // 8
    uint16x8x4_t t;
    uint16x8_t neon_q;
    unsigned k = 0;
    // TODO: add reduction between NTT level to avoid overflows
    for (unsigned j = 0; j < FALCON_N; j += 128)
    {
        vload_u16_x4(v0, &a[j]);
        vload_u16_x4(v1, &a[j + 32]);
        vload_u16_x4(v2, &a[j + 64]);
        vload_u16_x4(v3, &a[j + 96]);

        // Layer 0
        // v0.val[0]: 0, 4, 8,  12 | 16, 20, 24, 28
        // v0.val[1]: 1, 5, 9,  13 | 17, 21, 25, 29
        // v0.val[2]: 2, 6, 10, 14 | 18, 22, 26, 30
        // v0.val[3]: 3, 7, 11, 15 | 19, 23, 27, 31
        vload_u16_x4(zl, &ntt[k]);
        vload_u16_x4(zh, &ntt[k]);
        k += 32;

        // 0 - 1, 2 - 3
        gsbf(v0.val[0], v0.val[1], zl.val[0], zh.val[0], neon_q, t.val[0]);
        gsbf(v0.val[2], v0.val[3], zl.val[0], zh.val[0], neon_q, t.val[1]);
        gsbf(v1.val[0], v1.val[1], zl.val[0], zh.val[0], neon_q, t.val[2]);
        gsbf(v1.val[2], v1.val[3], zl.val[0], zh.val[0], neon_q, t.val[3]);

        gsbf(v2.val[0], v2.val[1], zl.val[0], zh.val[0], neon_q, t.val[0]);
        gsbf(v2.val[2], v2.val[3], zl.val[0], zh.val[0], neon_q, t.val[1]);
        gsbf(v3.val[0], v3.val[1], zl.val[0], zh.val[0], neon_q, t.val[2]);
        gsbf(v3.val[2], v3.val[3], zl.val[0], zh.val[0], neon_q, t.val[3]);

        // Layer 1
        // v0.val[0]: 0, 4, 8,  12 | 16, 20, 24, 28
        // v0.val[1]: 1, 5, 9,  13 | 17, 21, 25, 29
        // v0.val[2]: 2, 6, 10, 14 | 18, 22, 26, 30
        // v0.val[3]: 3, 7, 11, 15 | 19, 23, 27, 31
        // 0 - 2, 1 - 3
        gsbf(v0.val[0], v0.val[2], zl.val[0], zh.val[0], neon_q, t.val[0]);
        gsbf(v0.val[1], v0.val[3], zl.val[0], zh.val[0], neon_q, t.val[1]);
        gsbf(v1.val[0], v1.val[2], zl.val[0], zh.val[0], neon_q, t.val[2]);
        gsbf(v1.val[1], v1.val[3], zl.val[0], zh.val[0], neon_q, t.val[3]);

        gsbf(v2.val[0], v2.val[2], zl.val[0], zh.val[0], neon_q, t.val[0]);
        gsbf(v2.val[1], v2.val[3], zl.val[0], zh.val[0], neon_q, t.val[1]);
        gsbf(v3.val[0], v3.val[2], zl.val[0], zh.val[0], neon_q, t.val[2]);
        gsbf(v3.val[1], v3.val[3], zl.val[0], zh.val[0], neon_q, t.val[3]);

        // Layer 2
        // Transpose 
        transpose(v0, t);
        transpose(v1, t);
        transpose(v2, t);
        transpose(v3, t);

        // v0.val[0]: 0,  1,  2,  3  | 16,  17,  18,  19
        // v0.val[1]: 4,  5,  6,  7  | 20,  21,  22,  23
        // v0.val[2]: 8,  9,  10, 11 | 24,  25,  26,  27
        // v0.val[3]: 12, 13, 14, 15 | 28,  29,  30,  31
        // 0 - 1, 2 - 3
        gsbf(v0.val[0], v0.val[1], zl.val[0], zh.val[0], neon_q, t.val[0]);
        gsbf(v0.val[2], v0.val[3], zl.val[0], zh.val[0], neon_q, t.val[1]);
        gsbf(v1.val[0], v1.val[1], zl.val[0], zh.val[0], neon_q, t.val[2]);
        gsbf(v1.val[2], v1.val[3], zl.val[0], zh.val[0], neon_q, t.val[3]);

        gsbf(v2.val[0], v2.val[1], zl.val[0], zh.val[0], neon_q, t.val[0]);
        gsbf(v2.val[2], v2.val[3], zl.val[0], zh.val[0], neon_q, t.val[1]);
        gsbf(v3.val[0], v3.val[1], zl.val[0], zh.val[0], neon_q, t.val[2]);
        gsbf(v3.val[2], v3.val[3], zl.val[0], zh.val[0], neon_q, t.val[3]);

        // Layer 3
        // v0.val[0]: 0,  1,  2,  3  | 16,  17,  18,  19
        // v0.val[1]: 4,  5,  6,  7  | 20,  21,  22,  23
        // v0.val[2]: 8,  9,  10, 11 | 24,  25,  26,  27
        // v0.val[3]: 12, 13, 14, 15 | 28,  29,  30,  31
        // 0 - 2, 1 - 3
        gsbf(v0.val[0], v0.val[2], zl.val[0], zh.val[0], neon_q, t.val[0]);
        gsbf(v0.val[1], v0.val[3], zl.val[0], zh.val[0], neon_q, t.val[1]);
        gsbf(v1.val[0], v1.val[2], zl.val[0], zh.val[0], neon_q, t.val[2]);
        gsbf(v1.val[1], v1.val[3], zl.val[0], zh.val[0], neon_q, t.val[3]);

        gsbf(v2.val[0], v2.val[2], zl.val[0], zh.val[0], neon_q, t.val[0]);
        gsbf(v2.val[1], v2.val[3], zl.val[0], zh.val[0], neon_q, t.val[1]);
        gsbf(v3.val[0], v3.val[2], zl.val[0], zh.val[0], neon_q, t.val[2]);
        gsbf(v3.val[1], v3.val[3], zl.val[0], zh.val[0], neon_q, t.val[3]);

        // Layer 4
        // Re-arrange vector
        // v0.val[0]: 0,  1,  2,  3  | 4,  5,  6,  7  
        // v0.val[1]: 16, 17, 18, 19 | 20, 21, 22, 23
        // v0.val[2]: 8,  9,  10, 11 | 12, 13, 14, 15 
        // v0.val[3]: 24, 25, 26, 27 | 28, 29, 30, 31
        // 0 - 1, 2 - 3
        gsbf(v0.val[0], v0.val[1], zl.val[0], zh.val[0], neon_q, t.val[0]);
        gsbf(v0.val[2], v0.val[3], zl.val[0], zh.val[0], neon_q, t.val[1]);
        gsbf(v1.val[0], v1.val[1], zl.val[0], zh.val[0], neon_q, t.val[2]);
        gsbf(v1.val[2], v1.val[3], zl.val[0], zh.val[0], neon_q, t.val[3]);

        gsbf(v2.val[0], v2.val[1], zl.val[0], zh.val[0], neon_q, t.val[0]);
        gsbf(v2.val[2], v2.val[3], zl.val[0], zh.val[0], neon_q, t.val[1]);
        gsbf(v3.val[0], v3.val[1], zl.val[0], zh.val[0], neon_q, t.val[2]);
        gsbf(v3.val[2], v3.val[3], zl.val[0], zh.val[0], neon_q, t.val[3]);

        // Layer 5
        // Cross block
        // v0.0->3 - v1.0->3
        gsbf(v0.val[0], v1.val[0], zl.val[0], zh.val[0], neon_q, t.val[0]);
        gsbf(v0.val[1], v1.val[1], zl.val[0], zh.val[0], neon_q, t.val[1]);
        gsbf(v0.val[2], v1.val[2], zl.val[0], zh.val[0], neon_q, t.val[2]);
        gsbf(v0.val[3], v1.val[3], zl.val[0], zh.val[0], neon_q, t.val[3]);

        gsbf(v2.val[0], v3.val[0], zl.val[0], zh.val[0], neon_q, t.val[0]);
        gsbf(v2.val[1], v3.val[1], zl.val[0], zh.val[0], neon_q, t.val[1]);
        gsbf(v2.val[2], v3.val[2], zl.val[0], zh.val[0], neon_q, t.val[2]);
        gsbf(v2.val[3], v3.val[3], zl.val[0], zh.val[0], neon_q, t.val[3]);

        // Layer 6
        // Cross block
        // v0.0->3 - v2.0->3
        gsbf(v0.val[0], v2.val[0], zl.val[0], zh.val[0], neon_q, t.val[0]);
        gsbf(v0.val[1], v2.val[1], zl.val[0], zh.val[0], neon_q, t.val[1]);
        gsbf(v0.val[2], v2.val[2], zl.val[0], zh.val[0], neon_q, t.val[2]);
        gsbf(v0.val[3], v2.val[3], zl.val[0], zh.val[0], neon_q, t.val[3]);

        gsbf(v1.val[0], v3.val[0], zl.val[0], zh.val[0], neon_q, t.val[0]);
        gsbf(v1.val[1], v3.val[1], zl.val[0], zh.val[0], neon_q, t.val[1]);
        gsbf(v1.val[2], v3.val[2], zl.val[0], zh.val[0], neon_q, t.val[2]);
        gsbf(v1.val[3], v3.val[3], zl.val[0], zh.val[0], neon_q, t.val[3]);
    }

    // Layer 7, 8

    // Layer 9

    // Layer 10
}

/*
 * Compute the inverse NTT on a ring element, binary case.
 */
void mq_iNTT(uint16_t *a, unsigned logn)
{
    size_t n, t, m;
    uint32_t ni;

    n = (size_t)1 << logn;
    t = 1;
    m = n;
    while (m > 1)
    {
        size_t hm, dt, i, j1;

        hm = m >> 1;
        dt = t << 1;
        for (i = 0, j1 = 0; i < hm; i++, j1 += dt)
        {
            size_t j, j2;
            uint32_t s;

            j2 = j1 + t;
            s = iGMb[hm + i];
            for (j = j1; j < j2; j++)
            {
                uint32_t u, v, w;

                u = a[j];
                v = a[j + t];
                a[j] = (uint16_t)mq_add(u, v);
                w = mq_sub(u, v);
                a[j + t] = (uint16_t)
                    mq_montymul(w, s);
            }
        }
        t = dt;
        m = hm;
    }

    /*
	 * To complete the inverse NTT, we must now divide all values by
	 * n (the vector size). We thus need the inverse of n, i.e. we
	 * need to divide 1 by 2 logn times. But we also want it in
	 * Montgomery representation, i.e. we also want to multiply it
	 * by R = 2^16. In the common case, this should be a simple right
	 * shift. The loop below is generic and works also in corner cases;
	 * its computation time is negligible.
	 */
    ni = R;
    for (m = n; m > 1; m >>= 1)
    {
        ni = mq_rshift1(ni);
    }
    for (m = 0; m < n; m++)
    {
        a[m] = (uint16_t)mq_montymul(a[m], ni);
    }
}

/*
 * Reduce a small signed integer modulo q. The source integer MUST
 * be between -q/2 and +q/2.
 */
extern inline uint32_t
mq_conv_small(int x)
{
    /*
	 * If x < 0, the cast to uint32_t will set the high bit to 1.
	 */
    uint32_t y;

    y = (uint32_t)x;
    y += Q & -(y >> 31);
    return y;
}

/*
 * Addition modulo q. Operands must be in the 0..q-1 range.
 */
static inline uint32_t
mq_add(uint32_t x, uint32_t y)
{
    /*
	 * We compute x + y - q. If the result is negative, then the
	 * high bit will be set, and 'd >> 31' will be equal to 1;
	 * thus '-(d >> 31)' will be an all-one pattern. Otherwise,
	 * it will be an all-zero pattern. In other words, this
	 * implements a conditional addition of q.
	 */
    uint32_t d;

    d = x + y - Q;
    d += Q & -(d >> 31);
    return d;
}

/*
 * Subtraction modulo q. Operands must be in the 0..q-1 range.
 */
extern inline uint32_t
mq_sub(uint32_t x, uint32_t y)
{
    /*
	 * As in mq_add(), we use a conditional addition to ensure the
	 * result is in the 0..q-1 range.
	 */
    uint32_t d;

    d = x - y;
    d += Q & -(d >> 31);
    return d;
}

/*
 * Division by 2 modulo q. Operand must be in the 0..q-1 range.
 */
static inline uint32_t
mq_rshift1(uint32_t x)
{
    x += Q & -(x & 1);
    return (x >> 1);
}

/*
 * Montgomery multiplication modulo q. If we set R = 2^16 mod q, then
 * this function computes: x * y / R mod q
 * Operands must be in the 0..q-1 range.
 */
static inline uint32_t
mq_montymul(uint32_t x, uint32_t y)
{
    uint32_t z, w;

    /*
	 * We compute x*y + k*q with a value of k chosen so that the 16
	 * low bits of the result are 0. We can then shift the value.
	 * After the shift, result may still be larger than q, but it
	 * will be lower than 2*q, so a conditional subtraction works.
	 */

    z = x * y;
    w = ((z * Q0I) & 0xFFFF) * Q;

    /*
	 * When adding z and w, the result will have its low 16 bits
	 * equal to 0. Since x, y and z are lower than q, the sum will
	 * be no more than (2^15 - 1) * q + (q - 1)^2, which will
	 * fit on 29 bits.
	 */
    z = (z + w) >> 16;

    /*
	 * After the shift, analysis shows that the value will be less
	 * than 2q. We do a subtraction then conditional subtraction to
	 * ensure the result is in the expected range.
	 */
    z -= Q;
    z += Q & -(z >> 31);
    return z;
}

/*
 * Montgomery squaring (computes (x^2)/R).
 */
static inline uint32_t
mq_montysqr(uint32_t x)
{
    return mq_montymul(x, x);
}

/*
 * Divide x by y modulo q = 12289.
 */
extern inline uint32_t
mq_div_12289(uint32_t x, uint32_t y)
{
    /*
	 * We invert y by computing y^(q-2) mod q.
	 *
	 * We use the following addition chain for exponent e = 12287:
	 *
	 *   e0 = 1
	 *   e1 = 2 * e0 = 2
	 *   e2 = e1 + e0 = 3
	 *   e3 = e2 + e1 = 5
	 *   e4 = 2 * e3 = 10
	 *   e5 = 2 * e4 = 20
	 *   e6 = 2 * e5 = 40
	 *   e7 = 2 * e6 = 80
	 *   e8 = 2 * e7 = 160
	 *   e9 = e8 + e2 = 163
	 *   e10 = e9 + e8 = 323
	 *   e11 = 2 * e10 = 646
	 *   e12 = 2 * e11 = 1292
	 *   e13 = e12 + e9 = 1455
	 *   e14 = 2 * e13 = 2910
	 *   e15 = 2 * e14 = 5820
	 *   e16 = e15 + e10 = 6143
	 *   e17 = 2 * e16 = 12286
	 *   e18 = e17 + e0 = 12287
	 *
	 * Additions on exponents are converted to Montgomery
	 * multiplications. We define all intermediate results as so
	 * many local variables, and let the C compiler work out which
	 * must be kept around.
	 */
    uint32_t y0, y1, y2, y3, y4, y5, y6, y7, y8, y9;
    uint32_t y10, y11, y12, y13, y14, y15, y16, y17, y18;

    y0 = mq_montymul(y, R2);
    y1 = mq_montysqr(y0);
    y2 = mq_montymul(y1, y0);
    y3 = mq_montymul(y2, y1);
    y4 = mq_montysqr(y3);
    y5 = mq_montysqr(y4);
    y6 = mq_montysqr(y5);
    y7 = mq_montysqr(y6);
    y8 = mq_montysqr(y7);
    y9 = mq_montymul(y8, y2);
    y10 = mq_montymul(y9, y8);
    y11 = mq_montysqr(y10);
    y12 = mq_montysqr(y11);
    y13 = mq_montymul(y12, y9);
    y14 = mq_montysqr(y13);
    y15 = mq_montysqr(y14);
    y16 = mq_montymul(y15, y10);
    y17 = mq_montysqr(y16);
    y18 = mq_montymul(y17, y0);

    /*
	 * Final multiplication with x, which is not in Montgomery
	 * representation, computes the correct division result.
	 */
    return mq_montymul(y18, x);
}

/*
 * Convert a polynomial (mod q) to Montgomery representation.
 */
void mq_poly_tomonty(uint16_t *f, unsigned logn)
{
    size_t u, n;

    n = (size_t)1 << logn;
    for (u = 0; u < n; u++)
    {
        f[u] = (uint16_t)mq_montymul(f[u], R2);
    }
}

/*
 * Multiply two polynomials together (NTT representation, and using
 * a Montgomery multiplication). Result f*g is written over f.
 */
void mq_poly_montymul_ntt(uint16_t *f, const uint16_t *g, unsigned logn)
{
    size_t u, n;

    n = (size_t)1 << logn;
    for (u = 0; u < n; u++)
    {
        f[u] = (uint16_t)mq_montymul(f[u], g[u]);
    }
}

/*
 * Subtract polynomial g from polynomial f.
 */
void mq_poly_sub(uint16_t *f, const uint16_t *g, unsigned logn)
{
    size_t u, n;

    n = (size_t)1 << logn;
    for (u = 0; u < n; u++)
    {
        f[u] = (uint16_t)mq_sub(f[u], g[u]);
    }
}

/* ===================================================================== */