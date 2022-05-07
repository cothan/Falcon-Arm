/*
 * Poly FFT
 *
 * =============================================================================
 * Copyright (c) 2021 by Cryptographic Engineering Research Group (CERG)
 * ECE Department, George Mason University
 * Fairfax, VA, U.S.A.
 * Author: Duc Tri Nguyen
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *     http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * =============================================================================
 * @author   Duc Tri Nguyen <dnguye69@gmu.edu>
 */

#include "inner.h"
#include "macrofx4.h"
#include "macrof.h"
// #include <assert.h>

/* see inner.h */
void ZfN(poly_add)(fpr *c, const fpr *restrict a, const fpr *restrict b, unsigned logn)
{
    float64x2x4_t neon_a, neon_b, neon_c;
    float64x2x2_t neon_a2, neon_b2, neon_c2;
    const int falcon_n = 1 << logn;
    switch (logn)
    {
    case 1:
        // n = 2;
        vload(neon_a.val[0], &a[0]);
        vload(neon_b.val[0], &b[0]);

        vfadd(neon_c.val[0], neon_a.val[0], neon_b.val[0]);

        vstore(&c[0], neon_c.val[0]);
        break;

    case 2:
        // n = 4
        vloadx2(neon_a2, &a[0]);
        vloadx2(neon_b2, &b[0]);

        vfadd(neon_c2.val[0], neon_a2.val[0], neon_b2.val[0]);
        vfadd(neon_c2.val[1], neon_a2.val[1], neon_b2.val[1]);

        vstorex2(&c[0], neon_c2);
        break;

    default:
        for (int i = 0; i < falcon_n; i += 8)
        {
            vloadx4(neon_a, &a[i]);
            vloadx4(neon_b, &b[i]);

            vfaddx4(neon_c, neon_a, neon_b);

            vstorex4(&c[i], neon_c);
        }
        break;
    }
}

/* see inner.h */
/*
 * c = a - b
 */
void ZfN(poly_sub)(fpr *c, const fpr *restrict a, const fpr *restrict b, unsigned logn)
{
    float64x2x4_t neon_a, neon_b, neon_c;
    float64x2x2_t neon_a2, neon_b2, neon_c2;
    const int falcon_n = 1 << logn;
    switch (logn)
    {
    case 1:
        vload(neon_a.val[0], &a[0]);
        vload(neon_b.val[0], &b[0]);

        vfsub(neon_c.val[0], neon_a.val[0], neon_b.val[0]);

        vstore(&c[0], neon_c.val[0]);
        break;

    case 2:
        vloadx2(neon_a2, &a[0]);
        vloadx2(neon_b2, &b[0]);

        vfsub(neon_c2.val[0], neon_a2.val[0], neon_b2.val[0]);
        vfsub(neon_c2.val[1], neon_a2.val[1], neon_b2.val[1]);

        vstorex2(&c[0], neon_c2);
        break;

    default:
        for (int i = 0; i < falcon_n; i += 8)
        {
            vloadx4(neon_a, &a[i]);
            vloadx4(neon_b, &b[i]);

            vfsubx4(neon_c, neon_a, neon_b);

            vstorex4(&c[i], neon_c);
        }
        break;
    }
}

/* see inner.h */
/*
 * c = -a
 */
void ZfN(poly_neg)(fpr *c, const fpr *restrict a, unsigned logn)
{
    float64x2x4_t neon_a, neon_c;
    float64x2x2_t neon_a2, neon_c2;
    const int falcon_n = 1 << logn;

    switch (logn)
    {
    case 1:
        vload(neon_a.val[0], &a[0]);

        vfneg(neon_c.val[0], neon_a.val[0]);

        vstore(&c[0], neon_c.val[0]);
        break;

    case 2:
        vloadx2(neon_a2, &a[0]);

        vfneg(neon_c2.val[0], neon_a2.val[0]);
        vfneg(neon_c2.val[1], neon_a2.val[1]);

        vstorex2(&c[0], neon_c2);
        break;

    default:
        for (int i = 0; i < falcon_n; i += 8)
        {
            vloadx4(neon_a, &a[i]);

            vfnegx4(neon_c, neon_a);

            vstorex4(&c[i], neon_c);
        }
        break;
    }
}

/* see inner.h */
void ZfN(poly_adj_fft)(fpr *c, const fpr *restrict a, unsigned logn)
{

    float64x2x4_t neon_a, neon_c;
    float64x2x2_t neon_a2, neon_c2;
    const int falcon_n = 1 << logn;
    const int hn = falcon_n >> 1;

    switch (logn)
    {
    case 1:
        // n = 2; hn = 1;
        c[1] = fpr_neg(a[1]);
        break;

    case 2:
        // n = 4; hn = 2
        vload(neon_a.val[0], &a[2]);
        vfneg(neon_c.val[0], neon_a.val[0]);
        vstore(&c[2], neon_c.val[0]);
        break;

    case 3:
        // n = 8; hn = 4
        vloadx2(neon_a2, &a[4]);
        vfneg(neon_c2.val[0], neon_a2.val[0]);
        vfneg(neon_c2.val[1], neon_a2.val[1]);
        vstorex2(&c[4], neon_c2);
        break;

    default:
        for (int i = hn; i < falcon_n; i += 8)
        {
            vloadx4(neon_a, &a[i]);

            vfnegx4(neon_c, neon_a);

            vstorex4(&c[i], neon_c);
        }
        break;
    }
}

static inline void ZfN(poly_mul_fft_log1)(fpr *restrict c, const fpr *restrict a, const fpr *restrict b)
{
    // n = 2
    float64x2_t neon_a, neon_b, neon_c;
    // re | im
    vload(neon_a, &a[0]);
    vload(neon_b, &b[0]);

    vfmul_lane(neon_c, neon_b, neon_a, 0);
    vfcmla_90(neon_c, neon_a, neon_b);

    vstore(&c[0], neon_c);
}

static inline void ZfN(poly_mul_fft_log2)(fpr *restrict c, const fpr *restrict a, const fpr *restrict b)
{
    // n = 4
    float64x2x2_t neon_a, neon_b, neon_c;
    float64x2_t a_re, a_im, b_re, b_im, c_re, c_im;

    // 0: re, re
    // 1: im, im
    vloadx2(neon_a, &a[0]);
    vloadx2(neon_b, &b[0]);

    a_re = neon_a.val[0];
    a_im = neon_a.val[1];
    b_re = neon_b.val[0];
    b_im = neon_b.val[1];

    vfmul(c_re, a_re, b_re);
    vfmul(c_im, a_re, b_im);

    vfms(c_re, c_re, a_im, b_im);
    vfma(c_im, c_im, a_im, b_re);

    neon_c.val[0] = c_re;
    neon_c.val[1] = c_im;

    vstorex2(&c[0], neon_c);
}

static inline void ZfN(poly_mul_fft_log3)(fpr *restrict c, const fpr *restrict a, const fpr *restrict b)
{
    // n = 8
    float64x2x4_t neon_a, neon_b, neon_c;
    float64x2x2_t a_re, a_im, b_re, b_im, c_re, c_im;

    vloadx4(neon_a, &a[0]);
    vloadx4(neon_b, &b[0]);

    a_re.val[0] = neon_a.val[0];
    a_re.val[1] = neon_a.val[1];
    a_im.val[0] = neon_a.val[2];
    a_im.val[1] = neon_a.val[3];

    b_re.val[0] = neon_b.val[0];
    b_re.val[1] = neon_b.val[1];
    b_im.val[0] = neon_b.val[2];
    b_im.val[1] = neon_b.val[3];

    vfmul(c_re.val[0], a_re.val[0], b_re.val[0]);
    vfmul(c_re.val[1], a_re.val[1], b_re.val[1]);
    vfmul(c_im.val[0], a_re.val[0], b_im.val[0]);
    vfmul(c_im.val[1], a_re.val[1], b_im.val[1]);

    vfms(c_re.val[0], c_re.val[0], a_im.val[0], b_im.val[0]);
    vfms(c_re.val[1], c_re.val[1], a_im.val[1], b_im.val[1]);
    vfma(c_im.val[0], c_im.val[0], a_im.val[0], b_re.val[0]);
    vfma(c_im.val[1], c_im.val[1], a_im.val[1], b_re.val[1]);

    neon_c.val[0] = c_re.val[0];
    neon_c.val[1] = c_re.val[1];
    neon_c.val[2] = c_im.val[0];
    neon_c.val[3] = c_im.val[1];

    vstorex4(&c[0], neon_c);
}

/* see inner.h */
/*
 * c = a * b
 */
void ZfN(poly_mul_fft)(fpr *c, const fpr *a, const fpr *restrict b, unsigned logn)
{
    // Total 32 registers
    float64x2x4_t a_re, b_re, a_im, b_im; // 24
    float64x2x4_t c_re, c_im;             // 8
    const int falcon_n = 1 << logn;
    const int hn = falcon_n >> 1;
    switch (logn)
    {
    case 1:
        ZfN(poly_mul_fft_log1)(c, a, b);
        break;

    case 2:
        ZfN(poly_mul_fft_log2)(c, a, b);
        break;

    case 3:
        ZfN(poly_mul_fft_log3)(c, a, b);
        break;

    default:
        for (int i = 0; i < hn; i += 8)
        {
            vloadx4(a_re, &a[i]);
            vloadx4(a_im, &a[i + hn]);
            vloadx4(b_re, &b[i]);
            vloadx4(b_im, &b[i + hn]);

            vfmulx4(c_re, a_re, b_re);
            vfmulx4(c_im, a_re, b_im);

            vfmsx4(c_re, c_re, a_im, b_im);
            vfmax4(c_im, c_im, a_im, b_re);

            vstorex4(&c[i], c_re);
            vstorex4(&c[i + hn], c_im);
        }
        break;
    }
}

static inline void ZfN(poly_mul_fft_add_log1)(fpr *restrict c, const fpr *restrict d,
                                              const fpr *restrict a, const fpr *restrict b)
{
    // n = 2
    float64x2_t neon_a, neon_b, neon_c, neon_d;

    // re | im
    vload(neon_a, &a[0]);
    vload(neon_b, &b[0]);
    vload(neon_d, &d[0]);

    vfma_lane(neon_c, neon_d, neon_b, neon_a, 0);
    vfcmla_90(neon_c, neon_a, neon_b);

    vstore(&c[0], neon_c);
}

static inline void ZfN(poly_mul_fft_add_log2)(fpr *restrict c, const fpr *restrict d,
                                              const fpr *restrict a, const fpr *restrict b)
{
    // n = 4
    float64x2x2_t neon_a, neon_b, neon_c, neon_d;
    float64x2_t a_re, a_im, b_re, b_im, c_re, c_im, d_re, d_im;

    // 0: re, re
    // 1: im, im
    vloadx2(neon_a, &a[0]);
    vloadx2(neon_b, &b[0]);
    vloadx2(neon_d, &d[0]);

    a_re = neon_a.val[0];
    a_im = neon_a.val[1];
    b_re = neon_b.val[0];
    b_im = neon_b.val[1];
    d_re = neon_d.val[0];
    d_im = neon_d.val[1];

    vfma(c_re, d_re, a_re, b_re);
    vfma(c_im, d_im, a_re, b_im);

    vfms(c_re, c_re, a_im, b_im);
    vfma(c_im, c_im, a_im, b_re);

    neon_c.val[0] = c_re;
    neon_c.val[1] = c_im;

    vstorex2(&c[0], neon_c);
}

static inline void ZfN(poly_mul_fft_add_log3)(fpr *restrict c, const fpr *restrict d,
                                              const fpr *restrict a, const fpr *restrict b)
{
    // n = 8
    float64x2x4_t neon_a, neon_b, neon_c, neon_d;
    float64x2x2_t a_re, a_im, b_re, b_im, c_re, c_im, d_re, d_im;

    vloadx4(neon_a, &a[0]);
    vloadx4(neon_b, &b[0]);
    vloadx4(neon_d, &d[0]);

    a_re.val[0] = neon_a.val[0];
    a_re.val[1] = neon_a.val[1];
    a_im.val[0] = neon_a.val[2];
    a_im.val[1] = neon_a.val[3];

    b_re.val[0] = neon_b.val[0];
    b_re.val[1] = neon_b.val[1];
    b_im.val[0] = neon_b.val[2];
    b_im.val[1] = neon_b.val[3];

    d_re.val[0] = neon_d.val[0];
    d_re.val[1] = neon_d.val[1];
    d_im.val[0] = neon_d.val[2];
    d_im.val[1] = neon_d.val[3];

    vfma(c_re.val[0], d_re.val[0], a_re.val[0], b_re.val[0]);
    vfma(c_re.val[1], d_re.val[1], a_re.val[1], b_re.val[1]);
    vfma(c_im.val[0], d_im.val[0], a_re.val[0], b_im.val[0]);
    vfma(c_im.val[1], d_im.val[1], a_re.val[1], b_im.val[1]);

    vfms(c_re.val[0], c_re.val[0], a_im.val[0], b_im.val[0]);
    vfms(c_re.val[1], c_re.val[1], a_im.val[1], b_im.val[1]);
    vfma(c_im.val[0], c_im.val[0], a_im.val[0], b_re.val[0]);
    vfma(c_im.val[1], c_im.val[1], a_im.val[1], b_re.val[1]);

    neon_c.val[0] = c_re.val[0];
    neon_c.val[1] = c_re.val[1];
    neon_c.val[2] = c_im.val[0];
    neon_c.val[3] = c_im.val[1];

    vstorex4(&c[0], neon_c);
}

/* see inner.h */
/*
 * c = d + a * b
 */
void ZfN(poly_mul_add_fft)(fpr *c, const fpr *restrict d,
                           const fpr *a, const fpr *restrict b, unsigned logn)
{
    // Total 32 registers
    float64x2x4_t a_re, b_re, a_im, b_im, d_re, d_im; // 32
    float64x2x4_t c_re, c_im;                         // 8
    const int falcon_n = 1 << logn;
    const int hn = falcon_n >> 1;
    switch (logn)
    {
    case 1:
        ZfN(poly_mul_fft_add_log1)(c, d, a, b);
        break;

    case 2:
        ZfN(poly_mul_fft_add_log2)(c, d, a, b);
        break;

    case 3:
        ZfN(poly_mul_fft_add_log3)(c, d, a, b);
        break;

    default:
        for (int i = 0; i < hn; i += 8)
        {
            vloadx4(a_re, &a[i]);
            vloadx4(a_im, &a[i + hn]);
            vloadx4(b_re, &b[i]);
            vloadx4(b_im, &b[i + hn]);
            vloadx4(d_re, &d[i]);
            vloadx4(d_im, &d[i + hn]);

            vfmax4(c_re, d_re, a_re, b_re);
            vfmax4(c_im, d_im, a_re, b_im);

            vfmsx4(c_re, c_re, a_im, b_im);
            vfmax4(c_im, c_im, a_im, b_re);

            vstorex4(&c[i], c_re);
            vstorex4(&c[i + hn], c_im);
        }
        break;
    }
}

/* see inner.h */
void ZfN(poly_muladj_fft)(fpr *d, fpr *a, const fpr *restrict b, unsigned logn)
{
    // assert(logn >= 4);
    float64x2x4_t a_re, b_re, d_re, a_im, b_im, d_im; // 24
    const int falcon_n = 1 << logn;
    const int hn = falcon_n >> 1;
    for (int i = 0; i < hn; i += 8)
    {
        vloadx4(a_re, &a[i]);
        vloadx4(a_im, &a[i + hn]);
        vloadx4(b_re, &b[i]);
        vloadx4(b_im, &b[i + hn]);

        vfmulx4(d_re, a_re, b_re);
        vfmulx4(d_im, a_im, b_re);

        vfmax4(d_re, d_re, a_im, b_im);
        vfmsx4(d_im, d_im, a_re, b_im);

        vstorex4(&d[i], d_re);
        vstorex4(&d[i + hn], d_im);
    }
}

// c = d + a*b
void ZfN(poly_muladj_add_fft)(fpr *c, fpr *d,
                              const fpr *a, const fpr *restrict b, unsigned logn)
{
    // assert(logn >= 4);
    float64x2x4_t a_re, b_re, d_re, a_im, b_im, d_im; // 24
    const int falcon_n = 1 << logn;
    const int hn = falcon_n >> 1;
    for (int i = 0; i < hn; i += 8)
    {
        vloadx4(a_re, &a[i]);
        vloadx4(a_im, &a[i + hn]);
        vloadx4(b_re, &b[i]);
        vloadx4(b_im, &b[i + hn]);
        vloadx4(d_re, &d[i]);
        vloadx4(d_im, &d[i + hn]);

        vfmax4(d_re, d_re, a_re, b_re);
        vfmax4(d_im, d_im, a_im, b_re);

        vfmax4(d_re, d_re, a_im, b_im);
        vfmsx4(d_im, d_im, a_re, b_im);

        vstorex4(&c[i], d_re);
        vstorex4(&c[i + hn], d_im);
    }
}

/* see inner.h */
/*
 * c = a * adj(a)
 */
void ZfN(poly_mulselfadj_fft)(fpr *c, const fpr *restrict a, unsigned logn)
{
    // assert(logn >= 4);
    /*
     * Since each coefficient is multiplied with its own conjugate,
     * the result contains only real values.
     */
    float64x2x4_t a_re, a_im, c_re, c_im; // 16
    const int falcon_n = 1 << logn;
    const int hn = falcon_n >> 1;

    vfdupx4(c_im, 0);

    for (int i = 0; i < hn; i += 8)
    {
        vloadx4(a_re, &a[i]);
        vloadx4(a_im, &a[i + hn]);

        vfmulx4(c_re, a_re, a_re);
        vfmax4(c_re, c_re, a_im, a_im);

        vstorex4(&c[i], c_re);
        vstorex4(&c[i + hn], c_im);
    }
}

/*
 * c = d + a * adj(a)
 */
void ZfN(poly_mulselfadj_add_fft)(fpr *c, const fpr *restrict d, const fpr *restrict a, unsigned logn)
{
    // assert(logn >= 4);
    /*
     * Since each coefficient is multiplied with its own conjugate,
     * the result contains only real values.
     */
    float64x2x4_t a_re, a_im, c_re, d_re; // 16
    const int falcon_n = 1 << logn;
    const int hn = falcon_n >> 1;

    for (int i = 0; i < hn; i += 8)
    {
        vloadx4(a_re, &a[i]);
        vloadx4(a_im, &a[i + hn]);
        vloadx4(d_re, &d[i]);

        vfmax4(c_re, d_re, a_re, a_re);
        vfmax4(c_re, c_re, a_im, a_im);

        vstorex4(&c[i], c_re);
    }
}

/* see inner.h */
/*
 * c = a * scalar_x
 */
void ZfN(poly_mulconst)(fpr *c, const fpr *a, const fpr x, unsigned logn)
{
    // assert(logn >= 3);
    // Total 9 registers
    float64x2x4_t neon_a, neon_c;
    const int falcon_n = 1 << logn;
    for (int i = 0; i < falcon_n; i += 8)
    {
        vloadx4(neon_a, &a[i]);

        vfmulnx4(neon_c, neon_a, x);

        vstorex4(&c[i], neon_c);
    }
}

/* see inner.h */
void ZfN(poly_div_fft)(fpr *restrict c, const fpr *restrict a, const fpr *restrict b, unsigned logn)
{
    // assert(logn >= 4);
    const int falcon_n = 1 << logn;
    const int hn = falcon_n >> 1;
    float64x2x4_t a_re, a_im, b_re, b_im, c_re, c_im, m;
    for (int i = 0; i < hn; i += 8)
    {
        vloadx4(a_re, &a[i]);
        vloadx4(a_im, &a[i + hn]);
        vloadx4(b_re, &b[i]);
        vloadx4(b_im, &b[i + hn]);

        vfmulx4(m, b_re, b_re);
        vfmulx4(c_re, a_re, b_re);

        vfmax4(m, m, b_im, b_im);
        vfmax4(c_re, c_re, a_im, b_im);

        vfinvx4(m, m);

        vfmulx4(c_im, a_im, b_re);
        vfmsx4(c_im, c_im, a_re, b_im);

        vfmulx4(c_re, c_re, m);
        vfmulx4(c_im, c_im, m);

        vstorex4(&c[i], c_re);
        vstorex4(&c[i + hn], c_im);
    }
}

/* see inner.h */
void ZfN(poly_invnorm2_fft)(fpr *restrict d, const fpr *restrict a, const fpr *restrict b, unsigned logn)
{
    const int falcon_n = 1 << logn;
    const int hn = falcon_n >> 1;
    float64x2x4_t a_re, a_im, b_re, b_im, c_re, c_im, d_re;
    float64x2x2_t x_re_im, y_re_im, z_re_im;

    switch (logn)
    {
    case 1:
        // n = 2; hn = 1; i = 0
        /*
         * x_re = a[0];
         * x_im = a[1];
         * y_re = b[0];
         * y_im = b[1];
         * d[0] = 1.0/( (x_re*x_re) + (x_im*x_im) + (y_re*y_re) + (y_im*y_im) );
         */
        vload(a_re.val[0], &a[0]);
        vload(b_re.val[0], &b[0]);
        vfmul(a_re.val[0], a_re.val[0], a_re.val[0]);
        vfma(c_re.val[0], a_re.val[0], b_re.val[0], b_re.val[0]);
        d[0] = 1.0 / vaddvq_f64(c_re.val[0]);
        break;

    case 2:
        // n = 4; hn = 2; i = 0, 1
        vloadx2(x_re_im, &a[0]);
        vloadx2(y_re_im, &b[0]);

        vfmul(z_re_im.val[0], x_re_im.val[0], x_re_im.val[0]);
        vfmul(z_re_im.val[1], y_re_im.val[0], y_re_im.val[0]);

        vfma(z_re_im.val[0], z_re_im.val[0], x_re_im.val[1], x_re_im.val[1]);
        vfma(z_re_im.val[1], z_re_im.val[1], y_re_im.val[1], y_re_im.val[1]);
        vfadd(z_re_im.val[0], z_re_im.val[0], z_re_im.val[1]);

        vfinv(z_re_im.val[0], z_re_im.val[0]);

        vstore(&d[0], z_re_im.val[0]);
        break;

    case 3:
        // n = 8; hn = 4; i = 0,1,2,3
        vloadx4(a_re, &a[0]);
        vloadx4(b_re, &b[0]);

        vfmulx4(c_re, a_re, a_re);
        vfmax4(c_re, c_re, b_re, b_re);

        vfadd(z_re_im.val[0], c_re.val[0], c_re.val[2]);
        vfadd(z_re_im.val[1], c_re.val[1], c_re.val[3]);

        vfinv(z_re_im.val[0], z_re_im.val[0]);
        vfinv(z_re_im.val[1], z_re_im.val[1]);

        vstorex2(&d[0], z_re_im);
        break;

    default:
        // assert(logn >= 4);
        for (int i = 0; i < hn; i += 8)
        {
            vloadx4(a_re, &a[i]);
            vloadx4(a_im, &a[i + hn]);
            vloadx4(b_re, &b[i]);
            vloadx4(b_im, &b[i + hn]);

            vfmulx4(c_re, a_re, a_re);
            vfmulx4(c_im, b_re, b_re);

            vfmax4(c_re, c_re, a_im, a_im);
            vfmax4(c_im, c_im, b_im, b_im);

            vfaddx4(d_re, c_re, c_im);
            vfinvx4(d_re, d_re);

            vstorex4(&d[i], d_re);
        }
        break;
    }
}

/* see inner.h */
void ZfN(poly_add_muladj_fft)(fpr *restrict d,
                              const fpr *restrict F, const fpr *restrict G,
                              const fpr *restrict f, const fpr *restrict g, unsigned logn)
{
    // assert(logn >= 4);
    const int falcon_n = 1 << logn;
    const int hn = falcon_n >> 1;
    float64x2x4_t F_re, F_im, G_re, G_im;
    float64x2x4_t f_re, f_im, g_re, g_im;
    float64x2x4_t a_re, a_im, b_re, b_im;

    for (int i = 0; i < hn; i += 8)
    {
        vloadx4(F_re, &F[i]);
        vloadx4(F_im, &F[i + hn]);
        vloadx4(G_re, &G[i]);
        vloadx4(G_im, &G[i + hn]);
        vloadx4(f_re, &f[i]);
        vloadx4(f_im, &f[i + hn]);
        vloadx4(g_re, &g[i]);
        vloadx4(g_im, &g[i + hn]);

        vfmulx4(a_re, F_re, f_re);
        vfmulx4(a_im, F_im, f_re);
        vfmulx4(b_re, G_re, g_re);
        vfmulx4(b_im, G_im, g_re);

        vfmax4(a_re, a_re, F_im, f_im);
        vfmax4(b_re, b_re, G_im, g_im);
        vfmsx4(a_im, a_im, F_re, f_im);
        vfmsx4(b_im, b_im, G_re, g_im);

        vfaddx4(a_re, a_re, b_re);
        vfaddx4(a_im, a_im, b_im);

        vstorex4(&d[i], a_re);
        vstorex4(&d[i + hn], a_im);
    }
}

/* see inner.h */
void ZfN(poly_mul_autoadj_fft)(fpr *c, const fpr *a, const fpr *restrict b, unsigned logn)
{
    const int falcon_n = 1 << logn;
    const int hn = falcon_n >> 1;
    float64x2x4_t a_re, a_im, b_re, c_re, c_im;
    float64x2x2_t a_re_im, b_re_im, c_re_im;
    switch (logn)
    {
    case 1:
        // n = 2; hn = 1; i = 0
        // c[0] = a[0] * b[0];
        // c[1] = a[1] * b[0];
        vload(a_re.val[0], &a[0]);
        vfmuln(a_re.val[0], a_re.val[0], b[0]);
        vstore(&c[0], a_re.val[0]);
        break;

    case 2:
        // n = 4; hn = 2; i = 0, 1
        // c[0] = a[0] * b[0];
        // c[2] = a[2] * b[0];
        // c[1] = a[1] * b[1];
        // c[3] = a[3] * b[1];
        vload2(a_re_im, &a[0]);
        vload(b_re_im.val[0], &b[0]);
        vfmul_lane(c_re_im.val[0], a_re_im.val[0], b_re_im.val[0], 0);
        vfmul_lane(c_re_im.val[1], a_re_im.val[1], b_re_im.val[0], 1);
        vstore2(&c[0], c_re_im);
        break;

    case 3:
        // n = 8; hn = 4; i = 0,1,2,3
        // c[0] = a[0] * b[0];
        // c[4] = a[4] * b[0];
        // c[1] = a[1] * b[1];
        // c[5] = a[5] * b[1];
        // c[2] = a[2] * b[2];
        // c[6] = a[6] * b[2];
        // c[3] = a[3] * b[3];
        // c[7] = a[7] * b[3];
        vload4(a_re, &a[0]);
        vloadx2(b_re_im, &b[0]);
        vfmul_lane(c_re.val[0], a_re.val[0], b_re_im.val[0], 0);
        vfmul_lane(c_re.val[1], a_re.val[1], b_re_im.val[0], 1);
        vfmul_lane(c_re.val[2], a_re.val[2], b_re_im.val[1], 0);
        vfmul_lane(c_re.val[3], a_re.val[3], b_re_im.val[1], 1);
        vstore4(&c[0], c_re);
        break;

    default:
        // assert(logn >= 4);
        for (int i = 0; i < hn; i += 8)
        {
            vloadx4(a_re, &a[i]);
            vloadx4(a_im, &a[i + hn]);
            vloadx4(b_re, &b[i]);

            vfmulx4(c_re, a_re, b_re);
            vfmulx4(c_im, a_im, b_re);

            vstorex4(&c[i], c_re);
            vstorex4(&c[i + hn], c_im);
        }
        break;
    }
}

/* see inner.h */
void ZfN(poly_div_autoadj_fft)(fpr *c, const fpr *a, const fpr *restrict b, unsigned logn)
{
    const int falcon_n = 1 << logn;
    const int hn = falcon_n >> 1;
    float64x2x4_t a_re, a_im, b_re, binv, c_re, c_im;

    for (int i = 0; i < hn; i += 8)
    {
        vloadx4(a_re, &a[i]);
        vloadx4(a_im, &a[i + hn]);
        vloadx4(b_re, &b[i]);

        vfinvx4(binv, b_re);

        vfmulx4(c_re, a_re, binv);
        vfmulx4(c_im, a_im, binv);

        vstorex4(&c[i], c_re);
        vstorex4(&c[i + hn], c_im);
    }
}

/* see inner.h */
void ZfN(poly_LDL_fft)(const fpr *restrict g00, fpr *restrict g01, fpr *restrict g11, unsigned logn)
{
    const int falcon_n = 1 << logn;
    const int hn = falcon_n >> 1;
    float64x2x4_t g00_re, g00_im, g01_re, g01_im, g11_re, g11_im;
    float64x2x4_t mu_re, mu_im, m, d_re, d_im;
    float64x2x2_t tmp;
    float64x2_t neon_1i2;

    const fpr imagine[2] = {1.0, -1.0};
    switch (logn)
    {
    case 1:
        // n = 2; hn = 1;
        vload(g00_re.val[0], &g00[0]);
        vload(g01_re.val[0], &g01[0]);
        vload(g11_re.val[0], &g11[0]);
        vload(neon_1i2, &imagine[0]);

        // g00_re^2 | g00_im^2
        vfmul(m.val[0], g00_re.val[0], g00_re.val[0]);
        // 1 / ( g00_re^2 + g00_im^2 )
        m.val[0] = vdupq_n_f64(1 / vaddvq_f64(m.val[0]));

        // re: g01_re * g00_re + g01_im * g00_im
        // im: g01_im * g00_re - g01_re * g00_im
        vfmul_lane(mu_re.val[0], g01_re.val[0], g00_re.val[0], 0);
        vfcmla_270(mu_re.val[0], g00_re.val[0], g01_re.val[0]);
        vfmul(mu_re.val[0], mu_re.val[0], m.val[0]);

        vswap(mu_re.val[1], mu_re.val[0]);
        // im: mu_im * g01_re - mu_re * g01_im
        // re: mu_im * g01_im + mu_re * g01_re
        vfmul_lane(g01_re.val[1], g01_re.val[0], mu_re.val[1], 0);
        vfcmla_90(g01_re.val[1], mu_re.val[1], g01_re.val[0]);

        vswap(g01_re.val[0], g01_re.val[1]);

        vfsub(g11_re.val[0], g11_re.val[0], g01_re.val[0]);
        vfmul(mu_re.val[0], mu_re.val[0], neon_1i2);

        vstore(&g11[0], g11_re.val[0]);
        vstore(&g01[0], mu_re.val[0]);

        break;

    case 2:
        // n = 4; hn = 2
        vloadx2(tmp, &g00[0]);
        g00_re.val[0] = tmp.val[0];
        g00_im.val[0] = tmp.val[1];
        vloadx2(tmp, &g01[0]);
        g01_re.val[0] = tmp.val[0];
        g01_im.val[0] = tmp.val[1];
        vloadx2(tmp, &g11[0]);
        g11_re.val[0] = tmp.val[0];
        g11_im.val[0] = tmp.val[1];

        vfmul(m.val[0], g00_re.val[0], g00_re.val[0]);
        vfma(m.val[0], m.val[0], g00_im.val[0], g00_im.val[0]);
        vfinv(m.val[0], m.val[0]);

        vfmul(mu_re.val[0], g01_re.val[0], g00_re.val[0]);
        vfmul(mu_im.val[0], g01_im.val[0], g00_re.val[0]);

        vfma(mu_re.val[0], mu_re.val[0], g01_im.val[0], g00_im.val[0]);
        vfms(mu_im.val[0], mu_im.val[0], g01_re.val[0], g00_im.val[0]);

        vfmul(mu_re.val[0], mu_re.val[0], m.val[0]);
        vfmul(mu_im.val[0], mu_im.val[0], m.val[0]);

        vfmul(d_re.val[0], mu_re.val[0], g01_re.val[0]);
        vfmul(d_im.val[0], mu_im.val[0], g01_re.val[0]);

        vfma(d_re.val[0], d_re.val[0], mu_im.val[0], g01_im.val[0]);
        vfms(d_im.val[0], d_im.val[0], mu_re.val[0], g01_im.val[0]);

        vfsub(g11_re.val[0], g11_re.val[0], d_re.val[0]);
        vfsub(g11_im.val[0], g11_im.val[0], d_im.val[0]);

        tmp.val[0] = g11_re.val[0];
        tmp.val[1] = g11_im.val[0];
        vstorex2(&g11[0], tmp);

        vfneg(mu_im.val[0], mu_im.val[0]);
        tmp.val[0] = mu_re.val[0];
        tmp.val[1] = mu_im.val[0];
        vstorex2(&g01[0], tmp);

        break;

    case 3:
        //  n = 8; hn = 4
        vloadx4(g00_re, &g00[0]);
        vloadx4(g01_re, &g01[0]);
        vloadx4(g11_re, &g11[0]);
        g00_im.val[0] = g00_re.val[2];
        g00_im.val[1] = g00_re.val[3];
        g01_im.val[0] = g01_re.val[2];
        g01_im.val[1] = g01_re.val[3];

        vfmul(m.val[0], g00_re.val[0], g00_re.val[0]);
        vfmul(m.val[1], g00_re.val[1], g00_re.val[1]);
        vfma(m.val[0], m.val[0], g00_im.val[0], g00_im.val[0]);
        vfma(m.val[1], m.val[1], g00_im.val[1], g00_im.val[1]);
        vfinv(m.val[0], m.val[0]);
        vfinv(m.val[1], m.val[1]);

        vfmul(mu_re.val[0], g01_re.val[0], g00_re.val[0]);
        vfmul(mu_re.val[1], g01_re.val[1], g00_re.val[1]);
        vfmul(mu_im.val[0], g01_im.val[0], g00_re.val[0]);
        vfmul(mu_im.val[1], g01_im.val[1], g00_re.val[1]);

        vfma(mu_re.val[0], mu_re.val[0], g01_im.val[0], g00_im.val[0]);
        vfma(mu_re.val[1], mu_re.val[1], g01_im.val[1], g00_im.val[1]);
        vfms(mu_im.val[0], mu_im.val[0], g01_re.val[0], g00_im.val[0]);
        vfms(mu_im.val[1], mu_im.val[1], g01_re.val[1], g00_im.val[1]);

        vfmul(mu_re.val[0], mu_re.val[0], m.val[0]);
        vfmul(mu_re.val[1], mu_re.val[1], m.val[1]);
        vfmul(mu_im.val[0], mu_im.val[0], m.val[0]);
        vfmul(mu_im.val[1], mu_im.val[1], m.val[1]);

        vfmul(d_re.val[0], mu_re.val[0], g01_re.val[0]);
        vfmul(d_re.val[1], mu_re.val[1], g01_re.val[1]);
        vfmul(d_im.val[0], mu_im.val[0], g01_re.val[0]);
        vfmul(d_im.val[1], mu_im.val[1], g01_re.val[1]);

        vfma(d_re.val[0], d_re.val[0], mu_im.val[0], g01_im.val[0]);
        vfma(d_re.val[1], d_re.val[1], mu_im.val[1], g01_im.val[1]);
        vfms(d_im.val[0], d_im.val[0], mu_re.val[0], g01_im.val[0]);
        vfms(d_im.val[1], d_im.val[1], mu_re.val[1], g01_im.val[1]);

        vfsub(g11_re.val[0], g11_re.val[0], d_re.val[0]);
        vfsub(g11_re.val[1], g11_re.val[1], d_re.val[1]);
        vfsub(g11_re.val[2], g11_re.val[2], d_im.val[0]);
        vfsub(g11_re.val[3], g11_re.val[3], d_im.val[1]);

        vfneg(mu_re.val[2], mu_im.val[0]);
        vfneg(mu_re.val[3], mu_im.val[1]);

        vstorex4(&g11[0], g11_re);
        vstorex4(&g01[0], mu_re);

        break;

    default:
        // assert(logn >= 4);
        for (int i = 0; i < hn; i += 8)
        {
            vloadx4(g00_re, &g00[i]);
            vloadx4(g00_im, &g00[i + hn]);
            vloadx4(g01_re, &g01[i]);
            vloadx4(g01_im, &g01[i + hn]);
            vloadx4(g11_re, &g11[i]);
            vloadx4(g11_im, &g11[i + hn]);

            vfmulx4(m, g00_re, g00_re);
            vfmax4(m, m, g00_im, g00_im);
            vfinvx4(m, m);

            vfmulx4(mu_re, g01_re, g00_re);
            vfmulx4(mu_im, g01_im, g00_re);

            vfmax4(mu_re, mu_re, g01_im, g00_im);
            vfmsx4(mu_im, mu_im, g01_re, g00_im);

            vfmulx4(mu_re, mu_re, m);
            vfmulx4(mu_im, mu_im, m);

            vfmulx4(d_re, mu_re, g01_re);
            vfmulx4(d_im, mu_im, g01_re);

            vfmax4(d_re, d_re, mu_im, g01_im);
            vfmsx4(d_im, d_im, mu_re, g01_im);

            vfsubx4(g11_re, g11_re, d_re);
            vfsubx4(g11_im, g11_im, d_im);

            vstorex4(&g11[i], g11_re);
            vstorex4(&g11[i + hn], g11_im);

            vfnegx4(mu_im, mu_im);

            vstorex4(&g01[i], mu_re);
            vstorex4(&g01[i + hn], mu_im);
        }
        break;
    }
}

void ZfN(poly_LDLmv_fft)(fpr *restrict d11, fpr *restrict l10,
                         const fpr *restrict g00, const fpr *restrict g01,
                         const fpr *restrict g11, unsigned logn)
{
    const int falcon_n = 1 << logn;
    const int hn = falcon_n >> 1;
    float64x2x4_t g00_re, g00_im, g01_re, g01_im, g11_re, g11_im;
    float64x2x4_t mu_re, mu_im, m, d_re, d_im;
    float64x2x2_t tmp;
    float64x2_t neon_1i2;

    const fpr imagine[2] = {1.0, -1.0};

    switch (logn)
    {
    case 1:
        // n = 2; hn = 1;
        vload(g00_re.val[0], &g00[0]);
        vload(g01_re.val[0], &g01[0]);
        vload(g11_re.val[0], &g11[0]);
        vload(neon_1i2, &imagine[0]);

        // g00_re^2 | g00_im^2
        vfmul(m.val[0], g00_re.val[0], g00_re.val[0]);
        // 1 / ( g00_re^2 + g00_im^2 )
        m.val[0] = vdupq_n_f64(1 / vaddvq_f64(m.val[0]));

        // re: g01_re * g00_re + g01_im * g00_im
        // im: g01_im * g00_re - g01_re * g00_im
        vfmul_lane(mu_re.val[0], g01_re.val[0], g00_re.val[0], 0);
        vfcmla_270(mu_re.val[0], g00_re.val[0], g01_re.val[0]);
        vfmul(mu_re.val[0], mu_re.val[0], m.val[0]);

        vswap(mu_re.val[1], mu_re.val[0]);
        // im: mu_im * g01_re - mu_re * g01_im
        // re: mu_im * g01_im + mu_re * g01_re
        vfmul_lane(g01_re.val[1], g01_re.val[0], mu_re.val[1], 0);
        vfcmla_90(g01_re.val[1], mu_re.val[1], g01_re.val[0]);

        vswap(g01_re.val[0], g01_re.val[1]);

        vfsub(g11_re.val[0], g11_re.val[0], g01_re.val[0]);
        vfmul(mu_re.val[0], mu_re.val[0], neon_1i2);

        vstore(&d11[0], g11_re.val[0]);
        vstore(&l10[0], mu_re.val[0]);

        break;

    case 2:
        // n = 4; hn = 2
        vloadx2(tmp, &g00[0]);
        g00_re.val[0] = tmp.val[0];
        g00_im.val[0] = tmp.val[1];
        vloadx2(tmp, &g01[0]);
        g01_re.val[0] = tmp.val[0];
        g01_im.val[0] = tmp.val[1];
        vloadx2(tmp, &g11[0]);
        g11_re.val[0] = tmp.val[0];
        g11_im.val[0] = tmp.val[1];

        // mu = g01 / g00
        vfmul(m.val[0], g00_re.val[0], g00_re.val[0]);
        vfma(m.val[0], m.val[0], g00_im.val[0], g00_im.val[0]);
        vfinv(m.val[0], m.val[0]);

        vfmul(mu_re.val[0], g01_re.val[0], g00_re.val[0]);
        vfmul(mu_im.val[0], g01_im.val[0], g00_re.val[0]);

        vfma(mu_re.val[0], mu_re.val[0], g01_im.val[0], g00_im.val[0]);
        vfms(mu_im.val[0], mu_im.val[0], g01_re.val[0], g00_im.val[0]);

        vfmul(mu_re.val[0], mu_re.val[0], m.val[0]);
        vfmul(mu_im.val[0], mu_im.val[0], m.val[0]);
        // d = g01 = mu*g01
        vfmul(d_re.val[0], mu_re.val[0], g01_re.val[0]);
        vfmul(d_im.val[0], mu_im.val[0], g01_re.val[0]);

        vfma(d_re.val[0], d_re.val[0], mu_im.val[0], g01_im.val[0]);
        vfms(d_im.val[0], d_im.val[0], mu_re.val[0], g01_im.val[0]);
        // d11 = g11 - g01 = g11 - d
        vfsub(g11_re.val[0], g11_re.val[0], d_re.val[0]);
        vfsub(g11_im.val[0], g11_im.val[0], d_im.val[0]);

        tmp.val[0] = g11_re.val[0];
        tmp.val[1] = g11_im.val[0];
        vstorex2(&d11[0], tmp);

        vfneg(mu_im.val[0], mu_im.val[0]);
        tmp.val[0] = mu_re.val[0];
        tmp.val[1] = mu_im.val[0];
        vstorex2(&l10[0], tmp);

        break;
    case 3:
        //  n = 8; hn = 4
        vloadx4(g00_re, &g00[0]);
        vloadx4(g01_re, &g01[0]);
        vloadx4(g11_re, &g11[0]);
        g00_im.val[0] = g00_re.val[2];
        g00_im.val[1] = g00_re.val[3];
        g01_im.val[0] = g01_re.val[2];
        g01_im.val[1] = g01_re.val[3];

        vfmul(m.val[0], g00_re.val[0], g00_re.val[0]);
        vfmul(m.val[1], g00_re.val[1], g00_re.val[1]);
        vfma(m.val[0], m.val[0], g00_im.val[0], g00_im.val[0]);
        vfma(m.val[1], m.val[1], g00_im.val[1], g00_im.val[1]);
        vfinv(m.val[0], m.val[0]);
        vfinv(m.val[1], m.val[1]);

        vfmul(mu_re.val[0], g01_re.val[0], g00_re.val[0]);
        vfmul(mu_re.val[1], g01_re.val[1], g00_re.val[1]);
        vfmul(mu_im.val[0], g01_im.val[0], g00_re.val[0]);
        vfmul(mu_im.val[1], g01_im.val[1], g00_re.val[1]);

        vfma(mu_re.val[0], mu_re.val[0], g01_im.val[0], g00_im.val[0]);
        vfma(mu_re.val[1], mu_re.val[1], g01_im.val[1], g00_im.val[1]);
        vfms(mu_im.val[0], mu_im.val[0], g01_re.val[0], g00_im.val[0]);
        vfms(mu_im.val[1], mu_im.val[1], g01_re.val[1], g00_im.val[1]);

        vfmul(mu_re.val[0], mu_re.val[0], m.val[0]);
        vfmul(mu_re.val[1], mu_re.val[1], m.val[1]);
        vfmul(mu_im.val[0], mu_im.val[0], m.val[0]);
        vfmul(mu_im.val[1], mu_im.val[1], m.val[1]);

        vfmul(d_re.val[0], mu_re.val[0], g01_re.val[0]);
        vfmul(d_re.val[1], mu_re.val[1], g01_re.val[1]);
        vfmul(d_im.val[0], mu_im.val[0], g01_re.val[0]);
        vfmul(d_im.val[1], mu_im.val[1], g01_re.val[1]);

        vfma(d_re.val[0], d_re.val[0], mu_im.val[0], g01_im.val[0]);
        vfma(d_re.val[1], d_re.val[1], mu_im.val[1], g01_im.val[1]);
        vfms(d_im.val[0], d_im.val[0], mu_re.val[0], g01_im.val[0]);
        vfms(d_im.val[1], d_im.val[1], mu_re.val[1], g01_im.val[1]);

        vfsub(g11_re.val[0], g11_re.val[0], d_re.val[0]);
        vfsub(g11_re.val[1], g11_re.val[1], d_re.val[1]);
        vfsub(g11_re.val[2], g11_re.val[2], d_im.val[0]);
        vfsub(g11_re.val[3], g11_re.val[3], d_im.val[1]);

        vfneg(mu_re.val[2], mu_im.val[0]);
        vfneg(mu_re.val[3], mu_im.val[1]);

        vstorex4(&d11[0], g11_re);
        vstorex4(&l10[0], mu_re);

        break;

    default:
        for (int i = 0; i < hn; i += 8)
        {
            vloadx4(g00_re, &g00[i]);
            vloadx4(g00_im, &g00[i + hn]);
            vloadx4(g01_re, &g01[i]);
            vloadx4(g01_im, &g01[i + hn]);
            vloadx4(g11_re, &g11[i]);
            vloadx4(g11_im, &g11[i + hn]);

            vfmulx4(m, g00_re, g00_re);
            vfmax4(m, m, g00_im, g00_im);
            vfinvx4(m, m);

            vfmulx4(mu_re, g01_re, g00_re);
            vfmulx4(mu_im, g01_im, g00_re);

            vfmax4(mu_re, mu_re, g01_im, g00_im);
            vfmsx4(mu_im, mu_im, g01_re, g00_im);

            vfmulx4(mu_re, mu_re, m);
            vfmulx4(mu_im, mu_im, m);

            vfmulx4(d_re, mu_re, g01_re);
            vfmulx4(d_im, mu_im, g01_re);

            vfmax4(d_re, d_re, mu_im, g01_im);
            vfmsx4(d_im, d_im, mu_re, g01_im);

            vfsubx4(g11_re, g11_re, d_re);
            vfsubx4(g11_im, g11_im, d_im);

            vstorex4(&d11[i], g11_re);
            vstorex4(&d11[i + hn], g11_im);

            vfnegx4(mu_im, mu_im);

            vstorex4(&l10[i], mu_re);
            vstorex4(&l10[i + hn], mu_im);
        }
        break;
    }
}

void ZfN(poly_fpr_of_s16)(fpr *t0, const uint16_t *hm, const unsigned falcon_n)
{
    float64x2x4_t neon_t0;
    uint16x8x4_t neon_hm;
    uint16x8_t neon_zero;
    uint32x4x4_t neon_hmu32[2];
    int64x2x4_t neon_hms64[4];
    neon_zero = vdupq_n_u16(0);
    for (unsigned u = 0; u < falcon_n; u += 32)
    {
        neon_hm = vld1q_u16_x4(&hm[u]);
        neon_hmu32[0].val[0] = (uint32x4_t)vzip1q_u16(neon_hm.val[0], neon_zero);
        neon_hmu32[0].val[1] = (uint32x4_t)vzip2q_u16(neon_hm.val[0], neon_zero);
        neon_hmu32[0].val[2] = (uint32x4_t)vzip1q_u16(neon_hm.val[1], neon_zero);
        neon_hmu32[0].val[3] = (uint32x4_t)vzip2q_u16(neon_hm.val[1], neon_zero);

        neon_hmu32[1].val[0] = (uint32x4_t)vzip1q_u16(neon_hm.val[2], neon_zero);
        neon_hmu32[1].val[1] = (uint32x4_t)vzip2q_u16(neon_hm.val[2], neon_zero);
        neon_hmu32[1].val[2] = (uint32x4_t)vzip1q_u16(neon_hm.val[3], neon_zero);
        neon_hmu32[1].val[3] = (uint32x4_t)vzip2q_u16(neon_hm.val[3], neon_zero);

        neon_hms64[0].val[0] = (int64x2_t)vzip1q_u32(neon_hmu32[0].val[0], (uint32x4_t)neon_zero);
        neon_hms64[0].val[1] = (int64x2_t)vzip2q_u32(neon_hmu32[0].val[0], (uint32x4_t)neon_zero);
        neon_hms64[0].val[2] = (int64x2_t)vzip1q_u32(neon_hmu32[0].val[1], (uint32x4_t)neon_zero);
        neon_hms64[0].val[3] = (int64x2_t)vzip2q_u32(neon_hmu32[0].val[1], (uint32x4_t)neon_zero);

        neon_hms64[1].val[0] = (int64x2_t)vzip1q_u32(neon_hmu32[0].val[2], (uint32x4_t)neon_zero);
        neon_hms64[1].val[1] = (int64x2_t)vzip2q_u32(neon_hmu32[0].val[2], (uint32x4_t)neon_zero);
        neon_hms64[1].val[2] = (int64x2_t)vzip1q_u32(neon_hmu32[0].val[3], (uint32x4_t)neon_zero);
        neon_hms64[1].val[3] = (int64x2_t)vzip2q_u32(neon_hmu32[0].val[3], (uint32x4_t)neon_zero);

        neon_hms64[2].val[0] = (int64x2_t)vzip1q_u32(neon_hmu32[1].val[0], (uint32x4_t)neon_zero);
        neon_hms64[2].val[1] = (int64x2_t)vzip2q_u32(neon_hmu32[1].val[0], (uint32x4_t)neon_zero);
        neon_hms64[2].val[2] = (int64x2_t)vzip1q_u32(neon_hmu32[1].val[1], (uint32x4_t)neon_zero);
        neon_hms64[2].val[3] = (int64x2_t)vzip2q_u32(neon_hmu32[1].val[1], (uint32x4_t)neon_zero);

        neon_hms64[3].val[0] = (int64x2_t)vzip1q_u32(neon_hmu32[1].val[2], (uint32x4_t)neon_zero);
        neon_hms64[3].val[1] = (int64x2_t)vzip2q_u32(neon_hmu32[1].val[2], (uint32x4_t)neon_zero);
        neon_hms64[3].val[2] = (int64x2_t)vzip1q_u32(neon_hmu32[1].val[3], (uint32x4_t)neon_zero);
        neon_hms64[3].val[3] = (int64x2_t)vzip2q_u32(neon_hmu32[1].val[3], (uint32x4_t)neon_zero);

        neon_t0.val[0] = vcvtq_f64_s64(neon_hms64[0].val[0]);
        neon_t0.val[1] = vcvtq_f64_s64(neon_hms64[0].val[1]);
        neon_t0.val[2] = vcvtq_f64_s64(neon_hms64[0].val[2]);
        neon_t0.val[3] = vcvtq_f64_s64(neon_hms64[0].val[3]);
        vstorex4(&t0[u], neon_t0);

        neon_t0.val[0] = vcvtq_f64_s64(neon_hms64[1].val[0]);
        neon_t0.val[1] = vcvtq_f64_s64(neon_hms64[1].val[1]);
        neon_t0.val[2] = vcvtq_f64_s64(neon_hms64[1].val[2]);
        neon_t0.val[3] = vcvtq_f64_s64(neon_hms64[1].val[3]);
        vstorex4(&t0[u + 8], neon_t0);

        neon_t0.val[0] = vcvtq_f64_s64(neon_hms64[2].val[0]);
        neon_t0.val[1] = vcvtq_f64_s64(neon_hms64[2].val[1]);
        neon_t0.val[2] = vcvtq_f64_s64(neon_hms64[2].val[2]);
        neon_t0.val[3] = vcvtq_f64_s64(neon_hms64[2].val[3]);
        vstorex4(&t0[u + 16], neon_t0);

        neon_t0.val[0] = vcvtq_f64_s64(neon_hms64[3].val[0]);
        neon_t0.val[1] = vcvtq_f64_s64(neon_hms64[3].val[1]);
        neon_t0.val[2] = vcvtq_f64_s64(neon_hms64[3].val[2]);
        neon_t0.val[3] = vcvtq_f64_s64(neon_hms64[3].val[3]);
        vstorex4(&t0[u + 24], neon_t0);
    }
}

fpr ZfN(compute_bnorm)(const fpr *rt1, const fpr *rt2)
{
    float64x2x4_t r1, r11, r2, r22;
    float64x2_t bnorm;

    bnorm = vdupq_n_f64(0);

    for (int i = 0; i < FALCON_N; i += 16)
    {
        vloadx4(r1, &rt1[i]);
        vloadx4(r11, &rt1[i + 8]);

        vfma(bnorm, bnorm, r1.val[0], r1.val[0]);
        vfma(bnorm, bnorm, r1.val[1], r1.val[1]);
        vfma(bnorm, bnorm, r1.val[2], r1.val[2]);
        vfma(bnorm, bnorm, r1.val[3], r1.val[3]);

        vfma(bnorm, bnorm, r11.val[0], r11.val[0]);
        vfma(bnorm, bnorm, r11.val[1], r11.val[1]);
        vfma(bnorm, bnorm, r11.val[2], r11.val[2]);
        vfma(bnorm, bnorm, r11.val[3], r11.val[3]);
    }

    for (int i = 0; i < FALCON_N; i += 16)
    {
        vloadx4(r2, &rt2[i]);
        vloadx4(r22, &rt2[i + 8]);

        vfma(bnorm, bnorm, r2.val[0], r2.val[0]);
        vfma(bnorm, bnorm, r2.val[1], r2.val[1]);
        vfma(bnorm, bnorm, r2.val[2], r2.val[2]);
        vfma(bnorm, bnorm, r2.val[3], r2.val[3]);

        vfma(bnorm, bnorm, r22.val[0], r22.val[0]);
        vfma(bnorm, bnorm, r22.val[1], r22.val[1]);
        vfma(bnorm, bnorm, r22.val[2], r22.val[2]);
        vfma(bnorm, bnorm, r22.val[3], r22.val[3]);
    }

    return vaddvq_f64(bnorm);
}
