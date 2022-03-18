/*
 * High-speed vectorize FFT tree for arbitrary `logn`.
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
#include "macrof.h"
#include "macrofx4.h"

/* 
 * Minimum logn: 5
 */
static void ZfN(poly_mergeFFT_log5)(fpr *f, const fpr *f0, const fpr *f1, unsigned logn)
{
    // n = 32; hn = 16; qn = 8
    // Total: 32 register
    float64x2x4_t f0_re, f0_im, f1_re, f1_im, d_re, d_im; // 16
    float64x2x2_t s_tmp[4];                               // 8
    float64x2x2_t x_tmp, y_tmp;                           // 4

    const unsigned falcon_n = 1 << logn;
    const unsigned hn = falcon_n >> 1;
    const unsigned qn = falcon_n >> 2;
    unsigned u1, u2;
    for (unsigned u = 0; u < qn; u += 8)
    {
        u1 = u << 1;
        u2 = u1 + falcon_n;
        vloadx4(f0_re, &f0[u]);
        vloadx4(f0_im, &f0[u + qn]);
        vloadx4(f1_re, &f1[u]);
        vloadx4(f1_im, &f1[u + qn]);

        vload2(s_tmp[0], &fpr_gm_tab[u2 + 0]);
        vload2(s_tmp[1], &fpr_gm_tab[u2 + 4]);
        vload2(s_tmp[2], &fpr_gm_tab[u2 + 8]);
        vload2(s_tmp[3], &fpr_gm_tab[u2 + 12]);

        // f0,f1_re: 0-> 7
        // f0,f1_im: qn -> qn + 7
        vfmul(d_re.val[0], f1_re.val[0], s_tmp[0].val[0]);
        vfmul(d_re.val[1], f1_re.val[1], s_tmp[1].val[0]);
        vfmul(d_re.val[2], f1_re.val[2], s_tmp[2].val[0]);
        vfmul(d_re.val[3], f1_re.val[3], s_tmp[3].val[0]);

        vfms(d_re.val[0], d_re.val[0], f1_im.val[0], s_tmp[0].val[1]);
        vfms(d_re.val[1], d_re.val[1], f1_im.val[1], s_tmp[1].val[1]);
        vfms(d_re.val[2], d_re.val[2], f1_im.val[2], s_tmp[2].val[1]);
        vfms(d_re.val[3], d_re.val[3], f1_im.val[3], s_tmp[3].val[1]);

        // f_re0
        // vfaddx4(v0, f0_re, tmp);
        vfadd(x_tmp.val[0], f0_re.val[0], d_re.val[0]);
        vfsub(x_tmp.val[1], f0_re.val[0], d_re.val[0]);
        vfadd(y_tmp.val[0], f0_re.val[1], d_re.val[1]);
        vfsub(y_tmp.val[1], f0_re.val[1], d_re.val[1]);

        vstore2(&f[u1], x_tmp);
        vstore2(&f[u1 + 4], y_tmp);

        // f_re1
        // vfsubx4(v1, f0_re, tmp);
        vfadd(x_tmp.val[0], f0_re.val[2], d_re.val[2]);
        vfsub(x_tmp.val[1], f0_re.val[2], d_re.val[2]);
        vfadd(y_tmp.val[0], f0_re.val[3], d_re.val[3]);
        vfsub(y_tmp.val[1], f0_re.val[3], d_re.val[3]);

        vstore2(&f[u1 + 8], x_tmp);
        vstore2(&f[u1 + 12], y_tmp);

        vfmul(d_im.val[0], f1_re.val[0], s_tmp[0].val[1]);
        vfmul(d_im.val[1], f1_re.val[1], s_tmp[1].val[1]);
        vfmul(d_im.val[2], f1_re.val[2], s_tmp[2].val[1]);
        vfmul(d_im.val[3], f1_re.val[3], s_tmp[3].val[1]);

        vfma(d_im.val[0], d_im.val[0], f1_im.val[0], s_tmp[0].val[0]);
        vfma(d_im.val[1], d_im.val[1], f1_im.val[1], s_tmp[1].val[0]);
        vfma(d_im.val[2], d_im.val[2], f1_im.val[2], s_tmp[2].val[0]);
        vfma(d_im.val[3], d_im.val[3], f1_im.val[3], s_tmp[3].val[0]);

        // f_re0
        // vfaddx4(v0, f0_im, tmp);
        vfadd(x_tmp.val[0], f0_im.val[0], d_im.val[0]);
        vfsub(x_tmp.val[1], f0_im.val[0], d_im.val[0]);
        vfadd(y_tmp.val[0], f0_im.val[1], d_im.val[1]);
        vfsub(y_tmp.val[1], f0_im.val[1], d_im.val[1]);

        vstore2(&f[u1 + hn], x_tmp);
        vstore2(&f[u1 + hn + 4], y_tmp);

        // f_re1
        // vfsubx4(v1, f0_im, tmp);
        vfadd(x_tmp.val[0], f0_im.val[2], d_im.val[2]);
        vfsub(x_tmp.val[1], f0_im.val[2], d_im.val[2]);
        vfadd(y_tmp.val[0], f0_im.val[3], d_im.val[3]);
        vfsub(y_tmp.val[1], f0_im.val[3], d_im.val[3]);

        vstore2(&f[u1 + hn + 8], x_tmp);
        vstore2(&f[u1 + hn + 12], y_tmp);
    }
}

/* 
 * Fix logn: 4
 */
static inline void ZfN(poly_mergeFFT_log4)(fpr *f, const fpr *f0, const fpr *f1)
{
    // n = 16; hn = 8; qn = 4
    // a_re = f0[0, 1, 2, 3]
    // a_im = f0[4, 5, 6, 7]
    // d_re = f1[0, 1, 2, 3]
    // d_im = f1[4, 5, 6, 7]

    // Total SIMD register: 24 register
    float64x2x4_t v0, v1;                             // 12
    float64x2x2_t s_tmp, s_re, s_im, d_re, d_im, tmp; // 12

    vloadx4(v0, &f0[0]);
    vloadx4(v1, &f1[0]);
    // 16, 18, 20, 22
    // 17, 19, 21, 23
    vload2(s_tmp, &fpr_gm_tab[16]);
    s_re.val[0] = s_tmp.val[0];
    s_im.val[0] = s_tmp.val[1];
    vload2(s_tmp, &fpr_gm_tab[20]);
    s_re.val[1] = s_tmp.val[0];
    s_im.val[1] = s_tmp.val[1];

    vfmul(d_re.val[0], v1.val[0], s_re.val[0]);
    vfmul(d_re.val[1], v1.val[1], s_re.val[1]);
    vfms(d_re.val[0], d_re.val[0], v1.val[2], s_im.val[0]);
    vfms(d_re.val[1], d_re.val[1], v1.val[3], s_im.val[1]);

    vfmul(d_im.val[0], v1.val[0], s_im.val[0]);
    vfmul(d_im.val[1], v1.val[1], s_im.val[1]);
    vfma(d_im.val[0], d_im.val[0], v1.val[2], s_re.val[0]);
    vfma(d_im.val[1], d_im.val[1], v1.val[3], s_re.val[1]);

    vfadd(tmp.val[0], v0.val[0], d_re.val[0]);
    vfsub(tmp.val[1], v0.val[0], d_re.val[0]);
    vstore2(&f[0], tmp);

    vfadd(tmp.val[0], v0.val[1], d_re.val[1]);
    vfsub(tmp.val[1], v0.val[1], d_re.val[1]);
    vstore2(&f[4], tmp);

    vfadd(tmp.val[0], v0.val[2], d_im.val[0]);
    vfsub(tmp.val[1], v0.val[2], d_im.val[0]);
    vstore2(&f[8], tmp);

    vfadd(tmp.val[0], v0.val[3], d_im.val[1]);
    vfsub(tmp.val[1], v0.val[3], d_im.val[1]);
    vstore2(&f[12], tmp);
}

/* 
 * Fix logn: 3
 */
static inline void ZfN(poly_mergeFFT_log3)(fpr *f, const fpr *f0, const fpr *f1)
{
    // n = 8, hn = 4, qn = 2
    // a_re = f0[0, 1]
    // a_im = f0[2, 3]
    // d_re = f1[0, 1]
    // d_im = f1[2, 3]
    // Total: 12 register
    float64x2_t a_re, a_im, d_re, d_im, b_re, b_im, t_re, t_im; // 6
    float64x2x2_t tmp, s_re_im;                                 // 2

    vloadx2(tmp, &f0[0]);
    a_re = tmp.val[0];
    a_im = tmp.val[1];
    vloadx2(tmp, &f1[0]);
    d_re = tmp.val[0];
    d_im = tmp.val[1];
    // 8, 10
    // 9, 11
    vload2(s_re_im, &fpr_gm_tab[8]);

    vfmul(b_re, d_re, s_re_im.val[0]);
    vfms(b_re, b_re, d_im, s_re_im.val[1]);

    vfmul(b_im, d_re, s_re_im.val[1]);
    vfma(b_im, b_im, d_im, s_re_im.val[0]);

    // 0, 2
    // 4, 6
    vfadd(t_re, a_re, b_re);
    vfadd(t_im, a_im, b_im);
    // 1, 3
    // 5, 7
    vfsub(d_re, a_re, b_re);
    vfsub(d_im, a_im, b_im);
    tmp.val[0] = t_re;
    tmp.val[1] = d_re;
    vstore2(&f[0], tmp);
    tmp.val[0] = t_im;
    tmp.val[1] = d_im;
    vstore2(&f[4], tmp);
}

/* 
 * see inner.h
 * Only support logn >= 3
 */
void ZfN(poly_merge_fft)(fpr *restrict f, const fpr *restrict f0,
                        const fpr *restrict f1, unsigned logn)
{
    float64x2_t a_re_im, b_re_im, d_re_im,
        s_re_im, s_im_re, neon_1i2, d_re, d_im;
    float64x2x2_t t_re_im;
    const fpr imagine[2] = {1.0, -1.0};

    switch (logn)
    {
    case 1:
        // n = 2; hn = 1;
        f[0] = f0[0];
        f[1] = f1[0];
        break;

    case 2:
        // n = 4; hn = 2; qn = 1;
        // a_re = f0[0];
        // a_im = f0[1];
        // d_re = f1[0];
        // d_im = f1[1];
        // s_re_im = 4, 5
        vload(a_re_im, &f0[0]);
        vload(neon_1i2, &imagine[0]);
        vload(d_re_im, &f1[0]);
        vload(s_re_im, &fpr_gm_tab[4]);
        s_im_re = vextq_f64(s_re_im, s_re_im, 1);
        vfmul(s_re_im, s_re_im, neon_1i2);

        vfmul(d_re, d_re_im, s_re_im);
        vfmul(d_im, d_re_im, s_im_re);

        b_re_im = vpaddq_f64(d_re, d_im);
        vfadd(t_re_im.val[0], a_re_im, b_re_im);
        vfsub(t_re_im.val[1], a_re_im, b_re_im);

        vstore2(&f[0], t_re_im);
        break;

    case 3:
        ZfN(poly_mergeFFT_log3)(f, f0, f1);
        break;

    case 4:
        ZfN(poly_mergeFFT_log4)(f, f0, f1);
        break;

    default:
        ZfN(poly_mergeFFT_log5)(f, f0, f1, logn);
        break;
    }
}

static inline void ZfN(poly_splitFFT_log3)(fpr *restrict f0, fpr *restrict f1, const fpr *f)
{
    // n = 8; hn = 4; qn = 2;
    // a_re = f[0];
    // b_re = f[1];
    // a_im = f[4];
    // b_im = f[5];

    // a_re = f[2];
    // b_re = f[3];
    // a_im = f[6];
    // b_im = f[7];

    // Total SIMD register: 10 = 8 + (4-2)
    float64x2x2_t tmp, s_re_im;                                 // 4
    float64x2_t a_re, b_re, a_im, b_im, t_re, t_im, d_re, d_im; // 8

    vload2(tmp, &f[0]);
    a_re = tmp.val[0];
    b_re = tmp.val[1];
    vload2(tmp, &f[4]);
    a_im = tmp.val[0];
    b_im = tmp.val[1];

    // 8, 10
    // 9, 11
    vload2(s_re_im, &fpr_gm_tab_half[8]);

    vfadd(t_re, a_re, b_re);
    vfadd(t_im, a_im, b_im);
    vfmuln(t_re, t_re, 0.5);
    vfmuln(t_im, t_im, 0.5);

    // 0, 2, 1, 3
    tmp.val[0] = t_re;
    tmp.val[1] = t_im;
    vstorex2(&f0[0], tmp);

    vfsub(t_re, a_re, b_re);
    vfsub(t_im, a_im, b_im);

    vfmul(d_re, t_re, s_re_im.val[0]);
    vfma(d_re, d_re, t_im, s_re_im.val[1]);

    vfmul(d_im, t_im, s_re_im.val[0]);
    vfms(d_im, d_im, t_re, s_re_im.val[1]);

    tmp.val[0] = d_re;
    tmp.val[1] = d_im;
    vstorex2(&f1[0], tmp);
}

static inline void ZfN(poly_splitFFT_log4)(fpr *restrict f0, fpr *restrict f1,
                                          const fpr *f, unsigned logn)
{
    // n = 16; hn = 8; qn = 4
    // a_re = f[0,  4,  2,  6];
    // b_re = f[1,  5,  3,  7];
    // a_im = f[8, 12, 10, 14];
    // b_im = f[9, 13, 11, 15];

    // Total SIMD register: 20
    float64x2x2_t a_re, b_re, a_im, b_im, s_re, s_im, t_re, t_im, d_re, d_im;
    float64x2x4_t tmp, s_re_im;

    const unsigned falcon_n = 1 << logn;
    const unsigned hn = falcon_n >> 1;
    const unsigned qn = falcon_n >> 2;
    unsigned u1;

    for (unsigned u = 0; u < qn; u += 4)
    {
        u1 = u << 1;
        vload4(tmp, &f[u1]);
        a_re.val[0] = tmp.val[0];
        b_re.val[0] = tmp.val[1];
        a_re.val[1] = tmp.val[2];
        b_re.val[1] = tmp.val[3];
        vload4(tmp, &f[u1 + hn]);
        a_im.val[0] = tmp.val[0];
        b_im.val[0] = tmp.val[1];
        a_im.val[1] = tmp.val[2];
        b_im.val[1] = tmp.val[3];
        // s_re = 16, 20, 18, 22
        // s_im = 17, 21, 19, 23
        vload4(s_re_im, &fpr_gm_tab_half[u1 + falcon_n]);
        s_re.val[0] = s_re_im.val[0];
        s_re.val[1] = s_re_im.val[2];
        s_im.val[0] = s_re_im.val[1];
        s_im.val[1] = s_re_im.val[3];

        vfadd(t_re.val[0], a_re.val[0], b_re.val[0]);
        vfadd(t_im.val[0], a_im.val[0], b_im.val[0]);
        vfadd(t_re.val[1], a_re.val[1], b_re.val[1]);
        vfadd(t_im.val[1], a_im.val[1], b_im.val[1]);

        // 0, 4
        vfmuln(t_re.val[0], t_re.val[0], 0.5);
        // 8, 12
        vfmuln(t_im.val[0], t_im.val[0], 0.5);
        // 2, 6
        vfmuln(t_re.val[1], t_re.val[1], 0.5);
        // 10, 14
        vfmuln(t_im.val[1], t_im.val[1], 0.5);

        vstore2(&f0[u], t_re);
        vstore2(&f0[u + qn], t_im);

        vfsub(t_re.val[0], a_re.val[0], b_re.val[0]);
        vfsub(t_im.val[0], a_im.val[0], b_im.val[0]);
        vfsub(t_re.val[1], a_re.val[1], b_re.val[1]);
        vfsub(t_im.val[1], a_im.val[1], b_im.val[1]);

        vfmul(d_re.val[0], t_re.val[0], s_re.val[0]);
        vfmul(d_re.val[1], t_re.val[1], s_re.val[1]);
        vfma(d_re.val[0], d_re.val[0], t_im.val[0], s_im.val[0]);
        vfma(d_re.val[1], d_re.val[1], t_im.val[1], s_im.val[1]);

        vfmul(d_im.val[0], t_im.val[0], s_re.val[0]);
        vfmul(d_im.val[1], t_im.val[1], s_re.val[1]);
        vfms(d_im.val[0], d_im.val[0], t_re.val[0], s_im.val[0]);
        vfms(d_im.val[1], d_im.val[1], t_re.val[1], s_im.val[1]);

        vstore2(&f1[u], d_re);
        vstore2(&f1[u + qn], d_im);
    }
}

/* 
 * Recursive Split FFT
 */
void ZfN(poly_split_fft)(fpr *restrict f0, fpr *restrict f1,
                        const fpr *restrict f, unsigned logn)
{
    float64x2x2_t tmp;
    float64x2_t a_re_im, b_re_im, t_re_im, d_re, d_im, neon_1i2;
    const fpr imagine[2] = {1.0, -1.0};

    switch (logn)
    {
    case 1:
        //  n = 2; hn = 1; qn = 0;
        f0[0] = f[0];
        f1[0] = f[1];
        break;
    case 2:
        // n = 4; hn = 2; qn = 1;
        // a_re = f[0];
        // a_im = f[2];
        // b_re = f[1];
        // b_im = f[3];
        vload2(tmp, &f[0]);
        a_re_im = tmp.val[0];
        b_re_im = tmp.val[1];
        vfadd(t_re_im, a_re_im, b_re_im);
        vfmuln(t_re_im, t_re_im, 0.5);
        vstore(&f0[0], t_re_im);

        vload(tmp.val[0], &fpr_gm_tab_half[4]);
        vload(neon_1i2, &imagine[0]);
        vfsub(t_re_im, a_re_im, b_re_im);

        vfmul(tmp.val[1], tmp.val[0], neon_1i2);
        tmp.val[1] = vextq_f64(tmp.val[1], tmp.val[1], 1);

        vfmul(d_re, t_re_im, tmp.val[0]);
        vfmul(d_im, t_re_im, tmp.val[1]);
        t_re_im = vpaddq_f64(d_re, d_im);

        vstore(&f1[0], t_re_im);
        break;

    case 3:
        ZfN(poly_splitFFT_log3)(f0, f1, f);
        break;

    default:
        ZfN(poly_splitFFT_log4)(f0, f1, f, logn);
        break;
    }
}
