/*
 * High-speed vectorize FFT code for arbitrary `logn`.
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

static void ZfN(FFT_log2)(fpr *f)
{
    /*
    x_re:   0 =   0 + (  1*  4 -   3*  5)
    x_im:   2 =   2 + (  1*  5 +   3*  4)
    y_re:   1 =   0 - (  1*  4 -   3*  5)
    y_im:   3 =   2 - (  1*  5 +   3*  4)
    
    Turn out this vectorize code is too short to be executed,
    the scalar version is consistently faster

    float64x2x2_t t1, t2;
    float64x2_t s_re_im, v;

    re: 0, 2
    im: 1, 3
    vload2(t1, &f[0]);
    vload(s_re_im, &fpr_gm_tab[4]);

    vfmul_lane(v, s_re_im, t1.val[1], 0);
    vfcmla_90(v, t1.val[1], s_re_im);

    vfadd(t2.val[0], t1.val[0], v);
    vfsub(t2.val[1], t1.val[0], v);

    vstore2(&f[0], t2);
    */

    fpr x_re, x_im, y_re, y_im, v_re, v_im, s_re, s_im;

    x_re = f[0];
    x_im = f[2];
    y_re = f[1];
    y_im = f[3];
    s_re = fpr_gm_tab[4];
    s_im = fpr_gm_tab[5];
    
    v_re = y_re*s_re - y_im*s_im;
    v_im = y_re*s_im + y_im*s_re;

    f[0] = x_re + v_re;
    f[2] = x_im + v_im;
    f[1] = x_re - v_re;
    f[3] = x_im - v_im;
}

static void ZfN(FFT_log3)(fpr *f)
{
    float64x2x4_t tmp;
    float64x2x2_t s_re_im, tmp2_0, tmp2_1;
    float64x2_t v_re, v_im, x_re, x_im, y_re, y_im;

    // 0: 0, 1
    // 1: 2, 3
    // 2: 4, 5
    // 3: 6, 7
    vloadx4(tmp, &f[0]);
    s_re_im.val[0] = vld1q_f64(&fpr_gm_tab[4]);

    /*
    level 1
    x_re:   0 =   0 + (  2*  4 -   6*  5)
    x_re:   1 =   1 + (  3*  4 -   7*  5)
    y_re:   2 =   0 - (  2*  4 -   6*  5)
    y_re:   3 =   1 - (  3*  4 -   7*  5)

    x_im:   4 =   4 + (  2*  5 +   6*  4)
    x_im:   5 =   5 + (  3*  5 +   7*  4)
    y_im:   6 =   4 - (  2*  5 +   6*  4)
    y_im:   7 =   5 - (  3*  5 +   7*  4)
    */

    vfmul_lane(v_re, tmp.val[1], s_re_im.val[0], 0);
    vfmul_lane(v_im, tmp.val[1], s_re_im.val[0], 1);

    vfms_lane(v_re, v_re, tmp.val[3], s_re_im.val[0], 1);
    vfma_lane(v_im, v_im, tmp.val[3], s_re_im.val[0], 0);

    vfsub(tmp.val[1], tmp.val[0], v_re);
    vfsub(tmp.val[3], tmp.val[2], v_im);

    vfadd(tmp.val[0], tmp.val[0], v_re);
    vfadd(tmp.val[2], tmp.val[2], v_im);

    /*
    x_re: 0, 2
    y_re: 1, 3
    x_im: 4, 6
    y_im: 5, 7
    */
    x_re = vtrn1q_f64(tmp.val[0], tmp.val[1]);
    y_re = vtrn2q_f64(tmp.val[0], tmp.val[1]);
    x_im = vtrn1q_f64(tmp.val[2], tmp.val[3]);
    y_im = vtrn2q_f64(tmp.val[2], tmp.val[3]);

    /*
    ----
    level 2
    x_re:   0 =   0 + (  1*  8 -   5*  9)
    x_re:   2 =   2 + (  3* 10 -   7* 11)
    y_re:   1 =   0 - (  1*  8 -   5*  9)
    y_re:   3 =   2 - (  3* 10 -   7* 11)

    x_im:   4 =   4 + (  1*  9 +   5*  8)
    x_im:   6 =   6 + (  3* 11 +   7* 10)
    y_im:   5 =   4 - (  1*  9 +   5*  8)
    y_im:   7 =   6 - (  3* 11 +   7* 10)
    */
    // 8, 10
    // 9, 11
    vload2(s_re_im, &fpr_gm_tab[8]);

    vfmul(v_re, y_re, s_re_im.val[0]);
    vfms(v_re, v_re, y_im, s_re_im.val[1]);

    vfmul(v_im, y_re, s_re_im.val[1]);
    vfma(v_im, v_im, y_im, s_re_im.val[0]);

    vfadd(tmp2_0.val[0], x_re, v_re);
    vfsub(tmp2_0.val[1], x_re, v_re);

    vfadd(tmp2_1.val[0], x_im, v_im);
    vfsub(tmp2_1.val[1], x_im, v_im);

    vstore2(&f[0], tmp2_0);
    vstore2(&f[4], tmp2_1);
}

static void ZfN(FFT_log4)(fpr *f)
{
    // Total SIMD register: 28 = 12 + 16
    float64x2x4_t t0, t1, s_re_im;                       // 12
    float64x2x2_t s_tmp, x_re, x_im, y_re, y_im, v1, v2; // 16

    /*
    level 1
    x_re:   0 =   0 + (  4*  4 -  12*  5)
    x_re:   1 =   1 + (  5*  4 -  13*  5)
    x_re:   2 =   2 + (  6*  4 -  14*  5)
    x_re:   3 =   3 + (  7*  4 -  15*  5)

    y_re:   4 =   0 - (  4*  4 -  12*  5)
    y_re:   5 =   1 - (  5*  4 -  13*  5)
    y_re:   6 =   2 - (  6*  4 -  14*  5)
    y_re:   7 =   3 - (  7*  4 -  15*  5)

    x_im:   8 =   8 + (  4*  5 +  12*  4)
    x_im:   9 =   9 + (  5*  5 +  13*  4)
    x_im:  10 =  10 + (  6*  5 +  14*  4)
    x_im:  11 =  11 + (  7*  5 +  15*  4)

    y_im:  12 =   8 - (  4*  5 +  12*  4)
    y_im:  13 =   9 - (  5*  5 +  13*  4)
    y_im:  14 =  10 - (  6*  5 +  14*  4)
    y_im:  15 =  11 - (  7*  5 +  15*  4)
     */

    vloadx4(t0, &f[0]);
    vloadx4(t1, &f[8]);
    s_re_im.val[0] = vld1q_f64(&fpr_gm_tab[4]);

    vfmul_lane(v1.val[0], t0.val[2], s_re_im.val[0], 0);
    vfmul_lane(v1.val[1], t0.val[3], s_re_im.val[0], 0);

    vfmul_lane(v2.val[0], t0.val[2], s_re_im.val[0], 1);
    vfmul_lane(v2.val[1], t0.val[3], s_re_im.val[0], 1);

    vfms_lane(v1.val[0], v1.val[0], t1.val[2], s_re_im.val[0], 1);
    vfms_lane(v1.val[1], v1.val[1], t1.val[3], s_re_im.val[0], 1);

    vfma_lane(v2.val[0], v2.val[0], t1.val[2], s_re_im.val[0], 0);
    vfma_lane(v2.val[1], v2.val[1], t1.val[3], s_re_im.val[0], 0);

    vfsub(t0.val[2], t0.val[0], v1.val[0]);
    vfsub(t0.val[3], t0.val[1], v1.val[1]);

    vfsub(t1.val[2], t1.val[0], v2.val[0]);
    vfsub(t1.val[3], t1.val[1], v2.val[1]);

    vfadd(t0.val[0], t0.val[0], v1.val[0]);
    vfadd(t0.val[1], t0.val[1], v1.val[1]);

    vfadd(t1.val[0], t1.val[0], v2.val[0]);
    vfadd(t1.val[1], t1.val[1], v2.val[1]);

    /*
    Level 2
    y_re:   2 =   0 - (  2*  8 -  10*  9)
    y_re:   3 =   1 - (  3*  8 -  11*  9)
    y_re:   6 =   4 - (  6* 10 -  14* 11)
    y_re:   7 =   5 - (  7* 10 -  15* 11)

    x_re:   0 =   0 + (  2*  8 -  10*  9)
    x_re:   1 =   1 + (  3*  8 -  11*  9)
    x_re:   4 =   4 + (  6* 10 -  14* 11)
    x_re:   5 =   5 + (  7* 10 -  15* 11)


    y_im:  10 =   8 - (  2*  9 +  10*  8)
    y_im:  11 =   9 - (  3*  9 +  11*  8)
    y_im:  14 =  12 - (  6* 11 +  14* 10)
    y_im:  15 =  13 - (  7* 11 +  15* 10)

    x_im:   8 =   8 + (  2*  9 +  10*  8)
    x_im:   9 =   9 + (  3*  9 +  11*  8)
    x_im:  12 =  12 + (  6* 11 +  14* 10)
    x_im:  13 =  13 + (  7* 11 +  15* 10)
     */
    vloadx2(s_tmp, &fpr_gm_tab[8]);

    vfmul_lane(v1.val[0], t0.val[1], s_tmp.val[0], 0);
    vfmul_lane(v1.val[1], t0.val[3], s_tmp.val[1], 0);
    vfmul_lane(v2.val[0], t0.val[1], s_tmp.val[0], 1);
    vfmul_lane(v2.val[1], t0.val[3], s_tmp.val[1], 1);

    vfms_lane(v1.val[0], v1.val[0], t1.val[1], s_tmp.val[0], 1);
    vfms_lane(v1.val[1], v1.val[1], t1.val[3], s_tmp.val[1], 1);
    vfma_lane(v2.val[0], v2.val[0], t1.val[1], s_tmp.val[0], 0);
    vfma_lane(v2.val[1], v2.val[1], t1.val[3], s_tmp.val[1], 0);

    vfsub(y_re.val[0], t0.val[0], v1.val[0]);
    vfsub(y_re.val[1], t0.val[2], v1.val[1]);
    vfsub(y_im.val[0], t1.val[0], v2.val[0]);
    vfsub(y_im.val[1], t1.val[2], v2.val[1]);

    vfadd(x_re.val[0], t0.val[0], v1.val[0]);
    vfadd(x_re.val[1], t0.val[2], v1.val[1]);
    vfadd(x_im.val[0], t1.val[0], v2.val[0]);
    vfadd(x_im.val[1], t1.val[2], v2.val[1]);

    // x_re: 0, 1 | 4, 5
    // y_re: 2, 3 | 6, 7
    // x_im: 8, 9 | 12, 13
    // y_im: 10, 11 | 14, 15

    t0.val[0] = vzip1q_f64(x_re.val[0], x_re.val[1]);
    t0.val[1] = vzip2q_f64(x_re.val[0], x_re.val[1]);
    t0.val[2] = vzip1q_f64(y_re.val[0], y_re.val[1]);
    t0.val[3] = vzip2q_f64(y_re.val[0], y_re.val[1]);

    t1.val[0] = vzip1q_f64(x_im.val[0], x_im.val[1]);
    t1.val[1] = vzip2q_f64(x_im.val[0], x_im.val[1]);
    t1.val[2] = vzip1q_f64(y_im.val[0], y_im.val[1]);
    t1.val[3] = vzip2q_f64(y_im.val[0], y_im.val[1]);

    // t0: 0, 4  | 1, 5  | 2, 6   | 3, 7
    // t1: 8, 12 | 9, 13 | 10, 14 | 11, 15
    /*
    level 3
    y_re:   1 =   0 - (  1* 16 -   9* 17)
    y_re:   5 =   4 - (  5* 20 -  13* 21)
    y_re:   3 =   2 - (  3* 18 -  11* 19)
    y_re:   7 =   6 - (  7* 22 -  15* 23)

    x_re:   0 =   0 + (  1* 16 -   9* 17)
    x_re:   4 =   4 + (  5* 20 -  13* 21)
    x_re:   2 =   2 + (  3* 18 -  11* 19)
    x_re:   6 =   6 + (  7* 22 -  15* 23)

    y_im:   9 =   8 - (  1* 17 +   9* 16)
    y_im:  13 =  12 - (  5* 21 +  13* 20)
    y_im:  11 =  10 - (  3* 19 +  11* 18)
    y_im:  15 =  14 - (  7* 23 +  15* 22)

    x_im:   8 =   8 + (  1* 17 +   9* 16)
    x_im:  12 =  12 + (  5* 21 +  13* 20)
    x_im:  10 =  10 + (  3* 19 +  11* 18)
    x_im:  14 =  14 + (  7* 23 +  15* 22)

    s_re_im: 16, 20 | 17, 21 | 18, 22 | 19, 23
     */
    vload4(s_re_im, &fpr_gm_tab[16]);

    vfmul(v1.val[0], t0.val[1], s_re_im.val[0]);
    vfmul(v1.val[1], t0.val[3], s_re_im.val[2]);

    vfmul(v2.val[0], t0.val[1], s_re_im.val[1]);
    vfmul(v2.val[1], t0.val[3], s_re_im.val[3]);

    vfms(v1.val[0], v1.val[0], t1.val[1], s_re_im.val[1]);
    vfms(v1.val[1], v1.val[1], t1.val[3], s_re_im.val[3]);

    vfma(v2.val[0], v2.val[0], t1.val[1], s_re_im.val[0]);
    vfma(v2.val[1], v2.val[1], t1.val[3], s_re_im.val[2]);

    vfsub(t0.val[1], t0.val[0], v1.val[0]);
    vfsub(t0.val[3], t0.val[2], v1.val[1]);
    vfsub(t1.val[1], t1.val[0], v2.val[0]);
    vfsub(t1.val[3], t1.val[2], v2.val[1]);

    vfadd(t0.val[0], t0.val[0], v1.val[0]);
    vfadd(t0.val[2], t0.val[2], v1.val[1]);
    vfadd(t1.val[0], t1.val[0], v2.val[0]);
    vfadd(t1.val[2], t1.val[2], v2.val[1]);

    vstore4(&f[0], t0);
    vstore4(&f[8], t1);
}

static void ZfN(FFT_log5)(fpr *f, unsigned logn)
{
    // Total SIMD register: 32 = 24 + 8
    float64x2x4_t s_re_im, tmp;                   // 8
    float64x2x4_t x_re, x_im, y_re, y_im, v1, v2; // 24
    float64x2x2_t s_tmp;
    const unsigned int falcon_n = 1 << logn;
    const unsigned int hn = falcon_n >> 1;
    unsigned i;

    for (unsigned j = 0; j < hn; j += 16)
    {
        // Level 3
        // x_re: 0->7
        // y_re: 8->15
        // x_im: 256->263
        // y_im: 264->271

        vloadx4(x_re, &f[j]);
        vloadx4(y_re, &f[j + 8]);
        vloadx4(x_im, &f[j + hn]);
        vloadx4(y_im, &f[j + hn + 8]);

        i = falcon_n + j;

        vload(s_re_im.val[0], &fpr_gm_tab[i >> 3]);

        vfmul_lane(v1.val[0], y_re.val[0], s_re_im.val[0], 0);
        vfmul_lane(v1.val[1], y_re.val[1], s_re_im.val[0], 0);
        vfmul_lane(v1.val[2], y_re.val[2], s_re_im.val[0], 0);
        vfmul_lane(v1.val[3], y_re.val[3], s_re_im.val[0], 0);

        vfms_lane(v1.val[0], v1.val[0], y_im.val[0], s_re_im.val[0], 1);
        vfms_lane(v1.val[1], v1.val[1], y_im.val[1], s_re_im.val[0], 1);
        vfms_lane(v1.val[2], v1.val[2], y_im.val[2], s_re_im.val[0], 1);
        vfms_lane(v1.val[3], v1.val[3], y_im.val[3], s_re_im.val[0], 1);

        vfmul_lane(v2.val[0], y_re.val[0], s_re_im.val[0], 1);
        vfmul_lane(v2.val[1], y_re.val[1], s_re_im.val[0], 1);
        vfmul_lane(v2.val[2], y_re.val[2], s_re_im.val[0], 1);
        vfmul_lane(v2.val[3], y_re.val[3], s_re_im.val[0], 1);

        vfma_lane(v2.val[0], v2.val[0], y_im.val[0], s_re_im.val[0], 0);
        vfma_lane(v2.val[1], v2.val[1], y_im.val[1], s_re_im.val[0], 0);
        vfma_lane(v2.val[2], v2.val[2], y_im.val[2], s_re_im.val[0], 0);
        vfma_lane(v2.val[3], v2.val[3], y_im.val[3], s_re_im.val[0], 0);

        vfsubx4(y_re, x_re, v1);
        vfsubx4(y_im, x_im, v2);

        vfaddx4(x_re, x_re, v1);
        vfaddx4(x_im, x_im, v2);

        // Level 2
        // x_re: 0->7
        // y_re: 8->15
        // x_im: 256->263
        // y_im: 264->271
        vloadx2(s_tmp, &fpr_gm_tab[i >> 2]);
        s_re_im.val[0] = s_tmp.val[0];
        s_re_im.val[1] = s_tmp.val[1];

        vfmul_lane(v1.val[0], x_re.val[2], s_re_im.val[0], 0);
        vfmul_lane(v1.val[1], x_re.val[3], s_re_im.val[0], 0);
        vfmul_lane(v1.val[2], y_re.val[2], s_re_im.val[1], 0);
        vfmul_lane(v1.val[3], y_re.val[3], s_re_im.val[1], 0);

        vfms_lane(v1.val[0], v1.val[0], x_im.val[2], s_re_im.val[0], 1);
        vfms_lane(v1.val[1], v1.val[1], x_im.val[3], s_re_im.val[0], 1);
        vfms_lane(v1.val[2], v1.val[2], y_im.val[2], s_re_im.val[1], 1);
        vfms_lane(v1.val[3], v1.val[3], y_im.val[3], s_re_im.val[1], 1);

        vfmul_lane(v2.val[0], x_re.val[2], s_re_im.val[0], 1);
        vfmul_lane(v2.val[1], x_re.val[3], s_re_im.val[0], 1);
        vfmul_lane(v2.val[2], y_re.val[2], s_re_im.val[1], 1);
        vfmul_lane(v2.val[3], y_re.val[3], s_re_im.val[1], 1);

        vfma_lane(v2.val[0], v2.val[0], x_im.val[2], s_re_im.val[0], 0);
        vfma_lane(v2.val[1], v2.val[1], x_im.val[3], s_re_im.val[0], 0);
        vfma_lane(v2.val[2], v2.val[2], y_im.val[2], s_re_im.val[1], 0);
        vfma_lane(v2.val[3], v2.val[3], y_im.val[3], s_re_im.val[1], 0);

        vfsub(y_re.val[2], y_re.val[0], v1.val[2]);
        vfsub(y_re.val[3], y_re.val[1], v1.val[3]);
        vfadd(x_re.val[2], y_re.val[0], v1.val[2]);
        vfadd(x_re.val[3], y_re.val[1], v1.val[3]);

        vfsub(y_re.val[0], x_re.val[0], v1.val[0]);
        vfsub(y_re.val[1], x_re.val[1], v1.val[1]);
        vfadd(x_re.val[0], x_re.val[0], v1.val[0]);
        vfadd(x_re.val[1], x_re.val[1], v1.val[1]);

        vfsub(y_im.val[2], y_im.val[0], v2.val[2]);
        vfsub(y_im.val[3], y_im.val[1], v2.val[3]);
        vfadd(x_im.val[2], y_im.val[0], v2.val[2]);
        vfadd(x_im.val[3], y_im.val[1], v2.val[3]);

        vfsub(y_im.val[0], x_im.val[0], v2.val[0]);
        vfsub(y_im.val[1], x_im.val[1], v2.val[1]);
        vfadd(x_im.val[0], x_im.val[0], v2.val[0]);
        vfadd(x_im.val[1], x_im.val[1], v2.val[1]);

        // Level 1
        // x_re: 0->3, 8->11
        // y_re: 4->7, 12->15
        // x_im: 256->259, 264->267
        // y_im: 260->263, 268->271
        vloadx4(s_re_im, &fpr_gm_tab[i >> 1]);

        vfmul_lane(v1.val[0], x_re.val[1], s_re_im.val[0], 0);
        vfmul_lane(v1.val[1], y_re.val[1], s_re_im.val[1], 0);
        vfmul_lane(v1.val[2], x_re.val[3], s_re_im.val[2], 0);
        vfmul_lane(v1.val[3], y_re.val[3], s_re_im.val[3], 0);

        vfms_lane(v1.val[0], v1.val[0], x_im.val[1], s_re_im.val[0], 1);
        vfms_lane(v1.val[1], v1.val[1], y_im.val[1], s_re_im.val[1], 1);
        vfms_lane(v1.val[2], v1.val[2], x_im.val[3], s_re_im.val[2], 1);
        vfms_lane(v1.val[3], v1.val[3], y_im.val[3], s_re_im.val[3], 1);

        vfmul_lane(v2.val[0], x_re.val[1], s_re_im.val[0], 1);
        vfmul_lane(v2.val[1], y_re.val[1], s_re_im.val[1], 1);
        vfmul_lane(v2.val[2], x_re.val[3], s_re_im.val[2], 1);
        vfmul_lane(v2.val[3], y_re.val[3], s_re_im.val[3], 1);

        vfma_lane(v2.val[0], v2.val[0], x_im.val[1], s_re_im.val[0], 0);
        vfma_lane(v2.val[1], v2.val[1], y_im.val[1], s_re_im.val[1], 0);
        vfma_lane(v2.val[2], v2.val[2], x_im.val[3], s_re_im.val[2], 0);
        vfma_lane(v2.val[3], v2.val[3], y_im.val[3], s_re_im.val[3], 0);

        vfsub(y_re.val[1], y_re.val[0], v1.val[1]);
        vfsub(y_re.val[3], y_re.val[2], v1.val[3]);
        vfadd(x_re.val[1], y_re.val[0], v1.val[1]);
        vfadd(x_re.val[3], y_re.val[2], v1.val[3]);

        vfsub(y_re.val[0], x_re.val[0], v1.val[0]);
        vfsub(y_re.val[2], x_re.val[2], v1.val[2]);
        vfadd(x_re.val[0], x_re.val[0], v1.val[0]);
        vfadd(x_re.val[2], x_re.val[2], v1.val[2]);

        vfsub(y_im.val[1], y_im.val[0], v2.val[1]);
        vfsub(y_im.val[3], y_im.val[2], v2.val[3]);
        vfadd(x_im.val[1], y_im.val[0], v2.val[1]);
        vfadd(x_im.val[3], y_im.val[2], v2.val[3]);

        vfsub(y_im.val[0], x_im.val[0], v2.val[0]);
        vfsub(y_im.val[2], x_im.val[2], v2.val[2]);
        vfadd(x_im.val[0], x_im.val[0], v2.val[0]);
        vfadd(x_im.val[2], x_im.val[2], v2.val[2]);

        // Level 0
        // Before Transpose
        // x_re: 0,1 | 4,5 | 8,9   | 12,13
        // y_re: 2,3 | 6,7 | 10,11 | 14,15
        // x_im: 256,257 | 260,261 | 264,265 | 268,269
        // y_im: 258,259 | 262,263 | 266,267 | 270,271
        transpose(x_re, x_re, tmp, 0, 1, 0);
        transpose(x_re, x_re, tmp, 2, 3, 1);
        transpose(y_re, y_re, tmp, 0, 1, 2);
        transpose(y_re, y_re, tmp, 2, 3, 3);

        transpose(x_im, x_im, tmp, 0, 1, 0);
        transpose(x_im, x_im, tmp, 2, 3, 1);
        transpose(y_im, y_im, tmp, 0, 1, 2);
        transpose(y_im, y_im, tmp, 2, 3, 3);

        // After Transpose
        // x_re: 0,4 | 1,5 | 8,12  | 9,13
        // y_re: 2,6 | 3,7 | 10,14 | 11,15
        // x_im: 256,260 | 257,261 | 264,268 | 265,269
        // y_im: 258,262 | 259,263 | 266,270 | 267,271
        vload4(s_re_im, &fpr_gm_tab[i]);

        vfmul(v1.val[0], x_re.val[1], s_re_im.val[0]);
        vfmul(v1.val[1], y_re.val[1], s_re_im.val[2]);
        vfmul(v2.val[0], x_re.val[1], s_re_im.val[1]);
        vfmul(v2.val[1], y_re.val[1], s_re_im.val[3]);

        vfms(v1.val[0], v1.val[0], x_im.val[1], s_re_im.val[1]);
        vfms(v1.val[1], v1.val[1], y_im.val[1], s_re_im.val[3]);
        vfma(v2.val[0], v2.val[0], x_im.val[1], s_re_im.val[0]);
        vfma(v2.val[1], v2.val[1], y_im.val[1], s_re_im.val[2]);

        vload4(s_re_im, &fpr_gm_tab[i + 8]);

        vfmul(v1.val[2], x_re.val[3], s_re_im.val[0]);
        vfmul(v1.val[3], y_re.val[3], s_re_im.val[2]);
        vfmul(v2.val[2], x_re.val[3], s_re_im.val[1]);
        vfmul(v2.val[3], y_re.val[3], s_re_im.val[3]);

        vfms(v1.val[2], v1.val[2], x_im.val[3], s_re_im.val[1]);
        vfms(v1.val[3], v1.val[3], y_im.val[3], s_re_im.val[3]);
        vfma(v2.val[2], v2.val[2], x_im.val[3], s_re_im.val[0]);
        vfma(v2.val[3], v2.val[3], y_im.val[3], s_re_im.val[2]);

        // x_re: 0,4 | 2,6 | 8,12 | 10,14
        // y_re: 1,5 | 3,7 | 9,13 | 11,15
        // x_im: 256,260 | 258,262 | 264,268 | 266,270
        // y_im: 257,261 | 259,263 | 265,269 | 267,271

        vfadd(tmp.val[0], x_re.val[0], v1.val[0]);
        vfsub(tmp.val[1], x_re.val[0], v1.val[0]);
        vfadd(tmp.val[2], y_re.val[0], v1.val[1]);
        vfsub(tmp.val[3], y_re.val[0], v1.val[1]);

        vstore4(&f[j], tmp);

        vfadd(tmp.val[0], x_re.val[2], v1.val[2]);
        vfsub(tmp.val[1], x_re.val[2], v1.val[2]);
        vfadd(tmp.val[2], y_re.val[2], v1.val[3]);
        vfsub(tmp.val[3], y_re.val[2], v1.val[3]);

        vstore4(&f[j + 8], tmp);

        vfadd(tmp.val[0], x_im.val[0], v2.val[0]);
        vfsub(tmp.val[1], x_im.val[0], v2.val[0]);
        vfadd(tmp.val[2], y_im.val[0], v2.val[1]);
        vfsub(tmp.val[3], y_im.val[0], v2.val[1]);

        vstore4(&f[j + hn], tmp);

        vfadd(tmp.val[0], x_im.val[2], v2.val[2]);
        vfsub(tmp.val[1], x_im.val[2], v2.val[2]);
        vfadd(tmp.val[2], y_im.val[2], v2.val[3]);
        vfsub(tmp.val[3], y_im.val[2], v2.val[3]);

        vstore4(&f[j + hn + 8], tmp);
    }
}

static void ZfN(FFT_logn1)(fpr *f, unsigned logn, const uint8_t level)
{
    // Total SIMD register: 26
    float64x2_t s_re_im;                          // 2
    float64x2x4_t x_re, x_im, y_re, y_im, v1, v2; // 24

    const unsigned int falcon_n = 1 << logn;
    const unsigned int hn = falcon_n >> 1;
    const unsigned l = level - 1;
    const unsigned distance = 1 << (l - 1);
    for (unsigned i = 0; i < hn; i += 1 << l)
    {
        vload(s_re_im, &fpr_gm_tab[(falcon_n + i) >> (l - 1)]);
        for (unsigned j = 0; j < distance; j += 8)
        {
            // Level 6
            // Distance: 16
            // hn = 32
            // falcon_n = 64
            // x_re: 0->3, 4->7
            // x_im: 32->35, 36-> 39
            // y_re: 16->19, 20->23
            // y_im: 48->51, 52 -> 55

            // Level 8
            // Distance: 64
            // hn = 128
            // x_re: 0->3, 4->7
            // x_im: 128->131, 132-> 135
            // y_re: 64->67, 68->71
            // y_im: 192->195, 196 -> 199
            vloadx4(x_re, &f[j]);
            vloadx4(y_re, &f[j + distance]);
            vloadx4(x_im, &f[j + hn]);
            vloadx4(y_im, &f[j + hn + distance]);

            vfmul_lane(v1.val[0], y_re.val[0], s_re_im, 0);
            vfmul_lane(v1.val[1], y_re.val[1], s_re_im, 0);
            vfmul_lane(v1.val[2], y_re.val[2], s_re_im, 0);
            vfmul_lane(v1.val[3], y_re.val[3], s_re_im, 0);

            vfms_lane(v1.val[0], v1.val[0], y_im.val[0], s_re_im, 1);
            vfms_lane(v1.val[1], v1.val[1], y_im.val[1], s_re_im, 1);
            vfms_lane(v1.val[2], v1.val[2], y_im.val[2], s_re_im, 1);
            vfms_lane(v1.val[3], v1.val[3], y_im.val[3], s_re_im, 1);

            vfmul_lane(v2.val[0], y_re.val[0], s_re_im, 1);
            vfmul_lane(v2.val[1], y_re.val[1], s_re_im, 1);
            vfmul_lane(v2.val[2], y_re.val[2], s_re_im, 1);
            vfmul_lane(v2.val[3], y_re.val[3], s_re_im, 1);

            vfma_lane(v2.val[0], v2.val[0], y_im.val[0], s_re_im, 0);
            vfma_lane(v2.val[1], v2.val[1], y_im.val[1], s_re_im, 0);
            vfma_lane(v2.val[2], v2.val[2], y_im.val[2], s_re_im, 0);
            vfma_lane(v2.val[3], v2.val[3], y_im.val[3], s_re_im, 0);

            vfsubx4(y_re, x_re, v1);
            vfsubx4(y_im, x_im, v2);

            vfaddx4(x_re, x_re, v1);
            vfaddx4(x_im, x_im, v2);

            vstorex4(&f[j], x_re);
            vstorex4(&f[j + distance], y_re);
            vstorex4(&f[j + hn], x_im);
            vstorex4(&f[j + hn + distance], y_im);
        }
    }
    // End function
}

static void ZfN(FFT_logn2)(fpr *f, unsigned logn, const uint8_t level)
{
    // Total: 26 = 16 + 8 + 2 register
    float64x2x4_t s_re_im;                        // 8
    float64x2x4_t x_re, x_im, y_re, y_im, v1, v2; // 16
    float64x2x2_t s_tmp;                          // 2
    // for storing instruction
    float64x2x2_t x_tmp, y_tmp;

    const unsigned int falcon_n = 1 << logn;
    const unsigned int hn = falcon_n >> 1;

    for (int l = level - 1; l > 4; l -= 2)
    {
        int distance = 1 << (l - 2);
        for (unsigned i = 0; i < hn; i += 1 << l)
        {
            vload(s_re_im.val[0], &fpr_gm_tab[(falcon_n + i) >> (l - 1)]);
            vloadx2(s_tmp, &fpr_gm_tab[(falcon_n + i) >> (l - 2)]);
            s_re_im.val[1] = s_tmp.val[0];
            s_re_im.val[2] = s_tmp.val[1];

            for (unsigned j = i; j < i + distance; j += 4)
            {
                // Level 7
                // x1_re: 0->3, 64->67
                // x1_im: 256->259, 320 -> 323
                // y_re: 128->131, 192->195
                // y_im: 384->387, 448 -> 451
                vloadx2(x_tmp, &f[j]);
                x_re.val[0] = x_tmp.val[0];
                x_re.val[1] = x_tmp.val[1];
                vloadx2(x_tmp, &f[j + distance]);
                x_re.val[2] = x_tmp.val[0];
                x_re.val[3] = x_tmp.val[1];

                vloadx2(y_tmp, &f[j + 2 * distance]);
                y_re.val[0] = y_tmp.val[0];
                y_re.val[1] = y_tmp.val[1];
                vloadx2(y_tmp, &f[j + 3 * distance]);
                y_re.val[2] = y_tmp.val[0];
                y_re.val[3] = y_tmp.val[1];

                vloadx2(x_tmp, &f[j + hn]);
                x_im.val[0] = x_tmp.val[0];
                x_im.val[1] = x_tmp.val[1];
                vloadx2(x_tmp, &f[j + hn + distance]);
                x_im.val[2] = x_tmp.val[0];
                x_im.val[3] = x_tmp.val[1];

                vloadx2(y_tmp, &f[j + hn + 2 * distance]);
                y_im.val[0] = y_tmp.val[0];
                y_im.val[1] = y_tmp.val[1];
                vloadx2(y_tmp, &f[j + hn + 3 * distance]);
                y_im.val[2] = y_tmp.val[0];
                y_im.val[3] = y_tmp.val[1];

                vfmul_lane(v1.val[0], y_re.val[0], s_re_im.val[0], 0);
                vfmul_lane(v1.val[1], y_re.val[1], s_re_im.val[0], 0);
                vfmul_lane(v1.val[2], y_re.val[2], s_re_im.val[0], 0);
                vfmul_lane(v1.val[3], y_re.val[3], s_re_im.val[0], 0);

                vfmul_lane(v2.val[0], y_re.val[0], s_re_im.val[0], 1);
                vfmul_lane(v2.val[1], y_re.val[1], s_re_im.val[0], 1);
                vfmul_lane(v2.val[2], y_re.val[2], s_re_im.val[0], 1);
                vfmul_lane(v2.val[3], y_re.val[3], s_re_im.val[0], 1);

                vfms_lane(v1.val[0], v1.val[0], y_im.val[0], s_re_im.val[0], 1);
                vfms_lane(v1.val[1], v1.val[1], y_im.val[1], s_re_im.val[0], 1);
                vfms_lane(v1.val[2], v1.val[2], y_im.val[2], s_re_im.val[0], 1);
                vfms_lane(v1.val[3], v1.val[3], y_im.val[3], s_re_im.val[0], 1);

                vfma_lane(v2.val[0], v2.val[0], y_im.val[0], s_re_im.val[0], 0);
                vfma_lane(v2.val[1], v2.val[1], y_im.val[1], s_re_im.val[0], 0);
                vfma_lane(v2.val[2], v2.val[2], y_im.val[2], s_re_im.val[0], 0);
                vfma_lane(v2.val[3], v2.val[3], y_im.val[3], s_re_im.val[0], 0);

                vfsubx4(y_re, x_re, v1);
                vfsubx4(y_im, x_im, v2);

                vfaddx4(x_re, x_re, v1);
                vfaddx4(x_im, x_im, v2);

                // Level 6
                // x_re: 0->3, 64->67
                // y_re: 128->131, 192 -> 195
                // x_im: 256->259, 320 -> 323
                // y_im: 384->387, 448 -> 451

                vfmul_lane(v1.val[0], x_re.val[2], s_re_im.val[1], 0);
                vfmul_lane(v1.val[1], x_re.val[3], s_re_im.val[1], 0);
                vfmul_lane(v1.val[2], y_re.val[2], s_re_im.val[2], 0);
                vfmul_lane(v1.val[3], y_re.val[3], s_re_im.val[2], 0);

                vfmul_lane(v2.val[0], x_re.val[2], s_re_im.val[1], 1);
                vfmul_lane(v2.val[1], x_re.val[3], s_re_im.val[1], 1);
                vfmul_lane(v2.val[2], y_re.val[2], s_re_im.val[2], 1);
                vfmul_lane(v2.val[3], y_re.val[3], s_re_im.val[2], 1);

                vfms_lane(v1.val[0], v1.val[0], x_im.val[2], s_re_im.val[1], 1);
                vfms_lane(v1.val[1], v1.val[1], x_im.val[3], s_re_im.val[1], 1);
                vfms_lane(v1.val[2], v1.val[2], y_im.val[2], s_re_im.val[2], 1);
                vfms_lane(v1.val[3], v1.val[3], y_im.val[3], s_re_im.val[2], 1);

                vfma_lane(v2.val[0], v2.val[0], x_im.val[2], s_re_im.val[1], 0);
                vfma_lane(v2.val[1], v2.val[1], x_im.val[3], s_re_im.val[1], 0);
                vfma_lane(v2.val[2], v2.val[2], y_im.val[2], s_re_im.val[2], 0);
                vfma_lane(v2.val[3], v2.val[3], y_im.val[3], s_re_im.val[2], 0);

                // Level 6
                // x1_re: 0->3, 128 -> 131
                // y1_re: 64->67, 192 -> 195
                // x1_im: 256 -> 259, 384->387
                // y1_im: 320->323, 448 -> 451
                vfadd(x_tmp.val[0], x_re.val[0], v1.val[0]);
                vfadd(x_tmp.val[1], x_re.val[1], v1.val[1]);
                vfsub(y_tmp.val[0], x_re.val[0], v1.val[0]);
                vfsub(y_tmp.val[1], x_re.val[1], v1.val[1]);

                vstorex2(&f[j], x_tmp);
                vstorex2(&f[j + distance], y_tmp);

                vfadd(x_tmp.val[0], y_re.val[0], v1.val[2]);
                vfadd(x_tmp.val[1], y_re.val[1], v1.val[3]);
                vfsub(y_tmp.val[0], y_re.val[0], v1.val[2]);
                vfsub(y_tmp.val[1], y_re.val[1], v1.val[3]);

                vstorex2(&f[j + 2 * distance], x_tmp);
                vstorex2(&f[j + 3 * distance], y_tmp);

                vfadd(x_tmp.val[0], x_im.val[0], v2.val[0]);
                vfadd(x_tmp.val[1], x_im.val[1], v2.val[1]);
                vfsub(y_tmp.val[0], x_im.val[0], v2.val[0]);
                vfsub(y_tmp.val[1], x_im.val[1], v2.val[1]);

                vstorex2(&f[j + hn], x_tmp);
                vstorex2(&f[j + hn + distance], y_tmp);

                vfadd(x_tmp.val[0], y_im.val[0], v2.val[2]);
                vfadd(x_tmp.val[1], y_im.val[1], v2.val[3]);
                vfsub(y_tmp.val[0], y_im.val[0], v2.val[2]);
                vfsub(y_tmp.val[1], y_im.val[1], v2.val[3]);

                vstorex2(&f[j + hn + 2 * distance], x_tmp);
                vstorex2(&f[j + hn + 3 * distance], y_tmp);
            }
        }
    }
}

static 
void ZfN(iFFT_log2)(fpr *f)
{
    /*
    y_re: 1 = (2 - 3) * 5 + (0 - 1) * 4
    y_im: 3 = (2 - 3) * 4 - (0 - 1) * 5
    x_re: 0 = 0 + 1
    x_im: 2 = 2 + 3

    Turn out this vectorize code is too short to be executed,
    the scalar version is consistently faster
     
    float64x2x2_t tmp;
    float64x2_t v, s, t;

    // 0: 0, 2
    // 1: 1, 3

    vload2(tmp, &f[0]);
    vload(s, &fpr_gm_tab[4]);

    vfsub(v, tmp.val[0], tmp.val[1]);
    vfadd(tmp.val[0], tmp.val[0], tmp.val[1]);

    // y_im: 3 = (2 - 3) * 4  - (0 - 1) * 5
    // y_re: 1 = (2 - 3) * 5  + (0 - 1) * 4
    vswap(t, v);

    vfmul_lane(tmp.val[1], s, t, 0);
    vfcmla_90(tmp.val[1], t, s);

    vfmuln(tmp.val[0], tmp.val[0], 0.5);
    vfmuln(tmp.val[1], tmp.val[1], 0.5);

    vswap(tmp.val[1], tmp.val[1]);

    vstore2(&f[0], tmp);
    */

    fpr x_re, x_im, y_re, y_im, s_re, s_im;
    x_re = f[0];
    y_re = f[1];
    x_im = f[2];
    y_im = f[3];
    s_re = fpr_gm_tab[4]/2;
    s_im = fpr_gm_tab[5]/2;

    f[0] = (x_re + y_re)/2;
    f[2] = (x_im + y_im)/2;

    x_re = x_re - y_re;
    x_im = x_im - y_im;

    f[1] = x_im * s_im + x_re * s_re;
    f[3] = x_im * s_re - x_re * s_im;
}

static void ZfN(iFFT_log3)(fpr *f)
{
    /*
    y_im: 5 = (4 - 5) *  8 - (0 - 1) *  9
    y_re: 1 = (4 - 5) *  9 + (0 - 1) *  8
    y_im: 7 = (6 - 7) * 10 - (2 - 3) * 11
    y_re: 3 = (6 - 7) * 11 + (2 - 3) * 10
    x_re: 0 = 0 + 1
    x_im: 4 = 4 + 5
    x_re: 2 = 2 + 3
    x_im: 6 = 6 + 7
     */

    float64x2x4_t t;
    float64x2x2_t v, s;

    // 0, 4
    // 1, 5
    // 2, 6
    // 3, 7
    vload4(t, &f[0]);
    vloadx2(s, &fpr_gm_tab[8]);

    vfsub(v.val[0], t.val[0], t.val[1]);
    vfsub(v.val[1], t.val[2], t.val[3]);

    vfadd(t.val[0], t.val[0], t.val[1]);
    vfadd(t.val[2], t.val[2], t.val[3]);

    // 5,1 <- 1,5
    // 7,3 <- 3,7
    vswap(v.val[0], v.val[0]);
    vswap(v.val[1], v.val[1]);

    vfmul_lane(t.val[1], s.val[0], v.val[0], 0);
    vfmul_lane(t.val[3], s.val[1], v.val[1], 0);

    vfcmla_90(t.val[1], v.val[0], s.val[0]);
    vfcmla_90(t.val[3], v.val[1], s.val[1]);

    /*
    y_im: 6 = (4 - 6) * 4 - (0 - 2) * 5
    y_re: 2 = (4 - 6) * 5 + (0 - 2) * 4
    y_im: 7 = (5 - 7) * 4 - (1 - 3) * 5
    y_re: 3 = (5 - 7) * 5 + (1 - 3) * 4

    x_re: 0 = 0 + 2
    x_im: 4 = 4 + 6
    x_re: 1 = 1 + 3
    x_im: 5 = 5 + 7
    */

    vfsub(v.val[0], t.val[0], t.val[2]);
    vfsub(v.val[1], t.val[1], t.val[3]);

    vfadd(t.val[0], t.val[0], t.val[2]);
    vfadd(t.val[1], t.val[1], t.val[3]);

    vload(s.val[0], &fpr_gm_tab[4]);

    vswap(v.val[0], v.val[0]);

    vfmuln(s.val[0], s.val[0], 0.25);
    vfmul_lane(t.val[2], s.val[0], v.val[0], 0);
    vfmul_lane(t.val[3], s.val[0], v.val[1], 0);

    vfcmla_90(t.val[2], v.val[0], s.val[0]);
    vfcmla_90(t.val[3], v.val[1], s.val[0]);

    vfmuln(t.val[0], t.val[0], 0.25);
    vfmuln(t.val[1], t.val[1], 0.25);
    

    vswap(t.val[1], t.val[1]);
    vswap(t.val[2], t.val[2]);
    vswap(t.val[3], t.val[3]);

    vstore4(&f[0], t);
}

static void ZfN(iFFT_log4)(fpr *f)
{
    // 0: 0, 4 | 8 , 12
    // 1: 1, 5 | 9 , 13
    // 2: 2, 6 | 10, 14
    // 3: 3, 7 | 11, 15
    /*
    y_re: 1 = (8 - 9) * 17 + (0 - 1) * 16
    y_re: 5 = (12 - 13) * 21 + (4 - 5) * 20
    y_re: 3 = (10 - 11) * 19 + (2 - 3) * 18
    y_re: 7 = (14 - 15) * 23 + (6 - 7) * 22

    y_im: 9 = (8 - 9) * 16 - (0 - 1) * 17
    y_im: 13 = (12 - 13) * 20 - (4 - 5) * 21
    y_im: 11 = (10 - 11) * 18 - (2 - 3) * 19
    y_im: 15 = (14 - 15) * 22 - (6 - 7) * 23

    x_re: 0 = 0 + 1
    x_re: 4 = 4 + 5
    x_re: 2 = 2 + 3
    x_re: 6 = 6 + 7

    x_im:  8 = 8 + 9
    x_im: 12 = 12 + 13
    x_im: 10 = 10 + 11
    x_im: 14 = 14 + 15
     */
    float64x2x4_t x_re_im, y_re_im, s_re_im, v, tmp; // 16
    float64x2x2_t s_tmp;
    float64x2_t s;

    vload4(x_re_im, &f[0]);
    vload4(y_re_im, &f[8]);

    // x: 0,4  | 1,5  |  2,6  |  3,7
    // y: 8,12 | 9,13 | 10,14 | 11,15

    vfsub(v.val[0], y_re_im.val[0], y_re_im.val[1]);
    vfsub(v.val[1], y_re_im.val[2], y_re_im.val[3]);
    vfsub(v.val[2], x_re_im.val[0], x_re_im.val[1]);
    vfsub(v.val[3], x_re_im.val[2], x_re_im.val[3]);

    vfadd(x_re_im.val[0], x_re_im.val[0], x_re_im.val[1]);
    vfadd(x_re_im.val[1], x_re_im.val[2], x_re_im.val[3]);
    vfadd(x_re_im.val[2], y_re_im.val[0], y_re_im.val[1]);
    vfadd(x_re_im.val[3], y_re_im.val[2], y_re_im.val[3]);

    vload4(s_re_im, &fpr_gm_tab[16]);

    vfmul(y_re_im.val[0], v.val[0], s_re_im.val[1]);
    vfmul(y_re_im.val[1], v.val[1], s_re_im.val[3]);
    vfmul(y_re_im.val[2], v.val[0], s_re_im.val[0]);
    vfmul(y_re_im.val[3], v.val[1], s_re_im.val[2]);

    vfma(y_re_im.val[0], y_re_im.val[0], v.val[2], s_re_im.val[0]);
    vfma(y_re_im.val[1], y_re_im.val[1], v.val[3], s_re_im.val[2]);
    vfms(y_re_im.val[2], y_re_im.val[2], v.val[2], s_re_im.val[1]);
    vfms(y_re_im.val[3], y_re_im.val[3], v.val[3], s_re_im.val[3]);

    /*
    y_re: 2 = (8 - 10) * 9 + (0 - 2) * 8
    y_re: 6 = (12 - 14) * 11 + (4 - 6) * 10
    y_re: 3 = (9 - 11) * 9 + (1 - 3) * 8
    y_re: 7 = (13 - 15) * 11 + (5 - 7) * 10

    y_im: 10 = (8 - 10) * 8 - (0 - 2) * 9
    y_im: 14 = (12 - 14) * 10 - (4 - 6) * 11
    y_im: 11 = (9 - 11) * 8 - (1 - 3) * 9
    y_im: 15 = (13 - 15) * 10 - (5 - 7) * 11

    x_re: 0 = 0 + 2
    x_re: 4 = 4 + 6
    x_re: 1 = 1 + 3
    x_re: 5 = 5 + 7

    x_im: 8 = 8 + 10
    x_im: 12 = 12 + 14
    x_im: 9 = 9 + 11
    x_im: 13 = 13 + 15
     */
    // x: 0,4 | 2,6 | 8,12 | 10,14
    // y: 1,5 | 3,7 | 9,13 | 11,15

    vload2(s_tmp, &fpr_gm_tab[8]);

    vfsub(v.val[0], x_re_im.val[2], x_re_im.val[3]);
    vfsub(v.val[1], y_re_im.val[2], y_re_im.val[3]);
    vfsub(v.val[2], x_re_im.val[0], x_re_im.val[1]);
    vfsub(v.val[3], y_re_im.val[0], y_re_im.val[1]);

    vfadd(x_re_im.val[0], x_re_im.val[0], x_re_im.val[1]);
    vfadd(x_re_im.val[1], y_re_im.val[0], y_re_im.val[1]);
    vfadd(x_re_im.val[2], x_re_im.val[2], x_re_im.val[3]);
    vfadd(x_re_im.val[3], y_re_im.val[2], y_re_im.val[3]);

    vfmul(y_re_im.val[0], v.val[0], s_tmp.val[1]);
    vfmul(y_re_im.val[1], v.val[1], s_tmp.val[1]);
    vfmul(y_re_im.val[2], v.val[0], s_tmp.val[0]);
    vfmul(y_re_im.val[3], v.val[1], s_tmp.val[0]);

    vfma(y_re_im.val[0], y_re_im.val[0], v.val[2], s_tmp.val[0]);
    vfma(y_re_im.val[1], y_re_im.val[1], v.val[3], s_tmp.val[0]);
    vfms(y_re_im.val[2], y_re_im.val[2], v.val[2], s_tmp.val[1]);
    vfms(y_re_im.val[3], y_re_im.val[3], v.val[3], s_tmp.val[1]);

    /*
    y_re: 4 = (8 - 12) * 5 + (0 - 4) * 4
    y_re: 5 = (9 - 13) * 5 + (1 - 5) * 4
    y_re: 6 = (10 - 14) * 5 + (2 - 6) * 4
    y_re: 7 = (11 - 15) * 5 + (3 - 7) * 4

    y_im: 12 = (8 - 12) * 4 - (0 - 4) * 5
    y_im: 13 = (9 - 13) * 4 - (1 - 5) * 5
    y_im: 14 = (10 - 14) * 4 - (2 - 6) * 5
    y_im: 15 = (11 - 15) * 4 - (3 - 7) * 5

    x_re: 0 = 0 + 4
    x_re: 1 = 1 + 5
    x_re: 2 = 2 + 6
    x_re: 3 = 3 + 7

    x_im: 8 = 8 + 12
    x_im: 9 = 9 + 13
    x_im: 10 = 10 + 14
    x_im: 11 = 11 + 15
     */
    // x: 0,4 | 1,5 |  8,12 |  9,13
    // y: 2,6 | 3,7 | 10,14 | 11,15

    transpose(x_re_im, x_re_im, tmp, 0, 1, 0);
    transpose(x_re_im, x_re_im, tmp, 2, 3, 1);
    transpose(y_re_im, y_re_im, tmp, 0, 1, 2);
    transpose(y_re_im, y_re_im, tmp, 2, 3, 3);

    // x: 0,1 | 4,5 |  8,9  | 12,13
    // y: 2,3 | 6,7 | 10,11 | 14,15

    vload(s, &fpr_gm_tab[4]);

    vfsub(v.val[0], x_re_im.val[2], x_re_im.val[3]);
    vfsub(v.val[1], y_re_im.val[2], y_re_im.val[3]);
    vfsub(v.val[2], x_re_im.val[0], x_re_im.val[1]);
    vfsub(v.val[3], y_re_im.val[0], y_re_im.val[1]);

    vfadd(x_re_im.val[0], x_re_im.val[0], x_re_im.val[1]);
    vfadd(x_re_im.val[1], y_re_im.val[0], y_re_im.val[1]);
    vfadd(x_re_im.val[2], x_re_im.val[2], x_re_im.val[3]);
    vfadd(x_re_im.val[3], y_re_im.val[2], y_re_im.val[3]);

    vfmuln(s, s, 0.12500000000);
    vfmuln(x_re_im.val[0], x_re_im.val[0], 0.12500000000);
    vfmuln(x_re_im.val[1], x_re_im.val[1], 0.12500000000);
    vfmuln(x_re_im.val[2], x_re_im.val[2], 0.12500000000);
    vfmuln(x_re_im.val[3], x_re_im.val[3], 0.12500000000);

    vfmul_lane(y_re_im.val[0], v.val[0], s, 1);
    vfmul_lane(y_re_im.val[1], v.val[1], s, 1);
    vfmul_lane(y_re_im.val[2], v.val[0], s, 0);
    vfmul_lane(y_re_im.val[3], v.val[1], s, 0);

    vfma_lane(y_re_im.val[0], y_re_im.val[0], v.val[2], s, 0);
    vfma_lane(y_re_im.val[1], y_re_im.val[1], v.val[3], s, 0);

    vfms_lane(y_re_im.val[2], y_re_im.val[2], v.val[2], s, 1);
    vfms_lane(y_re_im.val[3], y_re_im.val[3], v.val[3], s, 1);

    tmp.val[0] = x_re_im.val[0];
    tmp.val[1] = x_re_im.val[1];
    tmp.val[2] = y_re_im.val[0];
    tmp.val[3] = y_re_im.val[1];
    vstorex4(&f[0], tmp);
    tmp.val[0] = x_re_im.val[2];
    tmp.val[1] = x_re_im.val[3];
    tmp.val[2] = y_re_im.val[2];
    tmp.val[3] = y_re_im.val[3];
    vstorex4(&f[8], tmp);
}

static void ZfN(iFFT_log5)(fpr *f, const unsigned logn, const uint8_t last)
{
    // Total SIMD register: 28 = 16 + 8 + 4
    float64x2x4_t s_re_im, tmp1, tmp2;    // 8
    float64x2x4_t x_re, x_im, y_re, y_im; // 16
    float64x2x4_t v1, v2;                 // 8
    float64x2x2_t s_tmp;                  // 2

    const unsigned falcon_n = 1 << logn;
    const unsigned hn = falcon_n >> 1;
    unsigned i;

    // Level 0, 1, 2, 3
    for (unsigned j = 0; j < hn; j += 16)
    {
        i = falcon_n + j;

        vload4(x_re, &f[j]);
        vload4(y_re, &f[j + 8]);
        vload4(x_im, &f[j + hn]);
        vload4(y_im, &f[j + hn + 8]);
        vload4(s_re_im, &fpr_gm_tab[i]);

        // x_re: 0,4  | 1,5  |  2,6  |  3,7
        // y_re: 8,12 | 9,13 | 10,14 | 11,15
        // x_im: 64,68 | 65,69 | 66,70 | 67,71
        // y_im: 72,76 | 73,77 | 74,78 | 75,79
        vfsubx4_swap(v1, x_im, y_im, 0, 1, 2, 3);
        vfsubx4_swap(v2, x_re, y_re, 0, 1, 2, 3);

        vfaddx4_swap(x_re, x_re, y_re, 0, 1, 2, 3);
        vfaddx4_swap(x_im, x_im, y_im, 0, 1, 2, 3);

        vfmul(y_re.val[0], v1.val[0], s_re_im.val[1]);
        vfmul(y_re.val[1], v1.val[1], s_re_im.val[3]);
        vfmul(y_im.val[0], v1.val[0], s_re_im.val[0]);
        vfmul(y_im.val[1], v1.val[1], s_re_im.val[2]);

        vfma(y_re.val[0], y_re.val[0], v2.val[0], s_re_im.val[0]);
        vfma(y_re.val[1], y_re.val[1], v2.val[1], s_re_im.val[2]);
        vfms(y_im.val[0], y_im.val[0], v2.val[0], s_re_im.val[1]);
        vfms(y_im.val[1], y_im.val[1], v2.val[1], s_re_im.val[3]);

        vload4(s_re_im, &fpr_gm_tab[i + 8]);

        vfmul(y_re.val[2], v1.val[2], s_re_im.val[1]);
        vfmul(y_re.val[3], v1.val[3], s_re_im.val[3]);
        vfmul(y_im.val[2], v1.val[2], s_re_im.val[0]);
        vfmul(y_im.val[3], v1.val[3], s_re_im.val[2]);

        vfma(y_re.val[2], y_re.val[2], v2.val[2], s_re_im.val[0]);
        vfma(y_re.val[3], y_re.val[3], v2.val[3], s_re_im.val[2]);
        vfms(y_im.val[2], y_im.val[2], v2.val[2], s_re_im.val[1]);
        vfms(y_im.val[3], y_im.val[3], v2.val[3], s_re_im.val[3]);

        // x_re: 0,4 | 2,6 | 8,12 | 10,14
        // y_re: 1,5 | 3,7 | 9,13 | 11,15
        // x_im: 64,68 | 66,70 | 72,76 | 74,78
        // y_im: 65,69 | 67,71 | 73,77 | 75,79

        vfsubx4_swap(v1, x_im, y_im, 0, 1, 2, 3);
        vfsubx4_swap(v2, x_re, y_re, 0, 1, 2, 3);

        vfaddx4_swap(x_re, x_re, y_re, 0, 1, 2, 3);
        vfaddx4_swap(x_im, x_im, y_im, 0, 1, 2, 3);

        vload2(s_tmp, &fpr_gm_tab[i >> 1]);
        s_re_im.val[0] = s_tmp.val[0];
        s_re_im.val[1] = s_tmp.val[1];
        vload2(s_tmp, &fpr_gm_tab[(i + 8) >> 1]);
        s_re_im.val[2] = s_tmp.val[0];
        s_re_im.val[3] = s_tmp.val[1];

        vfmul(y_re.val[0], v1.val[0], s_re_im.val[1]);
        vfmul(y_re.val[1], v1.val[1], s_re_im.val[3]);
        vfmul(y_re.val[2], v1.val[2], s_re_im.val[1]);
        vfmul(y_re.val[3], v1.val[3], s_re_im.val[3]);

        vfmul(y_im.val[0], v1.val[0], s_re_im.val[0]);
        vfmul(y_im.val[1], v1.val[1], s_re_im.val[2]);
        vfmul(y_im.val[2], v1.val[2], s_re_im.val[0]);
        vfmul(y_im.val[3], v1.val[3], s_re_im.val[2]);

        vfma(y_re.val[0], y_re.val[0], v2.val[0], s_re_im.val[0]);
        vfma(y_re.val[1], y_re.val[1], v2.val[1], s_re_im.val[2]);
        vfma(y_re.val[2], y_re.val[2], v2.val[2], s_re_im.val[0]);
        vfma(y_re.val[3], y_re.val[3], v2.val[3], s_re_im.val[2]);

        vfms(y_im.val[0], y_im.val[0], v2.val[0], s_re_im.val[1]);
        vfms(y_im.val[1], y_im.val[1], v2.val[1], s_re_im.val[3]);
        vfms(y_im.val[2], y_im.val[2], v2.val[2], s_re_im.val[1]);
        vfms(y_im.val[3], y_im.val[3], v2.val[3], s_re_im.val[3]);

        // x_re: 0,4 | 8,12 | 1,5 | 9,13
        // y_re: 2,6 | 10,14 | 3,7 | 11,15
        // x_im: 64,68 | 72,76 | 65,69 | 73,77
        // y_im: 66,70 | 74,78 | 67,71 | 75,79

        transpose(x_re, x_re, tmp1, 0, 2, 0);
        transpose(x_re, x_re, tmp1, 1, 3, 1);
        transpose(y_re, y_re, tmp1, 0, 2, 2);
        transpose(y_re, y_re, tmp1, 1, 3, 3);

        transpose(x_im, x_im, tmp2, 0, 2, 0);
        transpose(x_im, x_im, tmp2, 1, 3, 1);
        transpose(y_im, y_im, tmp2, 0, 2, 2);
        transpose(y_im, y_im, tmp2, 1, 3, 3);

        // ZIP
        // x_re = 0,1 | 8,9 | 4,5 | 12,13
        // y_re = 2,3 | 10,11 | 6,7 | 14,15
        // x_im = 64,65 | 72,73 | 68,69 | 76,77
        // y_im = 66,67 | 74,75 | 70,71 | 78,79

        vfsubx4_swap(v1, x_im, y_im, 0, 2, 1, 3);
        vfsubx4_swap(v2, x_re, y_re, 0, 2, 1, 3);

        vfaddx4_swap(x_re, x_re, y_re, 0, 2, 1, 3);
        vfaddx4_swap(x_im, x_im, y_im, 0, 2, 1, 3);

        vloadx2(s_tmp, &fpr_gm_tab[i >> 2]);

        vfmul_lane(y_re.val[0], v1.val[0], s_tmp.val[0], 1);
        vfmul_lane(y_re.val[1], v1.val[1], s_tmp.val[1], 1);
        vfmul_lane(y_re.val[2], v1.val[2], s_tmp.val[0], 1);
        vfmul_lane(y_re.val[3], v1.val[3], s_tmp.val[1], 1);

        vfmul_lane(y_im.val[0], v1.val[0], s_tmp.val[0], 0);
        vfmul_lane(y_im.val[1], v1.val[1], s_tmp.val[1], 0);
        vfmul_lane(y_im.val[2], v1.val[2], s_tmp.val[0], 0);
        vfmul_lane(y_im.val[3], v1.val[3], s_tmp.val[1], 0);

        vfma_lane(y_re.val[0], y_re.val[0], v2.val[0], s_tmp.val[0], 0);
        vfma_lane(y_re.val[1], y_re.val[1], v2.val[1], s_tmp.val[1], 0);
        vfma_lane(y_re.val[2], y_re.val[2], v2.val[2], s_tmp.val[0], 0);
        vfma_lane(y_re.val[3], y_re.val[3], v2.val[3], s_tmp.val[1], 0);

        vfms_lane(y_im.val[0], y_im.val[0], v2.val[0], s_tmp.val[0], 1);
        vfms_lane(y_im.val[1], y_im.val[1], v2.val[1], s_tmp.val[1], 1);
        vfms_lane(y_im.val[2], y_im.val[2], v2.val[2], s_tmp.val[0], 1);
        vfms_lane(y_im.val[3], y_im.val[3], v2.val[3], s_tmp.val[1], 1);

        // x_re: 0,1 |  8,9  | 2,3 | 10,11
        // y_re: 4,5 | 12,13 | 6,7 | 14,15
        // x_im: 64,65 | 72,73 | 66,67 | 74,75
        // y_im: 68,69 | 76,77 | 70,71 | 78,79
        vfsubx4_swap(v1, x_im, y_im, 0, 1, 2, 3);
        vfsubx4_swap(v2, x_re, y_re, 0, 1, 2, 3);

        vfaddx4_swap(x_re, x_re, y_re, 0, 1, 2, 3);
        vfaddx4_swap(x_im, x_im, y_im, 0, 1, 2, 3);

        vload(s_re_im.val[0], &fpr_gm_tab[i >> 3]);

        if (last)
        {
            vfmulnx4(x_re, x_re, fpr_p2_tab[logn]);
            vfmulnx4(x_im, x_im, fpr_p2_tab[logn]);
            vfmuln(s_re_im.val[0], s_re_im.val[0], fpr_p2_tab[logn]);
        }

        vfmulx4_lane(y_re, v1, s_re_im.val[0], 1);
        vfmulx4_lane(y_im, v1, s_re_im.val[0], 0);

        vfmax4_lane(y_re, y_re, v2, s_re_im.val[0], 0);
        vfmsx4_lane(y_im, y_im, v2, s_re_im.val[0], 1);

        // x_re: 0,1 | 2,3 | 4,5 | 6,7
        // y_re: 8,9 | 10,11 | 12,13 | 14,15
        // x_im: 64,65 | 66,67 | 68,69 | 70,71
        // y_im: 72,73 | 74,75 | 76,77 | 78,79

        vstorex4(&f[j], x_re);
        vstorex4(&f[j + 8], y_re);
        vstorex4(&f[j + hn], x_im);
        vstorex4(&f[j + hn + 8], y_im);
    }
}

static void ZfN(iFFT_logn1)(fpr *f, const unsigned logn, const uint8_t last)
{
    /* Level 6
    y_re: 16 = (32 - 48) * 5 + (0 - 16) * 4
    y_re: 17 = (33 - 49) * 5 + (1 - 17) * 4
    y_re: 18 = (34 - 50) * 5 + (2 - 18) * 4
    y_re: 19 = (35 - 51) * 5 + (3 - 19) * 4
    y_re: 20 = (36 - 52) * 5 + (4 - 20) * 4
    y_re: 21 = (37 - 53) * 5 + (5 - 21) * 4
    y_re: 22 = (38 - 54) * 5 + (6 - 22) * 4
    y_re: 23 = (39 - 55) * 5 + (7 - 23) * 4

    y_im: 48 = (32 - 48) * 4 - (0 - 16) * 5
    y_im: 49 = (33 - 49) * 4 - (1 - 17) * 5
    y_im: 50 = (34 - 50) * 4 - (2 - 18) * 5
    y_im: 51 = (35 - 51) * 4 - (3 - 19) * 5
    y_im: 52 = (36 - 52) * 4 - (4 - 20) * 5
    y_im: 53 = (37 - 53) * 4 - (5 - 21) * 5
    y_im: 54 = (38 - 54) * 4 - (6 - 22) * 5
    y_im: 55 = (39 - 55) * 4 - (7 - 23) * 5

    x_re: 0 = 0 + 16
    x_re: 1 = 1 + 17
    x_re: 2 = 2 + 18
    x_re: 3 = 3 + 19
    x_re: 4 = 4 + 20
    x_re: 5 = 5 + 21
    x_re: 6 = 6 + 22
    x_re: 7 = 7 + 23

    x_im: 32 = 32 + 48
    x_im: 33 = 33 + 49
    x_im: 34 = 34 + 50
    x_im: 35 = 35 + 51
    x_im: 36 = 36 + 52
    x_im: 37 = 37 + 53
    x_im: 38 = 38 + 54
    x_im: 39 = 39 + 55
     */

    float64x2x4_t x_re, x_im, y_re, y_im, v1, v2;
    float64x2_t s_re_im;
    const unsigned falcon_n = 1 << logn;
    const unsigned hn = falcon_n >> 1;
    const unsigned distance = hn >> 1;
    for (unsigned j = 0; j < distance; j += 8)
    {
        // x_re: 0 -> 7
        // y_re: 16->23
        // x_im: 32->39
        // y_im: 48->55
        vloadx4(x_re, &f[j]);
        vloadx4(y_re, &f[j + distance]);
        vloadx4(x_im, &f[j + hn]);
        vloadx4(y_im, &f[j + hn + distance]);
        vload(s_re_im, &fpr_gm_tab[4]);

        vfsubx4(v1, x_im, y_im);
        vfsubx4(v2, x_re, y_re);

        vfaddx4(x_re, x_re, y_re);
        vfaddx4(x_im, x_im, y_im);

        if (last)
        {
            vfmulnx4(x_re, x_re, fpr_p2_tab[logn]);
            vfmulnx4(x_im, x_im, fpr_p2_tab[logn]);
            vfmuln(s_re_im, s_re_im, fpr_p2_tab[logn]);
        }

        vfmulx4_lane(y_re, v1, s_re_im, 1);
        vfmulx4_lane(y_im, v1, s_re_im, 0);
        vfmax4_lane(y_re, y_re, v2, s_re_im, 0);
        vfmsx4_lane(y_im, y_im, v2, s_re_im, 1);

        vstorex4(&f[j], x_re);
        vstorex4(&f[j + distance], y_re);
        vstorex4(&f[j + hn], x_im);
        vstorex4(&f[j + hn + distance], y_im);
    }
}

static void ZfN(iFFT_logn2)(fpr *f, const unsigned logn, const uint8_t level, uint8_t last)
{
    // Total SIMD registers: 27 = 24 + 3
    float64x2x4_t x_re, y_re, x_im, y_im, v1, v2; // 24
    float64x2x3_t s_re_im;                        // 3
    float64x2x2_t x_tmp, y_tmp, s_tmp;            // 6
    unsigned distance;
    const unsigned falcon_n = 1 << logn;
    const unsigned hn = falcon_n >> 1;

    for (uint8_t l = 4; l < (logn - 1 - level); l += 2)
    {
        distance = 1 << l;
        last -= 1;
        for (unsigned i = 0; i < hn; i += 1 << (l + 2))
        {
            vloadx2(s_tmp, &fpr_gm_tab[(falcon_n + i) >> l]);
            s_re_im.val[0] = s_tmp.val[0];
            s_re_im.val[1] = s_tmp.val[1];
            vload(s_re_im.val[2], &fpr_gm_tab[(falcon_n + i) >> (l + 1)]);
            if (!last)
            {
                vfmuln(s_re_im.val[2], s_re_im.val[2], fpr_p2_tab[logn]);
            }
            for (unsigned j = i; j < i + distance; j += 4)
            {
                vloadx2(x_tmp, &f[j]);
                x_re.val[0] = x_tmp.val[0];
                x_re.val[1] = x_tmp.val[1];
                vloadx2(y_tmp, &f[j + distance]);
                y_re.val[0] = y_tmp.val[0];
                y_re.val[1] = y_tmp.val[1];

                vloadx2(x_tmp, &f[j + 2 * distance]);
                x_re.val[2] = x_tmp.val[0];
                x_re.val[3] = x_tmp.val[1];
                vloadx2(y_tmp, &f[j + 3 * distance]);
                y_re.val[2] = y_tmp.val[0];
                y_re.val[3] = y_tmp.val[1];

                vloadx2(x_tmp, &f[j + hn]);
                x_im.val[0] = x_tmp.val[0];
                x_im.val[1] = x_tmp.val[1];
                vloadx2(y_tmp, &f[j + hn + distance]);
                y_im.val[0] = y_tmp.val[0];
                y_im.val[1] = y_tmp.val[1];

                vloadx2(x_tmp, &f[j + hn + 2 * distance]);
                x_im.val[2] = x_tmp.val[0];
                x_im.val[3] = x_tmp.val[1];
                vloadx2(y_tmp, &f[j + hn + 3 * distance]);
                y_im.val[2] = y_tmp.val[0];
                y_im.val[3] = y_tmp.val[1];

                // x_re: 0,1 | 2,3 | 32,33 | 34,35
                // y_re: 16,17 | 18,19 | 48,49 | 50,51
                // x_im: 256 -> 259 | 288 -> 291
                // y_im: 272 -> 275 | 304 -> 307

                vfsubx4(v1, x_im, y_im);
                vfsubx4(v2, x_re, y_re);

                vfaddx4(x_re, x_re, y_re);
                vfaddx4(x_im, x_im, y_im);

                vfmul_lane(y_re.val[0], v1.val[0], s_re_im.val[0], 1);
                vfmul_lane(y_re.val[1], v1.val[1], s_re_im.val[0], 1);
                vfmul_lane(y_re.val[2], v1.val[2], s_re_im.val[1], 1);
                vfmul_lane(y_re.val[3], v1.val[3], s_re_im.val[1], 1);

                vfmul_lane(y_im.val[0], v1.val[0], s_re_im.val[0], 0);
                vfmul_lane(y_im.val[1], v1.val[1], s_re_im.val[0], 0);
                vfmul_lane(y_im.val[2], v1.val[2], s_re_im.val[1], 0);
                vfmul_lane(y_im.val[3], v1.val[3], s_re_im.val[1], 0);

                vfma_lane(y_re.val[0], y_re.val[0], v2.val[0], s_re_im.val[0], 0);
                vfma_lane(y_re.val[1], y_re.val[1], v2.val[1], s_re_im.val[0], 0);
                vfma_lane(y_re.val[2], y_re.val[2], v2.val[2], s_re_im.val[1], 0);
                vfma_lane(y_re.val[3], y_re.val[3], v2.val[3], s_re_im.val[1], 0);

                vfms_lane(y_im.val[0], y_im.val[0], v2.val[0], s_re_im.val[0], 1);
                vfms_lane(y_im.val[1], y_im.val[1], v2.val[1], s_re_im.val[0], 1);
                vfms_lane(y_im.val[2], y_im.val[2], v2.val[2], s_re_im.val[1], 1);
                vfms_lane(y_im.val[3], y_im.val[3], v2.val[3], s_re_im.val[1], 1);

                // print_vector(x_re, x_im, y_re, y_im);

                // x_re: 0 -> 3 | 32 -> 35
                // y_re: 16 -> 19 | 48 -> 51
                // x_im: 256 -> 259 | 288 -> 291
                // y_im: 272 -> 275 | 304 -> 307

                vfsubx4_swap(v1, x_im, y_im, 0, 2, 1, 3);
                vfsubx4_swap(v2, x_re, y_re, 0, 2, 1, 3);

                vfaddx4_swap(x_re, x_re, y_re, 0, 2, 1, 3);
                vfaddx4_swap(x_im, x_im, y_im, 0, 2, 1, 3);

                if (!last)
                {
                    // printf("div %u\n", last);
                    vfmulnx4(x_re, x_re, fpr_p2_tab[logn]);
                    vfmulnx4(x_im, x_im, fpr_p2_tab[logn]);
                }

                vfmulx4_lane(y_re, v1, s_re_im.val[2], 1);
                vfmulx4_lane(y_im, v1, s_re_im.val[2], 0);

                vfmax4_lane(y_re, y_re, v2, s_re_im.val[2], 0);
                vfmsx4_lane(y_im, y_im, v2, s_re_im.val[2], 1);

                // x_re: 0->3 | 16 -> 19
                // y_re: 32->35 | 48 -> 51
                // x_im: 256->259 | 272->275
                // y_im: 288->291 | 304->307

                x_tmp.val[0] = x_re.val[0];
                x_tmp.val[1] = x_re.val[1];
                vstorex2(&f[j], x_tmp);
                x_tmp.val[0] = x_re.val[2];
                x_tmp.val[1] = x_re.val[3];
                vstorex2(&f[j + distance], x_tmp);
                y_tmp.val[0] = y_re.val[0];
                y_tmp.val[1] = y_re.val[1];
                vstorex2(&f[j + 2 * distance], y_tmp);
                y_tmp.val[0] = y_re.val[2];
                y_tmp.val[1] = y_re.val[3];
                vstorex2(&f[j + 3 * distance], y_tmp);

                x_tmp.val[0] = x_im.val[0];
                x_tmp.val[1] = x_im.val[1];
                vstorex2(&f[j + hn], x_tmp);
                x_tmp.val[0] = x_im.val[2];
                x_tmp.val[1] = x_im.val[3];
                vstorex2(&f[j + hn + distance], x_tmp);
                y_tmp.val[0] = y_im.val[0];
                y_tmp.val[1] = y_im.val[1];
                vstorex2(&f[j + hn + 2 * distance], y_tmp);
                y_tmp.val[0] = y_im.val[2];
                y_tmp.val[1] = y_im.val[3];
                vstorex2(&f[j + hn + 3 * distance], y_tmp);
            }
        }
    }
    // End function
}

/*
 * Support logn from [1, 10]
 * Can be easily extended to logn > 10
 */
void ZfN(FFT)(fpr *f, const unsigned logn)
{
    uint8_t level = logn;
    switch (logn)
    {
    case 2:
        ZfN(FFT_log2)(f);
        break;

    case 3:
        ZfN(FFT_log3)(f);
        break;

    case 4:
        ZfN(FFT_log4)(f);
        break;

    case 5:
        ZfN(FFT_log5)(f, logn);
        break;

    case 6:
        ZfN(FFT_logn1)(f, logn, level--);
        ZfN(FFT_log5)(f, logn);
        break;

    case 7:
    case 9:
        ZfN(FFT_logn2)(f, logn, level);
        ZfN(FFT_log5)(f, logn);
        break;

    case 8:
    case 10:
        ZfN(FFT_logn1)(f, logn, level--);
        ZfN(FFT_logn2)(f, logn, level);
        ZfN(FFT_log5)(f, logn);
        break;
    
    default:
        break;
    }
}

void ZfN(iFFT)(fpr *f, const unsigned logn)
{
    const uint8_t level = (logn - 5) & 1;

    switch (logn)
    {
    case 2:
        ZfN(iFFT_log2)(f);
        break;

    case 3:
        ZfN(iFFT_log3)(f);
        break;

    case 4:
        ZfN(iFFT_log4)(f);
        break;

    case 5:
        ZfN(iFFT_log5)(f, logn, 1);
        break;

    case 6:
        ZfN(iFFT_log5)(f, logn, 0);
        ZfN(iFFT_logn1)(f, logn, 1);
        break;

    case 7:
    case 9:
        ZfN(iFFT_log5)(f, logn, 0);
        ZfN(iFFT_logn2)(f, logn, level, 1);
        break;

    case 8:
    case 10:
        ZfN(iFFT_log5)(f, logn, 0);
        ZfN(iFFT_logn2)(f, logn, level, 0);
        ZfN(iFFT_logn1)(f, logn, 1);
        break;
    
    default:
        break;
    }
}
