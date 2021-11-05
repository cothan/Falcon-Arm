#include "inner.h"
#include "vfpr.h"
#include <stdio.h>

// TODO: extend iFFT and FFT to support FALCON1024
void PQCLEAN_FALCON512_NEON_iFFT(fpr *f)
{
    // Total: 32 = 16 + 8 + 8 register
    float64x2x4_t s_re_im, tmp;           // 8
    float64x2x4_t x_re, x_im, y_re, y_im; // 16
    float64x2x4_t v_re, v_im;             // 8
    float64x2x2_t s_tmp;                  // 2
    // Level 4, 5
    float64x2x2_t x_tmp, y_tmp;
    // Level 6, 7
    float64x2_t div_n;

    const unsigned int hn = FALCON_N >> 1;

    // Level 0, 1, 2, 3
    for (int j = 0; j < FALCON_N / 2; j += 16)
    {
        // Level 0
        // x_re = 0, 4 | 2, 6 | 8, 12 | 10, 14
        // y_re = 1, 5 | 3, 7 | 9, 13 | 11, 15
        // x_im = 256, 260 | 258, 262 | 264, 268 | 266, 270
        // y_im = 257, 261 | 259, 263 | 265, 269 | 267, 271

        // This assignment is free
        vload4(tmp, &f[j]);
        x_re.val[0] = tmp.val[0];
        x_re.val[1] = tmp.val[2];
        y_re.val[0] = tmp.val[1];
        y_re.val[1] = tmp.val[3];

        vload4(tmp, &f[j + 8]);
        x_re.val[2] = tmp.val[0];
        x_re.val[3] = tmp.val[2];
        y_re.val[2] = tmp.val[1];
        y_re.val[3] = tmp.val[3];

        vload4(tmp, &f[j + hn]);
        x_im.val[0] = tmp.val[0];
        x_im.val[1] = tmp.val[2];
        y_im.val[0] = tmp.val[1];
        y_im.val[1] = tmp.val[3];

        vload4(tmp, &f[j + hn + 8]);
        x_im.val[2] = tmp.val[0];
        x_im.val[3] = tmp.val[2];
        y_im.val[2] = tmp.val[1];
        y_im.val[3] = tmp.val[3];

        ////////
        // x - y
        ////////

        //  1,5 <= 0,4  -  1,5 |   3,7 <=   2,6 -  3,7
        // 9,13 <= 8,12 - 9,13 | 11,15 <= 10,14 - 11,15
        vfsubx4(v_re, x_re, y_re);

        // 257,261 <= 256,260 - 257,261 | 259,263 <= 258,262 - 259,263
        // 265,269 <= 264,268 - 265,269 | 267,271 <= 266,270 - 267,271
        vfsubx4(v_im, x_im, y_im);

        ////////
        // x + y
        ////////

        // x_re: 0,4 | 2,6 | 8,12 | 10,14
        // 0,4  <= 0,4  +  1,5 |   2,6 <=   2,6 +   3,7
        // 8,12 <= 8,12 + 9,13 | 10,14 <= 10,14 + 11,15
        vfaddx4(x_re, x_re, y_re);

        // x_im: 256,260 | 258,262 | 264,268 | 266,270
        // 256,260 <= 256,260 + 257,261 | 258,262 <= 258,262 + 259,263
        // 264,268 <= 264,268 + 265,269 | 266,270 <= 266,270 + 267,271
        vfaddx4(x_im, x_im, y_im);

        // s * (x - y) = s*v = (s_re + i*s_im)(v_re + i*v_im)
        // y_re:  v_re*s_re + v_im*s_im
        // y_im: -v_re*s_im + v_im*s_re

        // s_re_im = 512 -> 519
        // s_re: 512 -> 518, step 2
        // s_im: 513 -> 519, step 2
        vload4(s_re_im, &fpr_gm_tab[FALCON_N + j]);

        // y_re: 1,5 | 3,7 | 9,13 | 11,15
        // y_re = v_im*s_im + v_re*s_re
        vfmul(tmp.val[0], v_im.val[0], s_re_im.val[1]);
        vfmul(tmp.val[1], v_im.val[1], s_re_im.val[3]);
        vfma(y_re.val[0], tmp.val[0], v_re.val[0], s_re_im.val[0]);
        vfma(y_re.val[1], tmp.val[1], v_re.val[1], s_re_im.val[2]);

        vfmul(tmp.val[0], v_im.val[0], s_re_im.val[0]);
        vfmul(tmp.val[1], v_im.val[1], s_re_im.val[2]);
        vfms(y_im.val[0], tmp.val[0], v_re.val[0], s_re_im.val[1]);
        vfms(y_im.val[1], tmp.val[1], v_re.val[1], s_re_im.val[3]);

        // y_im: 257,261 | 259,263 | 265,269 | 267,271
        // y_im = v_im*s_re - v_re*s_im
        vload4(s_re_im, &fpr_gm_tab[FALCON_N + j + 8]);

        vfmul(tmp.val[2], v_im.val[2], s_re_im.val[1]);
        vfmul(tmp.val[3], v_im.val[3], s_re_im.val[3]);
        vfma(y_re.val[2], tmp.val[2], v_re.val[2], s_re_im.val[0]);
        vfma(y_re.val[3], tmp.val[3], v_re.val[3], s_re_im.val[2]);

        vfmul(tmp.val[2], v_im.val[2], s_re_im.val[0]);
        vfmul(tmp.val[3], v_im.val[3], s_re_im.val[2]);
        vfms(y_im.val[2], tmp.val[2], v_re.val[2], s_re_im.val[1]);
        vfms(y_im.val[3], tmp.val[3], v_re.val[3], s_re_im.val[3]);

        // Level 1
        // x_re = 0,4 | 2,6 | 8,12 | 10,14
        // y_re = 1,5 | 3,7 | 9,13 | 11,15
        // x_im = 256,260 | 258,262 | 264,268 | 266,270
        // y_im = 257,261 | 259,263 | 265,269 | 267,271

        vload2(s_tmp, &fpr_gm_tab[(FALCON_N + j) >> 1]);
        s_re_im.val[0] = s_tmp.val[0];
        s_re_im.val[1] = s_tmp.val[1];
        vload2(s_tmp, &fpr_gm_tab[(FALCON_N + j + 8) >> 1]);
        s_re_im.val[2] = s_tmp.val[0];
        s_re_im.val[3] = s_tmp.val[1];

        ////////
        // x - y
        ////////

        // 0,4 - 2,6 | 8,12 - 10,14
        // 1,5 - 3,7 | 9,13 - 11,15
        vfsubx4_swap(v_re, x_re, y_re, 0, 1, 2, 3);
        // 256,260 - 258,262 | 264,268 - 266,270
        // 257,261 - 259,263 | 265,269 - 267,271
        vfsubx4_swap(v_im, x_im, y_im, 0, 1, 2, 3);

        ////////
        // x + y
        ////////

        // x_re: 0,4 | 8,12 | 1,5 | 9,13
        // 0,4 <= 0,4 + 2,6 | 8,12 <= 8,12 + 10,14
        // 1,5 <= 1,5 + 3,7 | 9,13 <= 9,13 + 11,15
        vfaddx4_swap(x_re, x_re, y_re, 0, 1, 2, 3);

        // x_im: 256, 260 | 264, 268 | 257, 261 | 265, 269
        // 256,260 <= 256,260 + 258,262 | 264,268 <= 264,268 + 266,270
        // 257,261 <= 257,261 + 259,263 | 265,269 <= 265,269 + 267,271
        vfaddx4_swap(x_im, x_im, y_im, 0, 1, 2, 3);

        // Calculate y
        // v_im*s_im
        vfmul(tmp.val[0], v_im.val[0], s_re_im.val[1]);
        vfmul(tmp.val[1], v_im.val[1], s_re_im.val[3]);
        vfmul(tmp.val[2], v_im.val[2], s_re_im.val[1]);
        vfmul(tmp.val[3], v_im.val[3], s_re_im.val[3]);

        // v_im*s_im + v_re*s_re
        // y_re: 2,6 | 10,14 | 3,7 | 11,15
        vfma(y_re.val[0], tmp.val[0], v_re.val[0], s_re_im.val[0]);
        vfma(y_re.val[1], tmp.val[1], v_re.val[1], s_re_im.val[2]);
        vfma(y_re.val[2], tmp.val[2], v_re.val[2], s_re_im.val[0]);
        vfma(y_re.val[3], tmp.val[3], v_re.val[3], s_re_im.val[2]);

        // v_im*s_re
        vfmul(tmp.val[0], v_im.val[0], s_re_im.val[0]);
        vfmul(tmp.val[1], v_im.val[1], s_re_im.val[2]);
        vfmul(tmp.val[2], v_im.val[2], s_re_im.val[0]);
        vfmul(tmp.val[3], v_im.val[3], s_re_im.val[2]);

        // v_im*s_re - v_re*s_im
        // y_im: 258,262 | 266,270 | 259,263 | 267,271
        vfms(y_im.val[0], tmp.val[0], v_re.val[0], s_re_im.val[1]);
        vfms(y_im.val[1], tmp.val[1], v_re.val[1], s_re_im.val[3]);
        vfms(y_im.val[2], tmp.val[2], v_re.val[2], s_re_im.val[1]);
        vfms(y_im.val[3], tmp.val[3], v_re.val[3], s_re_im.val[3]);

        // Level 2
        // Before Transpose
        // x_re = 0,4 | 8,12 | 1,5 | 9,13
        // y_re = 2,6 | 10,14 | 3,7 | 11,15
        // x_im = 256,260 | 264,268 | 257,261 | 265,269
        // y_im = 258,262 | 266,270 | 259,263 | 267,271

        transpose(x_re, x_re, tmp, 0, 2, 0);
        transpose(x_re, x_re, tmp, 1, 3, 1);
        transpose(y_re, y_re, tmp, 0, 2, 2);
        transpose(y_re, y_re, tmp, 1, 3, 3);

        transpose(x_im, x_im, tmp, 0, 2, 0);
        transpose(x_im, x_im, tmp, 1, 3, 1);
        transpose(y_im, y_im, tmp, 0, 2, 2);
        transpose(y_im, y_im, tmp, 1, 3, 3);

        // After Transpose
        // x_re = 0,1 | 8,9 | 4,5 | 12,13
        // y_re = 2,3 | 10,11 | 6,7 | 14,15
        // x_im = 256,257 | 264,265 | 260,261 | 268,269
        // y_im = 258,259 | 266,267 | 262,263 | 270,271

        vload2(s_tmp, &fpr_gm_tab[(FALCON_N + j) >> 2]);
        s_re_im.val[0] = s_tmp.val[0];
        s_re_im.val[1] = s_tmp.val[1];

        ////////
        // x - y
        ////////

        // 0,1 - 4,5   | 2,3 - 6,7
        // 8,9 - 12,13 | 10,11 - 14,15
        vfsub(v_re.val[0], x_re.val[0], x_re.val[2]);
        vfsub(v_re.val[2], x_re.val[1], x_re.val[3]);
        vfsub(v_re.val[1], y_re.val[0], y_re.val[2]);
        vfsub(v_re.val[3], y_re.val[1], y_re.val[3]);

        // 256,257 - 260,261 | 258,259 - 262,263
        // 264,265 - 268,269 | 266,267 - 270,271
        vfsub(v_im.val[0], x_im.val[0], x_im.val[2]);
        vfsub(v_im.val[2], x_im.val[1], x_im.val[3]);
        vfsub(v_im.val[1], y_im.val[0], y_im.val[2]);
        vfsub(v_im.val[3], y_im.val[1], y_im.val[3]);

        ////////
        // x + y
        ////////

        // x_re: 0,1 | 2,3 |  8,9 | 10,11
        // 0,1 <= 0,1 + 4,5  | 2,3 <= 2,3 + 6,7
        // 8,9 <= 8,9 + 12,13| 10,11 <= 10,11 + 14,15
        vfadd(x_re.val[0], x_re.val[0], x_re.val[2]);
        vfadd(x_re.val[2], x_re.val[1], x_re.val[3]);
        vfadd(x_re.val[1], y_re.val[0], y_re.val[2]);
        vfadd(x_re.val[3], y_re.val[1], y_re.val[3]);

        // x_im: 256,257 | 258,259 | 264,265 | 266,267
        // 256,257 <= 256,257 + 260,261 | 258,259 <= 258,259 + 262,263
        // 264,265 <= 264,265 + 268,269 | 266,267 <= 266,267 + 270,271
        vfadd(x_im.val[0], x_im.val[0], x_im.val[2]);
        vfadd(x_im.val[2], x_im.val[1], x_im.val[3]);
        vfadd(x_im.val[1], y_im.val[0], y_im.val[2]);
        vfadd(x_im.val[3], y_im.val[1], y_im.val[3]);

        // Calculate y
        // v_im*s_im
        vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[1], 0);
        vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[1], 0);
        vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[1], 1);
        vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[1], 1);

        // v_im*s_im + v_re*s_re
        // y_re: 4,5 | 6,7 | 12,13 | 14,15
        vfma_lane(y_re.val[0], tmp.val[0], v_re.val[0], s_re_im.val[0], 0);
        vfma_lane(y_re.val[1], tmp.val[1], v_re.val[1], s_re_im.val[0], 0);
        vfma_lane(y_re.val[2], tmp.val[2], v_re.val[2], s_re_im.val[0], 1);
        vfma_lane(y_re.val[3], tmp.val[3], v_re.val[3], s_re_im.val[0], 1);

        // v_im*s_re
        vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[0], 0);
        vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[0], 0);
        vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[0], 1);
        vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[0], 1);

        // v_im*s_re - v_re*s_im
        // y_im: 260,261 |  262,263 | 268,269 | 270,271
        vfms_lane(y_im.val[0], tmp.val[0], v_re.val[0], s_re_im.val[1], 0);
        vfms_lane(y_im.val[1], tmp.val[1], v_re.val[1], s_re_im.val[1], 0);
        vfms_lane(y_im.val[2], tmp.val[2], v_re.val[2], s_re_im.val[1], 1);
        vfms_lane(y_im.val[3], tmp.val[3], v_re.val[3], s_re_im.val[1], 1);

        // Level 3
        // x_re: 0,1 | 2,3 | 8,9 | 10,11
        // y_re: 4,5 | 6,7 | 12,13 | 14,15
        // x_im: 256,257 | 258,259 | 264,265 | 266,267
        // y_im: 260,261 | 262,263 | 268,269 | 270,271

        // Load s_re_im
        vload(s_re_im.val[0], &fpr_gm_tab[(FALCON_N + j) >> 3]);
        ////////
        // x - y
        ////////

        // 0,1 -   8,9 | 2,3 - 10,11
        // 4,5 - 12,13 | 6,7 - 14,15
        vfsubx4_swap(v_re, x_re, y_re, 0, 2, 1, 3);

        // 256,257 - 264,265 | 258,259 - 266,267
        // 260,261 - 268,269 | 262,263 - 270,271
        vfsubx4_swap(v_im, x_im, y_im, 0, 2, 1, 3);

        ////////
        // x + y
        ////////

        // x_re: 0,1 | 2,3 | 4,5 | 6,7
        // 0,1 <= 0,1 +   8,9 | 2,3 <= 2,3 + 10,11
        // 4,5 <= 4,5 + 12,13 | 6,7 <= 6,7 + 14,15
        vfaddx4_swap(x_re, x_re, y_re, 0, 2, 1, 3);

        // x_im: 256,257 | 258,259 | 260,261 | 262,263
        // 256,257 <= 256,257 + 264,265 | 258,259 <= 258,259 + 266,267
        // 260,261 <= 260,261 + 268,269 | 262,263 <= 262,263 + 270,271
        vfaddx4_swap(x_im, x_im, y_im, 0, 2, 1, 3);

        // Calculate y
        // v_im*s_im
        vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[0], 1);
        vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[0], 1);
        vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[0], 1);
        vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[0], 1);

        // v_im*s_im + v_re*s_re
        // y_re: 8,9 | 10,11 | 12,13 | 14,15
        vfma_lane(y_re.val[0], tmp.val[0], v_re.val[0], s_re_im.val[0], 0);
        vfma_lane(y_re.val[1], tmp.val[1], v_re.val[1], s_re_im.val[0], 0);
        vfma_lane(y_re.val[2], tmp.val[2], v_re.val[2], s_re_im.val[0], 0);
        vfma_lane(y_re.val[3], tmp.val[3], v_re.val[3], s_re_im.val[0], 0);

        // v_im*s_re
        vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[0], 0);
        vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[0], 0);
        vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[0], 0);
        vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[0], 0);

        // v_im*s_re - v_re*s_im
        // y_im: 264,265 | 266,267 | 268,269 | 270,271
        vfms_lane(y_im.val[0], tmp.val[0], v_re.val[0], s_re_im.val[0], 1);
        vfms_lane(y_im.val[1], tmp.val[1], v_re.val[1], s_re_im.val[0], 1);
        vfms_lane(y_im.val[2], tmp.val[2], v_re.val[2], s_re_im.val[0], 1);
        vfms_lane(y_im.val[3], tmp.val[3], v_re.val[3], s_re_im.val[0], 1);

        // x_re: 0,1 | 2,3 | 4,5 | 6,7
        // y_re: 8,9 | 10,11 | 12,13 | 14,15
        // x_im: 256,257 | 258,259 | 260,261 | 262,263
        // y_im: 264,265 | 266,267 | 268,269 | 270,271
        vstorex4(&f[j], x_re);
        vstorex4(&f[j + 8], y_re);
        vstorex4(&f[j + hn], x_im);
        vstorex4(&f[j + hn + 8], y_im);
    }

    // Level 5
    for (int i = 0; i < FALCON_N / 2; i += 1 << 6)
    {
        vloadx2(s_tmp, &fpr_gm_tab[(FALCON_N + i) >> 4]);
        s_re_im.val[0] = s_tmp.val[0];
        s_re_im.val[1] = s_tmp.val[1];
        vload(s_re_im.val[2], &fpr_gm_tab[(FALCON_N + i) >> 5]);

        for (int j = i; j < i + 16; j += 4)
        {
            // Layer 4
            // x_re: 0 ->3, 32->35
            // y_re: 16 -> 19, 48 -> 51
            // x_im: 256 -> 259, 288-> 291
            // y_im: 272 -> 275, 304 -> 307
            vloadx2(x_tmp, &f[j]);
            x_re.val[0] = x_tmp.val[0];
            x_re.val[1] = x_tmp.val[1];
            vloadx2(y_tmp, &f[j + 16]);
            y_re.val[0] = y_tmp.val[0];
            y_re.val[1] = y_tmp.val[1];

            vloadx2(x_tmp, &f[j + 32]);
            x_re.val[2] = x_tmp.val[0];
            x_re.val[3] = x_tmp.val[1];
            vloadx2(y_tmp, &f[j + 48]);
            y_re.val[2] = y_tmp.val[0];
            y_re.val[3] = y_tmp.val[1];

            vloadx2(x_tmp, &f[j + hn]);
            x_im.val[0] = x_tmp.val[0];
            x_im.val[1] = x_tmp.val[1];
            vloadx2(y_tmp, &f[j + hn + 16]);
            y_im.val[0] = y_tmp.val[0];
            y_im.val[1] = y_tmp.val[1];

            vloadx2(x_tmp, &f[j + hn + 32]);
            x_im.val[2] = x_tmp.val[0];
            x_im.val[3] = x_tmp.val[1];
            vloadx2(y_tmp, &f[j + hn + 48]);
            y_im.val[2] = y_tmp.val[0];
            y_im.val[3] = y_tmp.val[1];

            ////////
            // x - y
            ////////
            vfsubx4(v_re, x_re, y_re);
            vfsubx4(v_im, x_im, y_im);

            ////////
            // x + y
            ////////
            // x_re: 0->3, 32->35
            // x_im: 256->259, 288->291
            vfaddx4(x_re, x_re, y_re);
            vfaddx4(x_im, x_im, y_im);

            // Calculate y
            // v_im*s_im
            vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[0], 1);
            vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[0], 1);
            vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[1], 1);
            vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[1], 1);

            // v_im*s_im + v_re*s_re
            // y_re: 16->19, 48->51
            vfma_lane(y_re.val[0], tmp.val[0], v_re.val[0], s_re_im.val[0], 0);
            vfma_lane(y_re.val[1], tmp.val[1], v_re.val[1], s_re_im.val[0], 0);
            vfma_lane(y_re.val[2], tmp.val[2], v_re.val[2], s_re_im.val[1], 0);
            vfma_lane(y_re.val[3], tmp.val[3], v_re.val[3], s_re_im.val[1], 0);

            // v_im*s_re
            vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[0], 0);
            vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[0], 0);
            vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[1], 0);
            vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[1], 0);

            // v_im*s_re - v_re*s_im
            // y_im: 272->275, 304->307
            vfms_lane(y_im.val[0], tmp.val[0], v_re.val[0], s_re_im.val[0], 1);
            vfms_lane(y_im.val[1], tmp.val[1], v_re.val[1], s_re_im.val[0], 1);
            vfms_lane(y_im.val[2], tmp.val[2], v_re.val[2], s_re_im.val[1], 1);
            vfms_lane(y_im.val[3], tmp.val[3], v_re.val[3], s_re_im.val[1], 1);

            // Level 5:
            // x_re: 0->3, 32->35
            // y_re: 16->19, 48->51
            // x_im: 256->259, 288->291
            // y_im: 272->275, 304->307

            ////////
            // x - y
            ////////
            vfsubx4_swap(v_re, x_re, y_re, 0, 2, 1, 3);
            vfsubx4_swap(v_im, x_im, y_im, 0, 2, 1, 3);

            ////////
            // x + y
            ////////
            // x_re: 0->3, 16->19
            // x_im: 256->259, 272-> 275
            vfaddx4_swap(x_re, x_re, y_re, 0, 2, 1, 3);
            vfaddx4_swap(x_im, x_im, y_im, 0, 2, 1, 3);

            // Calculate y
            // v_im*s_im
            vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[2], 1);
            vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[2], 1);
            vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[2], 1);
            vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[2], 1);

            // v_im*s_im + v_re*s_re
            // y_re: 32->35, 48->51
            vfma_lane(y_re.val[0], tmp.val[0], v_re.val[0], s_re_im.val[2], 0);
            vfma_lane(y_re.val[1], tmp.val[1], v_re.val[1], s_re_im.val[2], 0);
            vfma_lane(y_re.val[2], tmp.val[2], v_re.val[2], s_re_im.val[2], 0);
            vfma_lane(y_re.val[3], tmp.val[3], v_re.val[3], s_re_im.val[2], 0);

            // v_im*s_re
            vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[2], 0);
            vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[2], 0);
            vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[2], 0);
            vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[2], 0);

            // v_im*s_re - v_re*s_im
            // y_im: 288->291, 304->307
            vfms_lane(y_im.val[0], tmp.val[0], v_re.val[0], s_re_im.val[2], 1);
            vfms_lane(y_im.val[1], tmp.val[1], v_re.val[1], s_re_im.val[2], 1);
            vfms_lane(y_im.val[2], tmp.val[2], v_re.val[2], s_re_im.val[2], 1);
            vfms_lane(y_im.val[3], tmp.val[3], v_re.val[3], s_re_im.val[2], 1);

            // Store
            x_tmp.val[0] = x_re.val[0];
            x_tmp.val[1] = x_re.val[1];
            vstorex2(&f[j], x_tmp);
            x_tmp.val[0] = x_re.val[2];
            x_tmp.val[1] = x_re.val[3];
            vstorex2(&f[j + 16], x_tmp);

            y_tmp.val[0] = y_re.val[0];
            y_tmp.val[1] = y_re.val[1];
            vstorex2(&f[j + 32], y_tmp);
            y_tmp.val[0] = y_re.val[2];
            y_tmp.val[1] = y_re.val[3];
            vstorex2(&f[j + 48], y_tmp);

            x_tmp.val[0] = x_im.val[0];
            x_tmp.val[1] = x_im.val[1];
            vstorex2(&f[j + hn], x_tmp);
            x_tmp.val[0] = x_im.val[2];
            x_tmp.val[1] = x_im.val[3];
            vstorex2(&f[j + hn + 16], x_tmp);

            y_tmp.val[0] = y_im.val[0];
            y_tmp.val[1] = y_im.val[1];
            vstorex2(&f[j + hn + 32], y_tmp);
            y_tmp.val[0] = y_im.val[2];
            y_tmp.val[1] = y_im.val[3];
            vstorex2(&f[j + hn + 48], y_tmp);
        }
    }

    vloadx2(s_tmp, &fpr_gm_tab[(FALCON_N) >> 6]);
    s_re_im.val[0] = s_tmp.val[0];
    s_re_im.val[1] = s_tmp.val[1];

    vload(s_re_im.val[2], &fpr_gm_tab[(FALCON_N) >> 7]);

    div_n = vdupq_n_f64(fpr_p2_tab[FALCON_LOGN]);
    vfmul(s_re_im.val[2], s_re_im.val[2], div_n);

    // Level 6
    for (int j = 0; j < 64; j += 4)
    {
        // Level 6:
        vloadx2(x_tmp, &f[j]);
        x_re.val[0] = x_tmp.val[0];
        x_re.val[1] = x_tmp.val[1];
        vloadx2(y_tmp, &f[j + 64]);
        y_re.val[0] = y_tmp.val[0];
        y_re.val[1] = y_tmp.val[1];

        vloadx2(x_tmp, &f[j + 128]);
        x_re.val[2] = x_tmp.val[0];
        x_re.val[3] = x_tmp.val[1];
        vloadx2(y_tmp, &f[j + 192]);
        y_re.val[2] = y_tmp.val[0];
        y_re.val[3] = y_tmp.val[1];

        vloadx2(x_tmp, &f[j + hn]);
        x_im.val[0] = x_tmp.val[0];
        x_im.val[1] = x_tmp.val[1];
        vloadx2(y_tmp, &f[j + hn + 64]);
        y_im.val[0] = y_tmp.val[0];
        y_im.val[1] = y_tmp.val[1];

        vloadx2(x_tmp, &f[j + hn + 128]);
        x_im.val[2] = x_tmp.val[0];
        x_im.val[3] = x_tmp.val[1];
        vloadx2(y_tmp, &f[j + hn + 192]);
        y_im.val[2] = y_tmp.val[0];
        y_im.val[3] = y_tmp.val[1];

        ////////
        // x - y
        ////////
        vfsubx4(v_re, x_re, y_re);
        vfsubx4(v_im, x_im, y_im);

        ////////
        // x + y
        ////////
        vfaddx4(x_re, x_re, y_re);
        vfaddx4(x_im, x_im, y_im);

        // Calculate y
        // v_im*s_im
        vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[0], 1);
        vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[0], 1);
        vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[1], 1);
        vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[1], 1);

        // v_im*s_im + v_re*s_re
        vfma_lane(y_re.val[0], tmp.val[0], v_re.val[0], s_re_im.val[0], 0);
        vfma_lane(y_re.val[1], tmp.val[1], v_re.val[1], s_re_im.val[0], 0);
        vfma_lane(y_re.val[2], tmp.val[2], v_re.val[2], s_re_im.val[1], 0);
        vfma_lane(y_re.val[3], tmp.val[3], v_re.val[3], s_re_im.val[1], 0);

        // v_im*s_re
        vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[0], 0);
        vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[0], 0);
        vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[1], 0);
        vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[1], 0);

        // v_im*s_re - v_re*s_im
        vfms_lane(y_im.val[0], tmp.val[0], v_re.val[0], s_re_im.val[0], 1);
        vfms_lane(y_im.val[1], tmp.val[1], v_re.val[1], s_re_im.val[0], 1);
        vfms_lane(y_im.val[2], tmp.val[2], v_re.val[2], s_re_im.val[1], 1);
        vfms_lane(y_im.val[3], tmp.val[3], v_re.val[3], s_re_im.val[1], 1);

        // Level 7:
        ////////
        // x - y
        ////////
        vfsubx4_swap(v_re, x_re, y_re, 0, 2, 1, 3);
        vfsubx4_swap(v_im, x_im, y_im, 0, 2, 1, 3);

        ////////
        // x + y
        ////////
        vfaddx4_swap(x_re, x_re, y_re, 0, 2, 1, 3);
        vfaddx4_swap(x_im, x_im, y_im, 0, 2, 1, 3);

        // Calculate y
        // v_im*s_im
        vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[2], 1);
        vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[2], 1);
        vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[2], 1);
        vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[2], 1);

        // v_im*s_im + v_re*s_re
        vfma_lane(y_re.val[0], tmp.val[0], v_re.val[0], s_re_im.val[2], 0);
        vfma_lane(y_re.val[1], tmp.val[1], v_re.val[1], s_re_im.val[2], 0);
        vfma_lane(y_re.val[2], tmp.val[2], v_re.val[2], s_re_im.val[2], 0);
        vfma_lane(y_re.val[3], tmp.val[3], v_re.val[3], s_re_im.val[2], 0);

        // v_im*s_re
        vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[2], 0);
        vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[2], 0);
        vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[2], 0);
        vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[2], 0);

        // v_im*s_re - v_re*s_im
        vfms_lane(y_im.val[0], tmp.val[0], v_re.val[0], s_re_im.val[2], 1);
        vfms_lane(y_im.val[1], tmp.val[1], v_re.val[1], s_re_im.val[2], 1);
        vfms_lane(y_im.val[2], tmp.val[2], v_re.val[2], s_re_im.val[2], 1);
        vfms_lane(y_im.val[3], tmp.val[3], v_re.val[3], s_re_im.val[2], 1);

        // Divide by N at the end
        vfmul(x_re.val[0], x_re.val[0], div_n);
        vfmul(x_re.val[1], x_re.val[1], div_n);
        vfmul(x_re.val[2], x_re.val[2], div_n);
        vfmul(x_re.val[3], x_re.val[3], div_n);

        vfmul(x_im.val[0], x_im.val[0], div_n);
        vfmul(x_im.val[1], x_im.val[1], div_n);
        vfmul(x_im.val[2], x_im.val[2], div_n);
        vfmul(x_im.val[3], x_im.val[3], div_n);

        // Store
        x_tmp.val[0] = x_re.val[0];
        x_tmp.val[1] = x_re.val[1];
        vstorex2(&f[j], x_tmp);
        x_tmp.val[0] = x_re.val[2];
        x_tmp.val[1] = x_re.val[3];
        vstorex2(&f[j + 64], x_tmp);

        y_tmp.val[0] = y_re.val[0];
        y_tmp.val[1] = y_re.val[1];
        vstorex2(&f[j + 128], y_tmp);
        y_tmp.val[0] = y_re.val[2];
        y_tmp.val[1] = y_re.val[3];
        vstorex2(&f[j + 192], y_tmp);

        x_tmp.val[0] = x_im.val[0];
        x_tmp.val[1] = x_im.val[1];
        vstorex2(&f[j + hn], x_tmp);
        x_tmp.val[0] = x_im.val[2];
        x_tmp.val[1] = x_im.val[3];
        vstorex2(&f[j + hn + 64], x_tmp);

        y_tmp.val[0] = y_im.val[0];
        y_tmp.val[1] = y_im.val[1];
        vstorex2(&f[j + hn + 128], y_tmp);
        y_tmp.val[0] = y_im.val[2];
        y_tmp.val[1] = y_im.val[3];
        vstorex2(&f[j + hn + 192], y_tmp);
    }
    // End function

    // Optional, combine two level 4-5, 6-7 loop,
    // but the compiler generate overhead loop, so I discard
    /* 
    int distance;
    const unsigned int LAST_L = 6; // Last level
    for (int l = 4; l < 8; l += 2)
    {
        distance = 1 << l;

        for (int i = 0; i < FALCON_N/2; i += 1 << (l+2) )
        {
            vloadx2(s_tmp, &fpr_gm_tab[(FALCON_N + i) >> l]);
            s_re_im.val[0] = s_tmp.val[0];
            s_re_im.val[1] = s_tmp.val[1];
            vload(s_re_im.val[2], &fpr_gm_tab[(FALCON_N + i) >> (l+1)]);

            if (l == LAST_L)
            {
                vfmul(s_re_im.val[2], s_re_im.val[2], div_n);
            }
            
            for (int j = i; j < i + distance; j += 4)
            {
                vloadx2(x_tmp, &f[j]);
                x_re.val[0] = x_tmp.val[0];
                x_re.val[1] = x_tmp.val[1];
                vloadx2(y_tmp, &f[j + distance]);
                y_re.val[0] = y_tmp.val[0];
                y_re.val[1] = y_tmp.val[1];

                vloadx2(x_tmp, &f[j + 2*distance]);
                x_re.val[2] = x_tmp.val[0];
                x_re.val[3] = x_tmp.val[1];
                vloadx2(y_tmp, &f[j + 3*distance]);
                y_re.val[2] = y_tmp.val[0];
                y_re.val[3] = y_tmp.val[1];

                vloadx2(x_tmp, &f[j + hn]);
                x_im.val[0] = x_tmp.val[0];
                x_im.val[1] = x_tmp.val[1];
                vloadx2(y_tmp, &f[j + hn + distance]);
                y_im.val[0] = y_tmp.val[0];
                y_im.val[1] = y_tmp.val[1];

                vloadx2(x_tmp, &f[j + hn + 2*distance]);
                x_im.val[2] = x_tmp.val[0];
                x_im.val[3] = x_tmp.val[1];
                vloadx2(y_tmp, &f[j + hn + 3*distance]);
                y_im.val[2] = y_tmp.val[0];
                y_im.val[3] = y_tmp.val[1];

                ////////
                // x - y
                ////////
                vfsubx4(v_re, x_re, y_re);
                vfsubx4(v_im, x_im, y_im);

                ////////
                // x + y
                ////////
                // x_re: 0->3, 32->35
                // x_im: 256->259, 288->291
                vfaddx4(x_re, x_re, y_re);
                vfaddx4(x_im, x_im, y_im);

                // Calculate y
                // v_im*s_im
                vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[0], 1);
                vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[0], 1);
                vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[1], 1);
                vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[1], 1);

                // v_im*s_im + v_re*s_re
                // y_re: 16->19, 48->51
                vfma_lane(y_re.val[0], tmp.val[0], v_re.val[0], s_re_im.val[0], 0);
                vfma_lane(y_re.val[1], tmp.val[1], v_re.val[1], s_re_im.val[0], 0);
                vfma_lane(y_re.val[2], tmp.val[2], v_re.val[2], s_re_im.val[1], 0);
                vfma_lane(y_re.val[3], tmp.val[3], v_re.val[3], s_re_im.val[1], 0);

                // v_im*s_re
                vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[0], 0);
                vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[0], 0);
                vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[1], 0);
                vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[1], 0);

                // v_im*s_re - v_re*s_im
                // y_im: 272->275, 304->307
                vfms_lane(y_im.val[0], tmp.val[0], v_re.val[0], s_re_im.val[0], 1);
                vfms_lane(y_im.val[1], tmp.val[1], v_re.val[1], s_re_im.val[0], 1);
                vfms_lane(y_im.val[2], tmp.val[2], v_re.val[2], s_re_im.val[1], 1);
                vfms_lane(y_im.val[3], tmp.val[3], v_re.val[3], s_re_im.val[1], 1);

                // Layer 6:
                // x_re: 0->3, 32->35
                // y_re: 16->19, 48->51
                // x_im: 256->259, 288->291
                // y_im: 272->275, 304->307

                ////////
                // x - y
                ////////
                vfsubx4_swap(v_re, x_re, y_re, 0, 2, 1, 3);
                vfsubx4_swap(v_im, x_im, y_im, 0, 2, 1, 3);

                ////////
                // x + y
                ////////
                // x_re: 0->3, 16->19
                // x_im: 256->259, 272-> 275
                vfaddx4_swap(x_re, x_re, y_re, 0, 2, 1, 3);
                vfaddx4_swap(x_im, x_im, y_im, 0, 2, 1, 3);

                // Calculate y
                // v_im*s_im
                vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[2], 1);
                vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[2], 1);
                vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[2], 1);
                vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[2], 1);

                // v_im*s_im + v_re*s_re
                // y_re: 32->35, 48->51
                vfma_lane(y_re.val[0], tmp.val[0], v_re.val[0], s_re_im.val[2], 0);
                vfma_lane(y_re.val[1], tmp.val[1], v_re.val[1], s_re_im.val[2], 0);
                vfma_lane(y_re.val[2], tmp.val[2], v_re.val[2], s_re_im.val[2], 0);
                vfma_lane(y_re.val[3], tmp.val[3], v_re.val[3], s_re_im.val[2], 0);

                // v_im*s_re
                vfmul_lane(tmp.val[0], v_im.val[0], s_re_im.val[2], 0);
                vfmul_lane(tmp.val[1], v_im.val[1], s_re_im.val[2], 0);
                vfmul_lane(tmp.val[2], v_im.val[2], s_re_im.val[2], 0);
                vfmul_lane(tmp.val[3], v_im.val[3], s_re_im.val[2], 0);

                // v_im*s_re - v_re*s_im
                // y_im: 288->291, 304->307
                vfms_lane(y_im.val[0], tmp.val[0], v_re.val[0], s_re_im.val[2], 1);
                vfms_lane(y_im.val[1], tmp.val[1], v_re.val[1], s_re_im.val[2], 1);
                vfms_lane(y_im.val[2], tmp.val[2], v_re.val[2], s_re_im.val[2], 1);
                vfms_lane(y_im.val[3], tmp.val[3], v_re.val[3], s_re_im.val[2], 1);

                if (l == LAST_L)
                {
                    vfmul(x_re.val[0], x_re.val[0], div_n);
                    vfmul(x_re.val[1], x_re.val[1], div_n);
                    vfmul(x_re.val[2], x_re.val[2], div_n);
                    vfmul(x_re.val[3], x_re.val[3], div_n);

                    vfmul(x_im.val[0], x_im.val[0], div_n);
                    vfmul(x_im.val[1], x_im.val[1], div_n);
                    vfmul(x_im.val[2], x_im.val[2], div_n);
                    vfmul(x_im.val[3], x_im.val[3], div_n);
                }

                // Store
                x_tmp.val[0] = x_re.val[0];
                x_tmp.val[1] = x_re.val[1];
                vstorex2(&f[j], x_tmp);
                x_tmp.val[0] = x_re.val[2];
                x_tmp.val[1] = x_re.val[3];
                vstorex2(&f[j + distance], x_tmp);

                y_tmp.val[0] = y_re.val[0];
                y_tmp.val[1] = y_re.val[1];
                vstorex2(&f[j + 2*distance], y_tmp);
                y_tmp.val[0] = y_re.val[2];
                y_tmp.val[1] = y_re.val[3];
                vstorex2(&f[j + 3*distance], y_tmp);

                x_tmp.val[0] = x_im.val[0];
                x_tmp.val[1] = x_im.val[1];
                vstorex2(&f[j + hn], x_tmp);
                x_tmp.val[0] = x_im.val[2];
                x_tmp.val[1] = x_im.val[3];
                vstorex2(&f[j + hn + distance], x_tmp);

                y_tmp.val[0] = y_im.val[0];
                y_tmp.val[1] = y_im.val[1];
                vstorex2(&f[j + hn + 2*distance], y_tmp);
                y_tmp.val[0] = y_im.val[2];
                y_tmp.val[1] = y_im.val[3];
                vstorex2(&f[j + hn + 3*distance], y_tmp);
            }
        }

    } */
}

void PQCLEAN_FALCON512_NEON_FFT(fpr *f, const bool negate_true)
{
    // Total: 32 = 16 + 8 + 8 register
    float64x2x4_t s_re_im, tmp;               // 8
    float64x2x4_t x_re, x_im, y_re, y_im;     // 16
    float64x2x4_t x1_re, x1_im, y1_re, y1_im; // 16
    float64x2x2_t s_tmp;                      // 2
    // Level 4, 5, 6, 7
    float64x2x2_t x_tmp, y_tmp;

    const unsigned int hn = FALCON_N >> 1;
    for (int l = 8; l > 4; l -= 2)
    {
        int distance = 1 << (l - 2);
        for (int i = 0; i < FALCON_N / 2; i += 1 << l)
        {
            vload(s_re_im.val[0], &fpr_gm_tab[(FALCON_N + i) >> (l - 1)]);
            vloadx2(s_tmp, &fpr_gm_tab[(FALCON_N + i) >> (l - 2)]);
            s_re_im.val[1] = s_tmp.val[0];
            s_re_im.val[2] = s_tmp.val[1];

            for (int j = i; j < i + distance; j += 4)
            {
                // Level 7
                // x1_re: 0->3, 64->67
                // x1_im: 256->259, 320 -> 323
                // y_re: 128->131, 192->195
                // y_im: 384->387, 448 -> 451
                vloadx2(x_tmp, &f[j]);
                x1_re.val[0] = x_tmp.val[0];
                x1_re.val[1] = x_tmp.val[1];
                vloadx2(x_tmp, &f[j + distance]);
                x1_re.val[2] = x_tmp.val[0];
                x1_re.val[3] = x_tmp.val[1];

                vloadx2(y_tmp, &f[j + 2 * distance]);
                y1_re.val[0] = y_tmp.val[0];
                y1_re.val[1] = y_tmp.val[1];
                vloadx2(y_tmp, &f[j + 3 * distance]);
                y1_re.val[2] = y_tmp.val[0];
                y1_re.val[3] = y_tmp.val[1];

                vloadx2(x_tmp, &f[j + hn]);
                x1_im.val[0] = x_tmp.val[0];
                x1_im.val[1] = x_tmp.val[1];
                vloadx2(x_tmp, &f[j + hn + distance]);
                x1_im.val[2] = x_tmp.val[0];
                x1_im.val[3] = x_tmp.val[1];

                vloadx2(y_tmp, &f[j + hn + 2 * distance]);
                y1_im.val[0] = y_tmp.val[0];
                y1_im.val[1] = y_tmp.val[1];
                vloadx2(y_tmp, &f[j + hn + 3 * distance]);
                y1_im.val[2] = y_tmp.val[0];
                y1_im.val[3] = y_tmp.val[1];

                vfmul_lane(y_re.val[0], y1_re.val[0], s_re_im.val[0], 0);
                vfmul_lane(y_re.val[1], y1_re.val[1], s_re_im.val[0], 0);
                vfmul_lane(y_re.val[2], y1_re.val[2], s_re_im.val[0], 0);
                vfmul_lane(y_re.val[3], y1_re.val[3], s_re_im.val[0], 0);

                vfms_lane(y_re.val[0], y_re.val[0], y1_im.val[0], s_re_im.val[0], 1);
                vfms_lane(y_re.val[1], y_re.val[1], y1_im.val[1], s_re_im.val[0], 1);
                vfms_lane(y_re.val[2], y_re.val[2], y1_im.val[2], s_re_im.val[0], 1);
                vfms_lane(y_re.val[3], y_re.val[3], y1_im.val[3], s_re_im.val[0], 1);

                vfmul_lane(y_im.val[0], y1_re.val[0], s_re_im.val[0], 1);
                vfmul_lane(y_im.val[1], y1_re.val[1], s_re_im.val[0], 1);
                vfmul_lane(y_im.val[2], y1_re.val[2], s_re_im.val[0], 1);
                vfmul_lane(y_im.val[3], y1_re.val[3], s_re_im.val[0], 1);

                vfma_lane(y_im.val[0], y_im.val[0], y1_im.val[0], s_re_im.val[0], 0);
                vfma_lane(y_im.val[1], y_im.val[1], y1_im.val[1], s_re_im.val[0], 0);
                vfma_lane(y_im.val[2], y_im.val[2], y1_im.val[2], s_re_im.val[0], 0);
                vfma_lane(y_im.val[3], y_im.val[3], y1_im.val[3], s_re_im.val[0], 0);

                vfaddx4(x_re, x1_re, y_re);
                vfaddx4(x_im, x1_im, y_im);

                vfsubx4(y_re, x1_re, y_re);
                vfsubx4(y_im, x1_im, y_im);

                // Level 6
                // x_re: 0->3, 64->67
                // y_re: 128->131, 192 -> 195
                // x_im: 256->259, 320 -> 323
                // y_im: 384->387, 448 -> 451

                vfmul_lane(y1_re.val[0], x_re.val[2], s_re_im.val[1], 0);
                vfmul_lane(y1_re.val[1], x_re.val[3], s_re_im.val[1], 0);
                vfmul_lane(y1_re.val[2], y_re.val[2], s_re_im.val[2], 0);
                vfmul_lane(y1_re.val[3], y_re.val[3], s_re_im.val[2], 0);

                vfms_lane(y1_re.val[0], y1_re.val[0], x_im.val[2], s_re_im.val[1], 1);
                vfms_lane(y1_re.val[1], y1_re.val[1], x_im.val[3], s_re_im.val[1], 1);
                vfms_lane(y1_re.val[2], y1_re.val[2], y_im.val[2], s_re_im.val[2], 1);
                vfms_lane(y1_re.val[3], y1_re.val[3], y_im.val[3], s_re_im.val[2], 1);

                vfmul_lane(y1_im.val[0], x_re.val[2], s_re_im.val[1], 1);
                vfmul_lane(y1_im.val[1], x_re.val[3], s_re_im.val[1], 1);
                vfmul_lane(y1_im.val[2], y_re.val[2], s_re_im.val[2], 1);
                vfmul_lane(y1_im.val[3], y_re.val[3], s_re_im.val[2], 1);

                vfma_lane(y1_im.val[0], y1_im.val[0], x_im.val[2], s_re_im.val[1], 0);
                vfma_lane(y1_im.val[1], y1_im.val[1], x_im.val[3], s_re_im.val[1], 0);
                vfma_lane(y1_im.val[2], y1_im.val[2], y_im.val[2], s_re_im.val[2], 0);
                vfma_lane(y1_im.val[3], y1_im.val[3], y_im.val[3], s_re_im.val[2], 0);

                vfadd(x1_re.val[0], x_re.val[0], y1_re.val[0]);
                vfadd(x1_re.val[1], x_re.val[1], y1_re.val[1]);
                vfadd(x1_re.val[2], y_re.val[0], y1_re.val[2]);
                vfadd(x1_re.val[3], y_re.val[1], y1_re.val[3]);

                vfadd(x1_im.val[0], x_im.val[0], y1_im.val[0]);
                vfadd(x1_im.val[1], x_im.val[1], y1_im.val[1]);
                vfadd(x1_im.val[2], y_im.val[0], y1_im.val[2]);
                vfadd(x1_im.val[3], y_im.val[1], y1_im.val[3]);

                vfsub(y1_re.val[0], x_re.val[0], y1_re.val[0]);
                vfsub(y1_re.val[1], x_re.val[1], y1_re.val[1]);
                vfsub(y1_re.val[2], y_re.val[0], y1_re.val[2]);
                vfsub(y1_re.val[3], y_re.val[1], y1_re.val[3]);

                vfsub(y1_im.val[0], x_im.val[0], y1_im.val[0]);
                vfsub(y1_im.val[1], x_im.val[1], y1_im.val[1]);
                vfsub(y1_im.val[2], y_im.val[0], y1_im.val[2]);
                vfsub(y1_im.val[3], y_im.val[1], y1_im.val[3]);

                // Level 6
                // x1_re: 0->3, 128 -> 131
                // y1_re: 64->67, 192 -> 195
                // x1_im: 256 -> 259, 384->387
                // y1_im: 320->323, 448 -> 451

                // Store
                x_tmp.val[0] = x1_re.val[0];
                x_tmp.val[1] = x1_re.val[1];
                vstorex2(&f[j], x_tmp);
                y_tmp.val[0] = y1_re.val[0];
                y_tmp.val[1] = y1_re.val[1];
                vstorex2(&f[j + distance], y_tmp);

                x_tmp.val[0] = x1_re.val[2];
                x_tmp.val[1] = x1_re.val[3];
                vstorex2(&f[j + 2 * distance], x_tmp);
                y_tmp.val[0] = y1_re.val[2];
                y_tmp.val[1] = y1_re.val[3];
                vstorex2(&f[j + 3 * distance], y_tmp);

                x_tmp.val[0] = x1_im.val[0];
                x_tmp.val[1] = x1_im.val[1];
                vstorex2(&f[j + hn], x_tmp);
                y_tmp.val[0] = y1_im.val[0];
                y_tmp.val[1] = y1_im.val[1];
                vstorex2(&f[j + hn + distance], y_tmp);

                x_tmp.val[0] = x1_im.val[2];
                x_tmp.val[1] = x1_im.val[3];
                vstorex2(&f[j + hn + 2 * distance], x_tmp);
                y_tmp.val[0] = y1_im.val[2];
                y_tmp.val[1] = y1_im.val[3];
                vstorex2(&f[j + hn + 3 * distance], y_tmp);
            }
        }
    }
    // End level 7, 6, 5, 4 loop

    // Level 3, 2, 1, 0
    for (int j = 0; j < FALCON_N / 2; j += 16)
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

        vload(s_re_im.val[0], &fpr_gm_tab[(FALCON_N + j) >> 3]);

        vfmul_lane(y1_re.val[0], y_re.val[0], s_re_im.val[0], 0);
        vfmul_lane(y1_re.val[1], y_re.val[1], s_re_im.val[0], 0);
        vfmul_lane(y1_re.val[2], y_re.val[2], s_re_im.val[0], 0);
        vfmul_lane(y1_re.val[3], y_re.val[3], s_re_im.val[0], 0);

        vfms_lane(y1_re.val[0], y1_re.val[0], y_im.val[0], s_re_im.val[0], 1);
        vfms_lane(y1_re.val[1], y1_re.val[1], y_im.val[1], s_re_im.val[0], 1);
        vfms_lane(y1_re.val[2], y1_re.val[2], y_im.val[2], s_re_im.val[0], 1);
        vfms_lane(y1_re.val[3], y1_re.val[3], y_im.val[3], s_re_im.val[0], 1);

        vfmul_lane(y1_im.val[0], y_re.val[0], s_re_im.val[0], 1);
        vfmul_lane(y1_im.val[1], y_re.val[1], s_re_im.val[0], 1);
        vfmul_lane(y1_im.val[2], y_re.val[2], s_re_im.val[0], 1);
        vfmul_lane(y1_im.val[3], y_re.val[3], s_re_im.val[0], 1);

        vfma_lane(y1_im.val[0], y1_im.val[0], y_im.val[0], s_re_im.val[0], 0);
        vfma_lane(y1_im.val[1], y1_im.val[1], y_im.val[1], s_re_im.val[0], 0);
        vfma_lane(y1_im.val[2], y1_im.val[2], y_im.val[2], s_re_im.val[0], 0);
        vfma_lane(y1_im.val[3], y1_im.val[3], y_im.val[3], s_re_im.val[0], 0);

        vfaddx4(x1_re, x_re, y1_re);
        vfaddx4(x1_im, x_im, y1_im);

        vfsubx4(y1_re, x_re, y1_re);
        vfsubx4(y1_im, x_im, y1_im);

        // Level 2
        // x_re: 0->7
        // y_re: 8->15
        // x_im: 256->263
        // y_im: 264->271
        vloadx2(s_tmp, &fpr_gm_tab[(FALCON_N + j) >> 2]);
        s_re_im.val[0] = s_tmp.val[0];
        s_re_im.val[1] = s_tmp.val[1];

        vfmul_lane(y_re.val[0], x1_re.val[2], s_re_im.val[0], 0);
        vfmul_lane(y_re.val[1], x1_re.val[3], s_re_im.val[0], 0);
        vfmul_lane(y_re.val[2], y1_re.val[2], s_re_im.val[1], 0);
        vfmul_lane(y_re.val[3], y1_re.val[3], s_re_im.val[1], 0);

        vfms_lane(y_re.val[0], y_re.val[0], x1_im.val[2], s_re_im.val[0], 1);
        vfms_lane(y_re.val[1], y_re.val[1], x1_im.val[3], s_re_im.val[0], 1);
        vfms_lane(y_re.val[2], y_re.val[2], y1_im.val[2], s_re_im.val[1], 1);
        vfms_lane(y_re.val[3], y_re.val[3], y1_im.val[3], s_re_im.val[1], 1);

        vfmul_lane(y_im.val[0], x1_re.val[2], s_re_im.val[0], 1);
        vfmul_lane(y_im.val[1], x1_re.val[3], s_re_im.val[0], 1);
        vfmul_lane(y_im.val[2], y1_re.val[2], s_re_im.val[1], 1);
        vfmul_lane(y_im.val[3], y1_re.val[3], s_re_im.val[1], 1);

        vfma_lane(y_im.val[0], y_im.val[0], x1_im.val[2], s_re_im.val[0], 0);
        vfma_lane(y_im.val[1], y_im.val[1], x1_im.val[3], s_re_im.val[0], 0);
        vfma_lane(y_im.val[2], y_im.val[2], y1_im.val[2], s_re_im.val[1], 0);
        vfma_lane(y_im.val[3], y_im.val[3], y1_im.val[3], s_re_im.val[1], 0);

        vfadd(x_re.val[0], x1_re.val[0], y_re.val[0]);
        vfadd(x_re.val[1], x1_re.val[1], y_re.val[1]);
        vfadd(x_re.val[2], y1_re.val[0], y_re.val[2]);
        vfadd(x_re.val[3], y1_re.val[1], y_re.val[3]);

        vfadd(x_im.val[0], x1_im.val[0], y_im.val[0]);
        vfadd(x_im.val[1], x1_im.val[1], y_im.val[1]);
        vfadd(x_im.val[2], y1_im.val[0], y_im.val[2]);
        vfadd(x_im.val[3], y1_im.val[1], y_im.val[3]);

        vfsub(y_re.val[0], x1_re.val[0], y_re.val[0]);
        vfsub(y_re.val[1], x1_re.val[1], y_re.val[1]);
        vfsub(y_re.val[2], y1_re.val[0], y_re.val[2]);
        vfsub(y_re.val[3], y1_re.val[1], y_re.val[3]);

        vfsub(y_im.val[0], x1_im.val[0], y_im.val[0]);
        vfsub(y_im.val[1], x1_im.val[1], y_im.val[1]);
        vfsub(y_im.val[2], y1_im.val[0], y_im.val[2]);
        vfsub(y_im.val[3], y1_im.val[1], y_im.val[3]);

        // Level 1
        // x_re: 0->3, 8->11
        // y_re: 4->7, 12->15
        // x_im: 256->259, 264->267
        // y_im: 260->263, 268->271
        vloadx4(s_re_im, &fpr_gm_tab[(FALCON_N + j) >> 1]);

        vfmul_lane(y1_re.val[0], x_re.val[1], s_re_im.val[0], 0);
        vfmul_lane(y1_re.val[1], y_re.val[1], s_re_im.val[1], 0);
        vfmul_lane(y1_re.val[2], x_re.val[3], s_re_im.val[2], 0);
        vfmul_lane(y1_re.val[3], y_re.val[3], s_re_im.val[3], 0);

        vfms_lane(y1_re.val[0], y1_re.val[0], x_im.val[1], s_re_im.val[0], 1);
        vfms_lane(y1_re.val[1], y1_re.val[1], y_im.val[1], s_re_im.val[1], 1);
        vfms_lane(y1_re.val[2], y1_re.val[2], x_im.val[3], s_re_im.val[2], 1);
        vfms_lane(y1_re.val[3], y1_re.val[3], y_im.val[3], s_re_im.val[3], 1);

        vfmul_lane(y1_im.val[0], x_re.val[1], s_re_im.val[0], 1);
        vfmul_lane(y1_im.val[1], y_re.val[1], s_re_im.val[1], 1);
        vfmul_lane(y1_im.val[2], x_re.val[3], s_re_im.val[2], 1);
        vfmul_lane(y1_im.val[3], y_re.val[3], s_re_im.val[3], 1);

        vfma_lane(y1_im.val[0], y1_im.val[0], x_im.val[1], s_re_im.val[0], 0);
        vfma_lane(y1_im.val[1], y1_im.val[1], y_im.val[1], s_re_im.val[1], 0);
        vfma_lane(y1_im.val[2], y1_im.val[2], x_im.val[3], s_re_im.val[2], 0);
        vfma_lane(y1_im.val[3], y1_im.val[3], y_im.val[3], s_re_im.val[3], 0);

        vfadd(x1_re.val[0], x_re.val[0], y1_re.val[0]);
        vfadd(x1_re.val[1], y_re.val[0], y1_re.val[1]);
        vfadd(x1_re.val[2], x_re.val[2], y1_re.val[2]);
        vfadd(x1_re.val[3], y_re.val[2], y1_re.val[3]);

        vfadd(x1_im.val[0], x_im.val[0], y1_im.val[0]);
        vfadd(x1_im.val[1], y_im.val[0], y1_im.val[1]);
        vfadd(x1_im.val[2], x_im.val[2], y1_im.val[2]);
        vfadd(x1_im.val[3], y_im.val[2], y1_im.val[3]);

        vfsub(y1_re.val[0], x_re.val[0], y1_re.val[0]);
        vfsub(y1_re.val[1], y_re.val[0], y1_re.val[1]);
        vfsub(y1_re.val[2], x_re.val[2], y1_re.val[2]);
        vfsub(y1_re.val[3], y_re.val[2], y1_re.val[3]);

        vfsub(y1_im.val[0], x_im.val[0], y1_im.val[0]);
        vfsub(y1_im.val[1], y_im.val[0], y1_im.val[1]);
        vfsub(y1_im.val[2], x_im.val[2], y1_im.val[2]);
        vfsub(y1_im.val[3], y_im.val[2], y1_im.val[3]);

        // Level 0
        // Before Transpose
        // x_re: 0,1 | 4,5 | 8,9   | 12,13
        // y_re: 2,3 | 6,7 | 10,11 | 14,15
        // x_im: 256,257 | 260,261 | 264,265 | 268,269
        // y_im: 258,259 | 262,263 | 266,267 | 270,271
        transpose(x1_re, x1_re, tmp, 0, 1, 0);
        transpose(x1_re, x1_re, tmp, 2, 3, 1);
        transpose(y1_re, y1_re, tmp, 0, 1, 2);
        transpose(y1_re, y1_re, tmp, 2, 3, 3);

        transpose(x1_im, x1_im, tmp, 0, 1, 0);
        transpose(x1_im, x1_im, tmp, 2, 3, 1);
        transpose(y1_im, y1_im, tmp, 0, 1, 2);
        transpose(y1_im, y1_im, tmp, 2, 3, 3);
        // After Transpose
        // x_re: 0,4 | 1,5 | 8,12  | 9,13
        // y_re: 2,6 | 3,7 | 10,14 | 11,15
        // x_im: 256,260 | 257,261 | 264,268 | 265,269
        // y_im: 258,262 | 259,263 | 266,270 | 267,271
        vload4(s_re_im, &fpr_gm_tab[FALCON_N + j]);

        vfmul(y_re.val[0], x1_re.val[1], s_re_im.val[0]);
        vfmul(y_re.val[1], y1_re.val[1], s_re_im.val[2]);
        vfms(y_re.val[0], y_re.val[0], x1_im.val[1], s_re_im.val[1]);
        vfms(y_re.val[1], y_re.val[1], y1_im.val[1], s_re_im.val[3]);

        vfmul(y_im.val[0], x1_re.val[1], s_re_im.val[1]);
        vfmul(y_im.val[1], y1_re.val[1], s_re_im.val[3]);
        vfma(y_im.val[0], y_im.val[0], x1_im.val[1], s_re_im.val[0]);
        vfma(y_im.val[1], y_im.val[1], y1_im.val[1], s_re_im.val[2]);

        vload4(s_re_im, &fpr_gm_tab[FALCON_N + j + 8]);

        vfmul(y_re.val[2], x1_re.val[3], s_re_im.val[0]);
        vfmul(y_re.val[3], y1_re.val[3], s_re_im.val[2]);
        vfms(y_re.val[2], y_re.val[2], x1_im.val[3], s_re_im.val[1]);
        vfms(y_re.val[3], y_re.val[3], y1_im.val[3], s_re_im.val[3]);

        vfmul(y_im.val[2], x1_re.val[3], s_re_im.val[1]);
        vfmul(y_im.val[3], y1_re.val[3], s_re_im.val[3]);
        vfma(y_im.val[2], y_im.val[2], x1_im.val[3], s_re_im.val[0]);
        vfma(y_im.val[3], y_im.val[3], y1_im.val[3], s_re_im.val[2]);

        vfadd(x_re.val[0], x1_re.val[0], y_re.val[0]);
        vfadd(x_re.val[1], y1_re.val[0], y_re.val[1]);
        vfadd(x_re.val[2], x1_re.val[2], y_re.val[2]);
        vfadd(x_re.val[3], y1_re.val[2], y_re.val[3]);

        vfadd(x_im.val[0], x1_im.val[0], y_im.val[0]);
        vfadd(x_im.val[1], y1_im.val[0], y_im.val[1]);
        vfadd(x_im.val[2], x1_im.val[2], y_im.val[2]);
        vfadd(x_im.val[3], y1_im.val[2], y_im.val[3]);

        // Constant propagation will optimize this loop
        if (negate_true)
        {
            vfsub(y_re.val[0], y_re.val[0], x1_re.val[0]);
            vfsub(y_re.val[1], y_re.val[1], y1_re.val[0]);

            vfsub(y_re.val[2], y_re.val[2], x1_re.val[2]);
            vfsub(y_re.val[3], y_re.val[3], y1_re.val[2]);

            vfsub(y_im.val[0], y_im.val[0], x1_im.val[0]);
            vfsub(y_im.val[1], y_im.val[1], y1_im.val[0]);

            vfsub(y_im.val[2], y_im.val[2], x1_im.val[2]);
            vfsub(y_im.val[3], y_im.val[3], y1_im.val[2]);

            vfneg(x_re.val[0], x_re.val[0]);
            vfneg(x_re.val[1], x_re.val[1]);
            vfneg(x_re.val[2], x_re.val[2]);
            vfneg(x_re.val[3], x_re.val[3]);

            vfneg(x_im.val[0], x_im.val[0]);
            vfneg(x_im.val[1], x_im.val[1]);
            vfneg(x_im.val[2], x_im.val[2]);
            vfneg(x_im.val[3], x_im.val[3]);
        }
        else
        {
            vfsub(y_re.val[0], x1_re.val[0], y_re.val[0]);
            vfsub(y_re.val[1], y1_re.val[0], y_re.val[1]);

            vfsub(y_re.val[2], x1_re.val[2], y_re.val[2]);
            vfsub(y_re.val[3], y1_re.val[2], y_re.val[3]);

            vfsub(y_im.val[0], x1_im.val[0], y_im.val[0]);
            vfsub(y_im.val[1], y1_im.val[0], y_im.val[1]);

            vfsub(y_im.val[2], x1_im.val[2], y_im.val[2]);
            vfsub(y_im.val[3], y1_im.val[2], y_im.val[3]);
        }

        // x_re: 0,4 | 2,6 | 8,12 | 10,14
        // y_re: 1,5 | 3,7 | 9,13 | 11,15
        // x_im: 256,260 | 258,262 | 264,268 | 266,270
        // y_im: 257,261 | 259,263 | 265,269 | 267,271
        x1_re.val[0] = x_re.val[0];
        x1_re.val[1] = y_re.val[0];
        x1_re.val[2] = x_re.val[1];
        x1_re.val[3] = y_re.val[1];

        y1_re.val[0] = x_re.val[2];
        y1_re.val[1] = y_re.val[2];
        y1_re.val[2] = x_re.val[3];
        y1_re.val[3] = y_re.val[3];

        x1_im.val[0] = x_im.val[0];
        x1_im.val[1] = y_im.val[0];
        x1_im.val[2] = x_im.val[1];
        x1_im.val[3] = y_im.val[1];

        y1_im.val[0] = x_im.val[2];
        y1_im.val[1] = y_im.val[2];
        y1_im.val[2] = x_im.val[3];
        y1_im.val[3] = y_im.val[3];

        vstore4(&f[j], x1_re);
        vstore4(&f[j + 8], y1_re);
        vstore4(&f[j + hn], x1_im);
        vstore4(&f[j + hn + 8], y1_im);
    }
    // End function
}

// /////////////////////////========================

/*
 * Addition of two complex numbers (d = a + b).
 */
#define FPC_ADD(d_re, d_im, a_re, a_im, b_re, b_im) \
    do                                              \
    {                                               \
        fpr fpct_re, fpct_im;                       \
        fpct_re = fpr_add(a_re, b_re);              \
        fpct_im = fpr_add(a_im, b_im);              \
        (d_re) = fpct_re;                           \
        (d_im) = fpct_im;                           \
    } while (0)

/*
 * Subtraction of two complex numbers (d = a - b).
 */
#define FPC_SUB(d_re, d_im, a_re, a_im, b_re, b_im) \
    do                                              \
    {                                               \
        fpr fpct_re, fpct_im;                       \
        fpct_re = fpr_sub(a_re, b_re);              \
        fpct_im = fpr_sub(a_im, b_im);              \
        (d_re) = fpct_re;                           \
        (d_im) = fpct_im;                           \
    } while (0)

/*
 * Multplication of two complex numbers (d = a * b).
 */
#define FPC_MUL(d_re, d_im, a_re, a_im, b_re, b_im) \
    do                                              \
    {                                               \
        fpr fpct_a_re, fpct_a_im;                   \
        fpr fpct_b_re, fpct_b_im;                   \
        fpr fpct_d_re, fpct_d_im;                   \
        fpct_a_re = (a_re);                         \
        fpct_a_im = (a_im);                         \
        fpct_b_re = (b_re);                         \
        fpct_b_im = (b_im);                         \
        fpct_d_re = fpr_sub(                        \
            fpr_mul(fpct_a_re, fpct_b_re),          \
            fpr_mul(fpct_a_im, fpct_b_im));         \
        fpct_d_im = fpr_add(                        \
            fpr_mul(fpct_a_re, fpct_b_im),          \
            fpr_mul(fpct_a_im, fpct_b_re));         \
        (d_re) = fpct_d_re;                         \
        (d_im) = fpct_d_im;                         \
    } while (0)

/*
 * Squaring of a complex number (d = a * a).
 */
#define FPC_SQR(d_re, d_im, a_re, a_im)                              \
    do                                                               \
    {                                                                \
        fpr fpct_a_re, fpct_a_im;                                    \
        fpr fpct_d_re, fpct_d_im;                                    \
        fpct_a_re = (a_re);                                          \
        fpct_a_im = (a_im);                                          \
        fpct_d_re = fpr_sub(fpr_sqr(fpct_a_re), fpr_sqr(fpct_a_im)); \
        fpct_d_im = fpr_double(fpr_mul(fpct_a_re, fpct_a_im));       \
        (d_re) = fpct_d_re;                                          \
        (d_im) = fpct_d_im;                                          \
    } while (0)

/*
 * Inversion of a complex number (d = 1 / a).
 */
#define FPC_INV(d_re, d_im, a_re, a_im)                           \
    do                                                            \
    {                                                             \
        fpr fpct_a_re, fpct_a_im;                                 \
        fpr fpct_d_re, fpct_d_im;                                 \
        fpr fpct_m;                                               \
        fpct_a_re = (a_re);                                       \
        fpct_a_im = (a_im);                                       \
        fpct_m = fpr_add(fpr_sqr(fpct_a_re), fpr_sqr(fpct_a_im)); \
        fpct_m = fpr_inv(fpct_m);                                 \
        fpct_d_re = fpr_mul(fpct_a_re, fpct_m);                   \
        fpct_d_im = fpr_mul(fpr_neg(fpct_a_im), fpct_m);          \
        (d_re) = fpct_d_re;                                       \
        (d_im) = fpct_d_im;                                       \
    } while (0)

/*
 * Division of complex numbers (d = a / b).
 */
#define FPC_DIV(d_re, d_im, a_re, a_im, b_re, b_im)               \
    do                                                            \
    {                                                             \
        fpr fpct_a_re, fpct_a_im;                                 \
        fpr fpct_b_re, fpct_b_im;                                 \
        fpr fpct_d_re, fpct_d_im;                                 \
        fpr fpct_m;                                               \
        fpct_a_re = (a_re);                                       \
        fpct_a_im = (a_im);                                       \
        fpct_b_re = (b_re);                                       \
        fpct_b_im = (b_im);                                       \
        fpct_m = fpr_add(fpr_sqr(fpct_b_re), fpr_sqr(fpct_b_im)); \
        fpct_m = fpr_inv(fpct_m);                                 \
        fpct_b_re = fpr_mul(fpct_b_re, fpct_m);                   \
        fpct_b_im = fpr_mul(fpr_neg(fpct_b_im), fpct_m);          \
        fpct_d_re = fpr_sub(                                      \
            fpr_mul(fpct_a_re, fpct_b_re),                        \
            fpr_mul(fpct_a_im, fpct_b_im));                       \
        fpct_d_im = fpr_add(                                      \
            fpr_mul(fpct_a_re, fpct_b_im),                        \
            fpr_mul(fpct_a_im, fpct_b_re));                       \
        (d_re) = fpct_d_re;                                       \
        (d_im) = fpct_d_im;                                       \
    } while (0)

void PQCLEAN_FALCON512_CLEAN_FFT_original(fpr *f, unsigned logn)
{
    /*
     * FFT algorithm in bit-reversal order uses the following
     * iterative algorithm:
     *
     *   t = N
     *   for m = 1; m < N; m *= 2:
     *       ht = t/2
     *       for i1 = 0; i1 < m; i1 ++:
     *           j1 = i1 * t
     *           s = GM[m + i1]
     *           for j = j1; j < (j1 + ht); j ++:
     *               x = f[j]
     *               y = s * f[j + ht]
     *               f[j] = x + y
     *               f[j + ht] = x - y
     *       t = ht
     *
     * GM[k] contains w^rev(k) for primitive root w = exp(i*pi/N).
     *
     * In the description above, f[] is supposed to contain complex
     * numbers. In our in-memory representation, the real and
     * imaginary parts of f[k] are in array slots k and k+N/2.
     *
     * We only keep the first half of the complex numbers. We can
     * see that after the first iteration, the first and second halves
     * of the array of complex numbers have separate lives, so we
     * simply ignore the second part.
     */

    unsigned u;
    size_t t, n, hn, m;

    /*
     * First iteration: compute f[j] + i * f[j+N/2] for all j < N/2
     * (because GM[1] = w^rev(1) = w^(N/2) = i).
     * In our chosen representation, this is a no-op: everything is
     * already where it should be.
     */

    /*
     * Subsequent iterations are truncated to use only the first
     * half of values.
     */
    n = (size_t)1 << logn;
    hn = n >> 1;
    t = hn;
    for (u = 1, m = 2; u < logn; u++, m <<= 1)
    {
        size_t ht, hm, i1, j1;

        ht = t >> 1;
        hm = m >> 1;
        for (i1 = 0, j1 = 0; i1 < hm; i1++, j1 += t)
        {
            size_t j, j2;

            j2 = j1 + ht;
            fpr s_re, s_im;

            s_re = fpr_gm_tab[((m + i1) << 1) + 0];
            s_im = fpr_gm_tab[((m + i1) << 1) + 1];
            for (j = j1; j < j2; j++)
            {
                fpr x_re, x_im, y_re, y_im;

                x_re = f[j];
                x_im = f[j + hn];
                y_re = f[j + ht];
                y_im = f[j + ht + hn];
                FPC_MUL(y_re, y_im, y_re, y_im, s_re, s_im);
                FPC_ADD(f[j], f[j + hn],
                        x_re, x_im, y_re, y_im);
                FPC_SUB(f[j + ht], f[j + ht + hn],
                        x_re, x_im, y_re, y_im);
            }
        }
        t = ht;
    }
}

static void Zf(FFT_log2)(fpr *f)
{
    /* 
    x_re:   0 =   0 + (  1*  4 -   3*  5)
    x_im:   2 =   2 + (  1*  5 +   3*  4)
    y_re:   1 =   0 - (  1*  4 -   3*  5)
    y_im:   3 =   2 - (  1*  5 +   3*  4)
    */
    float64x2x2_t tmp, tmp2;
    float64x2_t s_re_im, s_re_im_rev, neon_i21, v;
    const double imagine[2] = {-1.0, 1.0};
    vload(neon_i21, &imagine[0]);

    // re: 0, 2
    // im: 1, 3
    vload2(tmp, &f[0]);
    vload(s_re_im, &fpr_gm_tab[4]);

    // -5, 4
    s_re_im_rev = vextq_f64(s_re_im, s_re_im, 1);
    vfmul(s_re_im_rev, s_re_im_rev, neon_i21);
    vfmul_lane(v, s_re_im, tmp.val[1], 0);
    vfma_lane(v, v, s_re_im_rev, tmp.val[1], 1);

    vfadd(tmp2.val[0], tmp.val[0], v);
    vfsub(tmp2.val[1], tmp.val[0], v);

    vstore2(&f[0], tmp2);
}

static void Zf(FFT_log3)(fpr *f)
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
    vfms_lane(v_re, v_re, tmp.val[3], s_re_im.val[0], 1);

    vfmul_lane(v_im, tmp.val[1], s_re_im.val[0], 1);
    vfma_lane(v_im, v_im, tmp.val[3], s_re_im.val[0], 0);

    vfsub(tmp.val[1], tmp.val[0], v_re);
    vfadd(tmp.val[0], tmp.val[0], v_re);

    vfsub(tmp.val[3], tmp.val[2], v_im);
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

static void Zf(FFT_log4)(fpr *f)
{
    // Total SIMD register: 30 = 12 + 16 + 2
    float64x2x4_t tmp[2], s_re_im;                             // 12
    float64x2x2_t s_tmp, x_re, x_im, y_re, y_im, v1, v2, tmp2; // 16
    float64x2_t v_re, v_im;                                    // 2

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

    vloadx4(tmp[0], &f[0]);
    vloadx4(tmp[1], &f[8]);
    s_re_im.val[0] = vld1q_f64(&fpr_gm_tab[4]);

    vfmul_lane(v1.val[0], tmp[0].val[2], s_re_im.val[0], 0);
    vfmul_lane(v1.val[1], tmp[0].val[3], s_re_im.val[0], 0);
    vfms_lane(v1.val[0], v1.val[0], tmp[1].val[2], s_re_im.val[0], 1);
    vfms_lane(v1.val[1], v1.val[1], tmp[1].val[3], s_re_im.val[0], 1);

    vfmul_lane(v2.val[0], tmp[0].val[2], s_re_im.val[0], 1);
    vfmul_lane(v2.val[1], tmp[0].val[3], s_re_im.val[0], 1);
    vfma_lane(v2.val[0], v2.val[0], tmp[1].val[2], s_re_im.val[0], 0);
    vfma_lane(v2.val[1], v2.val[1], tmp[1].val[3], s_re_im.val[0], 0);

    vfsub(tmp[0].val[2], tmp[0].val[0], v1.val[0]);
    vfsub(tmp[0].val[3], tmp[0].val[1], v1.val[1]);
    vfadd(tmp[0].val[0], tmp[0].val[0], v1.val[0]);
    vfadd(tmp[0].val[1], tmp[0].val[1], v1.val[1]);

    vfsub(tmp[1].val[2], tmp[1].val[0], v2.val[0]);
    vfsub(tmp[1].val[3], tmp[1].val[1], v2.val[1]);
    vfadd(tmp[1].val[0], tmp[1].val[0], v2.val[0]);
    vfadd(tmp[1].val[1], tmp[1].val[1], v2.val[1]);

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

    vfmul_lane(v1.val[0], tmp[0].val[1], s_tmp.val[0], 0);
    vfmul_lane(v1.val[1], tmp[0].val[3], s_tmp.val[1], 0);
    vfms_lane(v1.val[0], v1.val[0], tmp[1].val[1], s_tmp.val[0], 1);
    vfms_lane(v1.val[1], v1.val[1], tmp[1].val[3], s_tmp.val[1], 1);

    vfmul_lane(v2.val[0], tmp[0].val[1], s_tmp.val[0], 1);
    vfmul_lane(v2.val[1], tmp[0].val[3], s_tmp.val[1], 1);
    vfma_lane(v2.val[0], v2.val[0], tmp[1].val[1], s_tmp.val[0], 0);
    vfma_lane(v2.val[1], v2.val[1], tmp[1].val[3], s_tmp.val[1], 0);

    vfsub(y_re.val[0], tmp[0].val[0], v1.val[0]);
    vfsub(y_re.val[1], tmp[0].val[2], v1.val[1]);
    vfadd(x_re.val[0], tmp[0].val[0], v1.val[0]);
    vfadd(x_re.val[1], tmp[0].val[2], v1.val[1]);

    vfsub(y_im.val[0], tmp[1].val[0], v2.val[0]);
    vfsub(y_im.val[1], tmp[1].val[2], v2.val[1]);
    vfadd(x_im.val[0], tmp[1].val[0], v2.val[0]);
    vfadd(x_im.val[1], tmp[1].val[2], v2.val[1]);

    // x_re: 0, 1 | 4, 5
    // y_re: 2, 3 | 6, 7
    // x_im: 8, 9 | 12, 13
    // y_im: 10, 11 | 14, 15

    tmp[0].val[0] = vtrn1q_f64(x_re.val[0], y_re.val[0]);
    tmp[0].val[1] = vtrn2q_f64(x_re.val[0], y_re.val[0]);
    tmp[0].val[2] = vtrn1q_f64(x_re.val[1], y_re.val[1]);
    tmp[0].val[3] = vtrn2q_f64(x_re.val[1], y_re.val[1]);

    tmp[1].val[0] = vtrn1q_f64(x_im.val[0], y_im.val[0]);
    tmp[1].val[1] = vtrn2q_f64(x_im.val[0], y_im.val[0]);
    tmp[1].val[2] = vtrn1q_f64(x_im.val[1], y_im.val[1]);
    tmp[1].val[3] = vtrn2q_f64(x_im.val[1], y_im.val[1]);

    // tmp[0]: 0, 2 | 1, 3 |  4, 6 |  5, 7
    // tmp[1]: 8,10 | 9,11 | 12,14 | 13,15
    /* 
    level 3
    y_re:   1 =   0 - (  1* 16 -   9* 17)
    y_re:   3 =   2 - (  3* 18 -  11* 19)
    y_re:   5 =   4 - (  5* 20 -  13* 21)
    y_re:   7 =   6 - (  7* 22 -  15* 23)

    x_re:   0 =   0 + (  1* 16 -   9* 17)
    x_re:   2 =   2 + (  3* 18 -  11* 19)
    x_re:   4 =   4 + (  5* 20 -  13* 21)
    x_re:   6 =   6 + (  7* 22 -  15* 23)

    y_im:   9 =   8 - (  1* 17 +   9* 16)
    y_im:  11 =  10 - (  3* 19 +  11* 18)
    y_im:  13 =  12 - (  5* 21 +  13* 20)
    y_im:  15 =  14 - (  7* 23 +  15* 22)

    x_im:   8 =   8 + (  1* 17 +   9* 16)
    x_im:  10 =  10 + (  3* 19 +  11* 18)
    x_im:  12 =  12 + (  5* 21 +  13* 20)
    x_im:  14 =  14 + (  7* 23 +  15* 22)
     */
    // s_re_im: 16,18 | 17,19 | 20,22 | 21,23
    vload2(s_tmp, &fpr_gm_tab[16]);
    s_re_im.val[0] = s_tmp.val[0];
    s_re_im.val[1] = s_tmp.val[1];
    vload2(s_tmp, &fpr_gm_tab[20]);
    s_re_im.val[2] = s_tmp.val[0];
    s_re_im.val[3] = s_tmp.val[1];

    vfmul(v1.val[0], tmp[0].val[1], s_re_im.val[0]);
    vfmul(v1.val[1], tmp[0].val[3], s_re_im.val[2]);
    vfms(v1.val[0], v1.val[0], tmp[1].val[1], s_re_im.val[1]);
    vfms(v1.val[1], v1.val[1], tmp[1].val[3], s_re_im.val[3]);

    vfmul(v2.val[0], tmp[0].val[1], s_re_im.val[1]);
    vfmul(v2.val[1], tmp[0].val[3], s_re_im.val[3]);
    vfma(v2.val[0], v2.val[0], tmp[1].val[1], s_re_im.val[0]);
    vfma(v2.val[1], v2.val[1], tmp[1].val[3], s_re_im.val[2]);

    vfsub(y_re.val[0], tmp[0].val[0], v1.val[0]);
    vfsub(y_re.val[1], tmp[0].val[2], v1.val[1]);
    vfadd(x_re.val[0], tmp[0].val[0], v1.val[0]);
    vfadd(x_re.val[1], tmp[0].val[2], v1.val[1]);

    vfsub(y_im.val[0], tmp[1].val[0], v2.val[0]);
    vfsub(y_im.val[1], tmp[1].val[2], v2.val[1]);
    vfadd(x_im.val[0], tmp[1].val[0], v2.val[0]);
    vfadd(x_im.val[1], tmp[1].val[2], v2.val[1]);

    tmp2.val[0] = x_re.val[0];
    tmp2.val[1] = y_re.val[0];
    vstore2(&f[0], tmp2);
    tmp2.val[0] = x_re.val[1];
    tmp2.val[1] = y_re.val[1];
    vstore2(&f[4], tmp2);
    tmp2.val[0] = x_im.val[0];
    tmp2.val[1] = y_im.val[0];
    vstore2(&f[8], tmp2);
    tmp2.val[0] = x_im.val[1];
    tmp2.val[1] = y_im.val[1];
    vstore2(&f[12], tmp2);
}

static void Zf(FFT_log5)(fpr *f, unsigned logn)
{
    // Total: 34 = 16 + 8 + 2 register
    float64x2x4_t s_re_im, tmp;                   // 8
    float64x2x4_t x_re, x_im, y_re, y_im, v1, v2; // 24
    float64x2x2_t s_tmp;                          // 2
    // Level 4, 5, 6, 7
    float64x2x2_t x_tmp, y_tmp;
    const unsigned int falcon_n = 1 << logn;
    const unsigned int hn = falcon_n >> 1;

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

        vload(s_re_im.val[0], &fpr_gm_tab[(falcon_n + j) >> 3]);

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
        vloadx2(s_tmp, &fpr_gm_tab[(falcon_n + j) >> 2]);
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
        vloadx4(s_re_im, &fpr_gm_tab[(falcon_n + j) >> 1]);

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
        vload4(s_re_im, &fpr_gm_tab[falcon_n + j]);

        vfmul(v1.val[0], x_re.val[1], s_re_im.val[0]);
        vfmul(v1.val[1], y_re.val[1], s_re_im.val[2]);
        vfms(v1.val[0], v1.val[0], x_im.val[1], s_re_im.val[1]);
        vfms(v1.val[1], v1.val[1], y_im.val[1], s_re_im.val[3]);

        vfmul(v2.val[0], x_re.val[1], s_re_im.val[1]);
        vfmul(v2.val[1], y_re.val[1], s_re_im.val[3]);
        vfma(v2.val[0], v2.val[0], x_im.val[1], s_re_im.val[0]);
        vfma(v2.val[1], v2.val[1], y_im.val[1], s_re_im.val[2]);

        vload4(s_re_im, &fpr_gm_tab[falcon_n + j + 8]);

        vfmul(v1.val[2], x_re.val[3], s_re_im.val[0]);
        vfmul(v1.val[3], y_re.val[3], s_re_im.val[2]);
        vfms(v1.val[2], v1.val[2], x_im.val[3], s_re_im.val[1]);
        vfms(v1.val[3], v1.val[3], y_im.val[3], s_re_im.val[3]);

        vfmul(v2.val[2], x_re.val[3], s_re_im.val[1]);
        vfmul(v2.val[3], y_re.val[3], s_re_im.val[3]);
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

static void Zf(FFT_logn1)(fpr *f, unsigned logn, const unsigned l)
{
    // Total SIMD register: 26
    float64x2_t s_re_im;                          // 2
    float64x2x4_t x_re, x_im, y_re, y_im, v1, v2; // 24

    const unsigned int falcon_n = 1 << logn;
    const unsigned int hn = falcon_n >> 1;
    const unsigned level = l - 1;
    const unsigned distance = 1 << (level - 1);
    for (unsigned i = 0; i < hn; i += 1 << level)
    {
        vload(s_re_im, &fpr_gm_tab[(falcon_n + i) >> (level - 1)]);
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

static void Zf(FFT_logn2)(fpr *f, unsigned logn, const unsigned level)
{
    // Total: 26 = 16 + 8 + 2 register
    float64x2x4_t s_re_im, tmp;                   // 8
    float64x2x4_t x_re, x_im, y_re, y_im, v1, v2; // 16
    float64x2x2_t s_tmp;                          // 2
    // for storing instruction
    float64x2x2_t x_tmp, y_tmp;

    const unsigned int falcon_n = 1 << logn;
    const unsigned int hn = falcon_n >> 1;

    for (int l = level - 1; l > 4; l -= 2)
    {
        int distance = 1 << (l - 2);
        for (int i = 0; i < hn; i += 1 << l)
        {
            vload(s_re_im.val[0], &fpr_gm_tab[(falcon_n + i) >> (l - 1)]);
            vloadx2(s_tmp, &fpr_gm_tab[(falcon_n + i) >> (l - 2)]);
            s_re_im.val[1] = s_tmp.val[0];
            s_re_im.val[2] = s_tmp.val[1];

            for (int j = i; j < i + distance; j += 4)
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

                // Level 6
                // x_re: 0->3, 64->67
                // y_re: 128->131, 192 -> 195
                // x_im: 256->259, 320 -> 323
                // y_im: 384->387, 448 -> 451

                vfmul_lane(v1.val[0], x_re.val[2], s_re_im.val[1], 0);
                vfmul_lane(v1.val[1], x_re.val[3], s_re_im.val[1], 0);
                vfmul_lane(v1.val[2], y_re.val[2], s_re_im.val[2], 0);
                vfmul_lane(v1.val[3], y_re.val[3], s_re_im.val[2], 0);

                vfms_lane(v1.val[0], v1.val[0], x_im.val[2], s_re_im.val[1], 1);
                vfms_lane(v1.val[1], v1.val[1], x_im.val[3], s_re_im.val[1], 1);
                vfms_lane(v1.val[2], v1.val[2], y_im.val[2], s_re_im.val[2], 1);
                vfms_lane(v1.val[3], v1.val[3], y_im.val[3], s_re_im.val[2], 1);

                vfmul_lane(v2.val[0], x_re.val[2], s_re_im.val[1], 1);
                vfmul_lane(v2.val[1], x_re.val[3], s_re_im.val[1], 1);
                vfmul_lane(v2.val[2], y_re.val[2], s_re_im.val[2], 1);
                vfmul_lane(v2.val[3], y_re.val[3], s_re_im.val[2], 1);

                vfma_lane(v2.val[0], v2.val[0], x_im.val[2], s_re_im.val[1], 0);
                vfma_lane(v2.val[1], v2.val[1], x_im.val[3], s_re_im.val[1], 0);
                vfma_lane(v2.val[2], v2.val[2], y_im.val[2], s_re_im.val[2], 0);
                vfma_lane(v2.val[3], v2.val[3], y_im.val[3], s_re_im.val[2], 0);

                // Level 6
                // x1_re: 0->3, 128 -> 131
                // y1_re: 64->67, 192 -> 195
                // x1_im: 256 -> 259, 384->387
                // y1_im: 320->323, 448 -> 451
                vfsub(y_tmp.val[0], x_re.val[0], v1.val[0]);
                vfsub(y_tmp.val[1], x_re.val[1], v1.val[1]);
                vfadd(x_tmp.val[0], x_re.val[0], v1.val[0]);
                vfadd(x_tmp.val[1], x_re.val[1], v1.val[1]);

                vstorex2(&f[j], x_tmp);
                vstorex2(&f[j + distance], y_tmp);

                vfsub(y_tmp.val[0], y_re.val[0], v1.val[2]);
                vfsub(y_tmp.val[1], y_re.val[1], v1.val[3]);
                vfadd(x_tmp.val[0], y_re.val[0], v1.val[2]);
                vfadd(x_tmp.val[1], y_re.val[1], v1.val[3]);

                vstorex2(&f[j + 2 * distance], x_tmp);
                vstorex2(&f[j + 3 * distance], y_tmp);

                vfsub(y_tmp.val[0], x_im.val[0], v2.val[0]);
                vfsub(y_tmp.val[1], x_im.val[1], v2.val[1]);
                vfadd(x_tmp.val[0], x_im.val[0], v2.val[0]);
                vfadd(x_tmp.val[1], x_im.val[1], v2.val[1]);

                vstorex2(&f[j + hn], x_tmp);
                vstorex2(&f[j + hn + distance], y_tmp);

                vfsub(y_tmp.val[0], y_im.val[0], v2.val[2]);
                vfsub(y_tmp.val[1], y_im.val[1], v2.val[3]);
                vfadd(x_tmp.val[0], y_im.val[0], v2.val[2]);
                vfadd(x_tmp.val[1], y_im.val[1], v2.val[3]);

                vstorex2(&f[j + hn + 2 * distance], x_tmp);
                vstorex2(&f[j + hn + 3 * distance], y_tmp);
            }
        }
    }
}
/* 
 * Support logn from [1, 10]
 * Can be easily extended to logn > 10
 */
void Zf(FFT_logn)(fpr *f, const unsigned logn)
{
    unsigned level = logn;
    switch (logn)
    {
    case 1:
        break;

    case 2:
        Zf(FFT_log2)(f);
        break;

    case 3:
        Zf(FFT_log3)(f);
        break;

    case 4:
        Zf(FFT_log4)(f);
        break;

    case 5:
        Zf(FFT_log5)(f, logn);
        break;

    case 6:
        Zf(FFT_logn1)(f, logn, level--);
        Zf(FFT_log5)(f, logn);
        break;

    case 7:
    case 9:
        Zf(FFT_logn2)(f, logn, level);
        Zf(FFT_log5)(f, logn);
        break;

    default:
        // case 8:
        // case 10:
        Zf(FFT_logn1)(f, logn, level--);
        Zf(FFT_logn2)(f, logn, level);
        Zf(FFT_log5)(f, logn);
        break;
    }
}