/*
 * poly_int.c
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

#include <arm_neon.h>
#include "macrous.h"
#include "config.h"
#include "poly.h"
#include "ntt_consts.h"

void ZfN(poly_smallints_to_bigints)(int16_t out[FALCON_N], const int8_t in[FALCON_N])
{
    // Total SIMD registers: 24 = 16 + 8
    int16x8x4_t a, b, e, f; // 16
    int8x16x4_t c, d;       // 8

    for (int i = 0; i < FALCON_N; i += 128)
    {
        c = vld1q_s8_x4(&in[i]);
        d = vld1q_s8_x4(&in[i + 64]);

        a.val[0] = vmovl_s8(vget_low_s8(c.val[0]));
        a.val[2] = vmovl_s8(vget_low_s8(c.val[1]));
        b.val[0] = vmovl_s8(vget_low_s8(c.val[2]));
        b.val[2] = vmovl_s8(vget_low_s8(c.val[3]));

        a.val[1] = vmovl_high_s8(c.val[0]);
        a.val[3] = vmovl_high_s8(c.val[1]);
        b.val[1] = vmovl_high_s8(c.val[2]);
        b.val[3] = vmovl_high_s8(c.val[3]);

        e.val[0] = vmovl_s8(vget_low_s8(d.val[0]));
        e.val[2] = vmovl_s8(vget_low_s8(d.val[1]));
        f.val[0] = vmovl_s8(vget_low_s8(d.val[2]));
        f.val[2] = vmovl_s8(vget_low_s8(d.val[3]));

        e.val[1] = vmovl_high_s8(d.val[0]);
        e.val[3] = vmovl_high_s8(d.val[1]);
        f.val[1] = vmovl_high_s8(d.val[2]);
        f.val[3] = vmovl_high_s8(d.val[3]);

        vst1q_s16_x4(&out[i], a);
        vst1q_s16_x4(&out[i + 32], b);
        vst1q_s16_x4(&out[i + 64], e);
        vst1q_s16_x4(&out[i + 96], f);
    }
}

/*
 * Return f[] = f[]/g[] % 12289
 * See assembly https://godbolt.org/z/od3Ex7Mbx
 */

void ZfN(poly_div_12289)(int16_t f[FALCON_N], const int16_t g[FALCON_N])
{
    // Total SIMD registers: 24 = 4 + 19 + 1
    int16x8x4_t src, dst, t, k; // 4
    int16x8x4_t y0, y1, y2, y3, y4, y5,
        y6, y7, y8, y9, y10, y11, y12,
        y13, y14, y15, y16, y17, y18; // 19
    int16x8_t neon_qmvm;              // 1

    neon_qmvm = vld1q_s16(qmvq);

    for (int i = 0; i < FALCON_N; i += 32)
    {
        // Find y0 = g^12287
        vload_s16_x4(y0, &g[i]);
        vload_s16_x4(src, &f[i]);

        // y0 is already in Montgomery domain

        montmul_x4(y1, y0, y0, neon_qmvm, t);
        montmul_x4(y2, y1, y0, neon_qmvm, k);
        montmul_x4(y3, y2, y1, neon_qmvm, t);
        montmul_x4(y4, y3, y3, neon_qmvm, k);
        montmul_x4(y5, y4, y4, neon_qmvm, t);
        montmul_x4(y6, y5, y5, neon_qmvm, k);
        montmul_x4(y7, y6, y6, neon_qmvm, t);
        montmul_x4(y8, y7, y7, neon_qmvm, k);
        montmul_x4(y9, y8, y2, neon_qmvm, t);
        montmul_x4(y10, y9, y8, neon_qmvm, k);
        montmul_x4(y11, y10, y10, neon_qmvm, t);
        montmul_x4(y12, y11, y11, neon_qmvm, k);
        montmul_x4(y13, y12, y9, neon_qmvm, t);
        montmul_x4(y14, y13, y13, neon_qmvm, k);
        montmul_x4(y15, y14, y14, neon_qmvm, t);
        montmul_x4(y16, y15, y10, neon_qmvm, k);
        montmul_x4(y17, y16, y16, neon_qmvm, t);
        montmul_x4(y18, y17, y0, neon_qmvm, k);
        montmul_x4(dst, y18, src, neon_qmvm, t);

        vstore_s16_x4(&f[i], dst);
    }
}


/* 
 * f = g - s
 */
void ZfN(poly_sub_barrett)(int16_t f[FALCON_N], const int16_t g[FALCON_N], const int16_t s[FALCON_N])
{
    // Total SIMD registers: 29 = 28 + 1
    int16x8x4_t a, b, c, d, e, h, t; // 28
    int16x8_t neon_qmvm;             // 1
    neon_qmvm = vld1q_s16(qmvq);

    for (int i = 0; i < FALCON_N; i += 64)
    {
        vload_s16_x4(a, &s[i]);
        vload_s16_x4(b, &g[i]);
        vload_s16_x4(c, &s[i + 32]);
        vload_s16_x4(d, &g[i + 32]);

        e.val[0] = vsubq_s16(a.val[0], b.val[0]);
        e.val[1] = vsubq_s16(a.val[1], b.val[1]);
        e.val[2] = vsubq_s16(a.val[2], b.val[2]);
        e.val[3] = vsubq_s16(a.val[3], b.val[3]);

        h.val[0] = vsubq_s16(c.val[0], d.val[0]);
        h.val[1] = vsubq_s16(c.val[1], d.val[1]);
        h.val[2] = vsubq_s16(c.val[2], d.val[2]);
        h.val[3] = vsubq_s16(c.val[3], d.val[3]);

        barrett_x4(e, neon_qmvm, t);
        barrett_x4(h, neon_qmvm, t);

        vstore_s16_x4(&f[i], e);
        vstore_s16_x4(&f[i + 32], h);
    }
}

/*
 * Check f[] has 0
 * Return:
 * 1 if 0 in f[]
 * otherwise, 0
 */
uint16_t ZfN(poly_compare_with_zero)(int16_t f[FALCON_N])
{
    // Total SIMD registers: 22 = 12 + 8 + 2
    int16x8x4_t a, b;      // 8
    uint16x8x4_t c, d, e1; // 12
    uint16x8x2_t e2;       // 2

    e2.val[1] = vdupq_n_u16(0);

    for (int i = 0; i < FALCON_N; i += 64)
    {
        vload_s16_x4(a, &f[i]);
        vload_s16_x4(b, &f[i + 32]);

        // Compare bitwise Equal to zero (vector)
        // a == 0 ? 1 : 0;
        c.val[0] = vceqzq_s16(a.val[0]);
        c.val[1] = vceqzq_s16(a.val[1]);
        c.val[2] = vceqzq_s16(a.val[2]);
        c.val[3] = vceqzq_s16(a.val[3]);

        d.val[0] = vceqzq_s16(b.val[0]);
        d.val[1] = vceqzq_s16(b.val[1]);
        d.val[2] = vceqzq_s16(b.val[2]);
        d.val[3] = vceqzq_s16(b.val[3]);

        e1.val[0] = vorrq_u16(d.val[0], c.val[0]);
        e1.val[1] = vorrq_u16(d.val[1], c.val[1]);
        e1.val[2] = vorrq_u16(d.val[2], c.val[2]);
        e1.val[3] = vorrq_u16(d.val[3], c.val[3]);

        e1.val[0] = vorrq_u16(e1.val[0], e1.val[2]);
        e1.val[1] = vorrq_u16(e1.val[1], e1.val[3]);

        e2.val[0] = vorrq_u16(e1.val[0], e1.val[1]);

        e2.val[1] = vorrq_u16(e2.val[1], e2.val[0]);
    }

    uint16_t ret = vmaxvq_u16(e2.val[1]);

    return ret;
}

/*
 * Branchless conditional addtion with FALCON_Q if coeffcient is < 0
 */
void ZfN(poly_convert_to_unsigned)(int16_t f[FALCON_N])
{
    // Total SIMD registers: 25 = 8 + 16 + 1
    uint16x8x4_t b0, b1;        // 8
    int16x8x4_t a0, a1, c0, c1; // 16
    uint16x8_t neon_q;          // 1

    neon_q = vdupq_n_u16(FALCON_Q);

    for (int i = 0; i < FALCON_N; i += 64)
    {
        vload_s16_x4(a0, &f[i]);
        vload_s16_x4(a1, &f[i + 32]);

        b0.val[0] = vcltzq_s16(a0.val[0]);
        b0.val[1] = vcltzq_s16(a0.val[1]);
        b0.val[2] = vcltzq_s16(a0.val[2]);
        b0.val[3] = vcltzq_s16(a0.val[3]);

        b1.val[0] = vcltzq_s16(a1.val[0]);
        b1.val[1] = vcltzq_s16(a1.val[1]);
        b1.val[2] = vcltzq_s16(a1.val[2]);
        b1.val[3] = vcltzq_s16(a1.val[3]);

        c0.val[0] = (int16x8_t)vandq_u16(b0.val[0], neon_q);
        c0.val[1] = (int16x8_t)vandq_u16(b0.val[1], neon_q);
        c0.val[2] = (int16x8_t)vandq_u16(b0.val[2], neon_q);
        c0.val[3] = (int16x8_t)vandq_u16(b0.val[3], neon_q);

        c1.val[0] = (int16x8_t)vandq_u16(b1.val[0], neon_q);
        c1.val[1] = (int16x8_t)vandq_u16(b1.val[1], neon_q);
        c1.val[2] = (int16x8_t)vandq_u16(b1.val[2], neon_q);
        c1.val[3] = (int16x8_t)vandq_u16(b1.val[3], neon_q);

        vadd_x4(c0, a0, c0);
        vadd_x4(c1, a1, c1);

        vstore_s16_x4(&f[i], c0);
        vstore_s16_x4(&f[i + 32], c1);
    }
}

int ZfN(bigints_to_smallints)(int8_t G[FALCON_N], const int16_t t[FALCON_N])
{
    // Total SIMD registers: 32
    int16x8x4_t a, f;              // 8
    uint16x8x4_t c0, c1, d0, d1;   // 16
    uint16x8x2_t e;                // 2
    int8x16x4_t g;                 // 4
    int16x8_t neon_127, neon__127; // 2
    neon_127 = vdupq_n_s16(127);
    neon__127 = vdupq_n_s16(-127);

    e.val[1] = vdupq_n_u16(0);

    for (int i = 0; i < FALCON_N; i += 64)
    {
        vload_s16_x4(a, &t[i]);
        vload_s16_x4(f, &t[i + 32]);

        g.val[0] = vmovn_high_s16(vmovn_s16(a.val[0]), a.val[1]);
        g.val[1] = vmovn_high_s16(vmovn_s16(a.val[2]), a.val[3]);
        g.val[2] = vmovn_high_s16(vmovn_s16(f.val[0]), f.val[1]);
        g.val[3] = vmovn_high_s16(vmovn_s16(f.val[2]), f.val[3]);

        vst1q_s8_x4(&G[i], g);

        // -127 > a ? 1 : 0
        c0.val[0] = vcgtq_s16(neon__127, a.val[0]);
        c0.val[1] = vcgtq_s16(neon__127, a.val[1]);
        c0.val[2] = vcgtq_s16(neon__127, a.val[2]);
        c0.val[3] = vcgtq_s16(neon__127, a.val[3]);
        // a > 127 ? 1 : 0
        c1.val[0] = vcgtq_s16(a.val[0], neon_127);
        c1.val[1] = vcgtq_s16(a.val[1], neon_127);
        c1.val[2] = vcgtq_s16(a.val[2], neon_127);
        c1.val[3] = vcgtq_s16(a.val[3], neon_127);

        // -127 > f ? 1 : 0
        d0.val[0] = vcgtq_s16(neon__127, f.val[0]);
        d0.val[1] = vcgtq_s16(neon__127, f.val[1]);
        d0.val[2] = vcgtq_s16(neon__127, f.val[2]);
        d0.val[3] = vcgtq_s16(neon__127, f.val[3]);
        // f > 127 ? 1 : 0
        d1.val[0] = vcgtq_s16(f.val[0], neon_127);
        d1.val[1] = vcgtq_s16(f.val[1], neon_127);
        d1.val[2] = vcgtq_s16(f.val[2], neon_127);
        d1.val[3] = vcgtq_s16(f.val[3], neon_127);

        c0.val[0] = vorrq_u16(c0.val[0], c1.val[0]);
        c0.val[1] = vorrq_u16(c0.val[1], c1.val[1]);
        c0.val[2] = vorrq_u16(c0.val[2], c1.val[2]);
        c0.val[3] = vorrq_u16(c0.val[3], c1.val[3]);

        d0.val[0] = vorrq_u16(d0.val[0], d1.val[0]);
        d0.val[1] = vorrq_u16(d0.val[1], d1.val[1]);
        d0.val[2] = vorrq_u16(d0.val[2], d1.val[2]);
        d0.val[3] = vorrq_u16(d0.val[3], d1.val[3]);

        c0.val[0] = vorrq_u16(c0.val[0], d0.val[0]);
        c0.val[2] = vorrq_u16(c0.val[2], d0.val[2]);
        c0.val[1] = vorrq_u16(c0.val[1], d0.val[1]);
        c0.val[3] = vorrq_u16(c0.val[3], d0.val[3]);

        c0.val[0] = vorrq_u16(c0.val[0], c0.val[2]);
        c0.val[1] = vorrq_u16(c0.val[1], c0.val[3]);

        e.val[0] = vorrq_u16(c0.val[0], c0.val[1]);

        e.val[1] = vorrq_u16(e.val[1], e.val[0]);
    }
    if (vmaxvq_u16(e.val[1]))
    {
        return 1;
    }
    return 0;
}


/*
 * Check if (t < low || t > high)
 * Return 1 if True
 * Otherwise 0
 */
int ZfN(poly_check_bound_int8)(const int8_t t[FALCON_N], 
                        const int8_t low, const int8_t high)
{
    // Total SIMD registers: 10
    int8x16x4_t a;                 // 4
    uint8x16x4_t c;                // 4
    uint8x16_t e;                  // 1
    int8x16_t neon_low, neon_high; // 1

    neon_high = vdupq_n_s8(high);
    neon_low = vdupq_n_s8(low);
    e = vdupq_n_u8(0);

    for (int i = 0; i < FALCON_N; i += 64)
    {
        a = vld1q_s8_x4(&t[i]);

        // low > a ? 1 : 0
        c.val[0] = vcgtq_s8(neon_low, a.val[0]);
        c.val[1] = vcgtq_s8(neon_low, a.val[1]);
        c.val[2] = vcgtq_s8(neon_low, a.val[2]);
        c.val[3] = vcgtq_s8(neon_low, a.val[3]);
        // a > high ? 1 : 0
        c.val[0] = vcgtq_s8(a.val[0], neon_high);
        c.val[1] = vcgtq_s8(a.val[1], neon_high);
        c.val[2] = vcgtq_s8(a.val[2], neon_high);
        c.val[3] = vcgtq_s8(a.val[3], neon_high);

        c.val[0] = vorrq_u8(c.val[0], c.val[0]);
        c.val[1] = vorrq_u8(c.val[1], c.val[1]);
        c.val[2] = vorrq_u8(c.val[2], c.val[2]);
        c.val[3] = vorrq_u8(c.val[3], c.val[3]);

        c.val[0] = vorrq_u8(c.val[0], c.val[2]);
        c.val[1] = vorrq_u8(c.val[1], c.val[3]);

        c.val[0] = vorrq_u8(c.val[0], c.val[1]);

        e = vorrq_u8(e, c.val[0]);

        if (vmaxvq_u8(e))
        {
            return 1;
        }
    }
    return 0;
}

/*
 * Check if (t < low || t > high)
 * Return 1 if True
 * Otherwise 0
 */
int ZfN(poly_check_bound_int16)(const int16_t t[FALCON_N], 
                    const int16_t low, const int16_t high)
{
    // Total SIMD registers: 10
    int16x8x4_t a;                 // 4
    uint16x8x4_t c;                // 4
    uint16x8_t e;                  // 1
    int16x8_t neon_low, neon_high; // 1

    neon_high = vdupq_n_s16(high);
    neon_low = vdupq_n_s16(low);
    e = vdupq_n_u16(0);

    for (int i = 0; i < FALCON_N; i += 32)
    {
        a = vld1q_s16_x4(&t[i]);

        // low > a ? 1 : 0
        c.val[0] = vcgtq_s16(neon_low, a.val[0]);
        c.val[1] = vcgtq_s16(neon_low, a.val[1]);
        c.val[2] = vcgtq_s16(neon_low, a.val[2]);
        c.val[3] = vcgtq_s16(neon_low, a.val[3]);
        // a > high ? 1 : 0
        c.val[0] = vcgtq_s16(a.val[0], neon_high);
        c.val[1] = vcgtq_s16(a.val[1], neon_high);
        c.val[2] = vcgtq_s16(a.val[2], neon_high);
        c.val[3] = vcgtq_s16(a.val[3], neon_high);

        c.val[0] = vorrq_u16(c.val[0], c.val[0]);
        c.val[1] = vorrq_u16(c.val[1], c.val[1]);
        c.val[2] = vorrq_u16(c.val[2], c.val[2]);
        c.val[3] = vorrq_u16(c.val[3], c.val[3]);

        c.val[0] = vorrq_u16(c.val[0], c.val[2]);
        c.val[1] = vorrq_u16(c.val[1], c.val[3]);

        c.val[0] = vorrq_u16(c.val[0], c.val[1]);

        e = vorrq_u16(e, c.val[0]);

        if (vmaxvq_u16(e))
        {
            return 1;
        }
    }
    return 0;
}
