#include "inner.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stddef.h>
#include "util.h"
#include "config.h"
#include "api.h"
#include "poly.h"

#if _APPLE_M1_ == 1
#include "m1cycles.h"

// Result is cycle per call
#define TIME(s) s = rdtsc();
#define CALC(start, stop, ntests) (stop - start) / ntests;
#else
// Result is nanosecond per call
#define TIME(s) clock_gettime(CLOCK_MONOTONIC_RAW, &s);
#define CALC(start, stop, ntests) ((double)((stop.tv_sec - start.tv_sec) * 1000000000 + (stop.tv_nsec - start.tv_nsec))) / ntests;
#endif

#define NTESTS 10000
#define FALCON_LOGN 10
#define FALCON_N (1 << FALCON_LOGN)
#define FALCON_Q 12289

void test_FFT(fpr *f, unsigned logn)
{
#if _APPLE_M1_ == 0
    struct timespec start, stop;
#else
    long long start_cc, stop_cc;
#endif
    long long fft, ifft;
    unsigned ntests = NTESTS;
    if (logn < 7)
    {
        ntests = NTESTS * 100;
    }
    /* =================================== */
    TIME(start);
    for (unsigned i = 0; i < ntests; i++)
    {
        ZfN(FFT)(f, logn);
    }
    TIME(stop);
    fft = CALC(start, stop, ntests);

    /* =================================== */
    TIME(start);
    for (unsigned i = 0; i < ntests; i++)
    {
        ZfN(iFFT)(f, logn);
    }
    TIME(stop);

    /* =================================== */
    ifft = CALC(start, stop, ntests);
    printf("FFT %u: %lld - %lld\n", logn, fft, ifft);
    printf("=======\n");
}

void test_NTT(uint16_t *a, unsigned logn)
{
#if _APPLE_M1_ == 0
    struct timespec start, stop;
#else
    long long start_cc, stop_cc;
#endif
    long long fft, ifft;
    unsigned ntests = NTESTS;
    if (logn < 7)
    {
        ntests = NTESTS * 100;
    }
    /* =================================== */
    TIME(start);
    for (unsigned i = 0; i < ntests; i++)
    {
        ZfN(poly_ntt)(a, 0);
    }
    TIME(stop);
    fft = CALC(start, stop, ntests);

    /* =================================== */
    TIME(start);
    for (unsigned i = 0; i < ntests; i++)
    {
        ZfN(poly_invntt)(a);
    }
    TIME(stop);

    /* =================================== */
    ifft = CALC(start, stop, ntests);
    printf("NTT %u: %lld - %lld\n", logn, fft, ifft);
    printf("=======\n");
}

void test_poly_add(fpr *c, fpr *a, fpr *b, unsigned logn, char *string)
{
#if _APPLE_M1_ == 0
    struct timespec start, stop;
#else
    long long start_cc, stop_cc;
#endif
    long long fft, ifft;
    unsigned ntests = NTESTS;
    if (logn < 10)
    {
        ntests = NTESTS * 1000;
    }
    /* =================================== */
    TIME(start);
    for (unsigned i = 0; i < ntests; i++)
    {
        ZfN(poly_add)(c, a, b, logn);
    }
    TIME(stop);
    fft = CALC(start, stop, ntests);

    printf("%s, %u, %lld\n", string, logn, fft);
}

void test_poly_sub(fpr *c, fpr *a, fpr *b, unsigned logn, char *string)
{
#if _APPLE_M1_ == 0
    struct timespec start, stop;
#else
    long long start_cc, stop_cc;
#endif
    long long fft, ifft;
    unsigned ntests = NTESTS;
    if (logn < 10)
    {
        ntests = NTESTS * 1000;
    }
    /* =================================== */
    TIME(start);
    for (unsigned i = 0; i < ntests; i++)
    {
        ZfN(poly_sub)(c, a, b, logn);
    }
    TIME(stop);
    fft = CALC(start, stop, ntests);

    printf("%s, %u, %lld\n", string, logn, fft);
}

void test_poly_neg(fpr *c, fpr *a, fpr *b, unsigned logn, char *string)
{
#if _APPLE_M1_ == 0
    struct timespec start, stop;
#else
    long long start_cc, stop_cc;
#endif
    long long fft, ifft;
    unsigned ntests = NTESTS;
    if (logn < 10)
    {
        ntests = NTESTS * 1000;
    }
    /* =================================== */
    TIME(start);
    for (unsigned i = 0; i < ntests; i++)
    {
        ZfN(poly_neg)(c, a, logn);
    }
    TIME(stop);
    fft = CALC(start, stop, ntests);

    printf("%s, %u, %lld\n", string, logn, fft);
}

void test_poly_adj_fft(fpr *c, fpr *a, fpr *b, unsigned logn, char *string)
{
#if _APPLE_M1_ == 0
    struct timespec start, stop;
#else
    long long start_cc, stop_cc;
#endif
    long long fft, ifft;
    unsigned ntests = NTESTS;
    if (logn < 10)
    {
        ntests = NTESTS * 1000;
    }
    /* =================================== */
    TIME(start);
    for (unsigned i = 0; i < ntests; i++)
    {
        ZfN(poly_adj_fft)(c, a, logn);
    }
    TIME(stop);
    fft = CALC(start, stop, ntests);

    printf("%s, %u, %lld\n", string, logn, fft);
}

void test_poly_mul_fft(fpr *c, fpr *a, fpr *b, unsigned logn, char *string)
{
#if _APPLE_M1_ == 0
    struct timespec start, stop;
#else
    long long start_cc, stop_cc;
#endif
    long long fft, ifft;
    unsigned ntests = NTESTS;
    if (logn < 7)
    {
        ntests = NTESTS * 100;
    }
    /* =================================== */
    TIME(start);
    for (unsigned i = 0; i < ntests; i++)
    {
        ZfN(poly_mul_fft)(c, a, b, logn);
    }
    TIME(stop);
    fft = CALC(start, stop, ntests);

    printf("%s, %u, %lld\n", string, logn, fft);
}

void test_poly_invnorm2_fft(fpr *c, fpr *a, fpr *b, unsigned logn, char *string)
{
#if _APPLE_M1_ == 0
    struct timespec start, stop;
#else
    long long start_cc, stop_cc;
#endif
    long long fft, ifft;
    unsigned ntests = NTESTS;
    if (logn < 7)
    {
        ntests = NTESTS * 100;
    }
    /* =================================== */
    TIME(start);
    for (unsigned i = 0; i < ntests; i++)
    {
        ZfN(poly_invnorm2_fft)(c, a, b, logn);
    }
    TIME(stop);
    fft = CALC(start, stop, ntests);

    printf("%s, %u, %lld\n", string, logn, fft);
}

void test_poly_mul_autoadj_fft(fpr *c, fpr *a, fpr *b, unsigned logn, char *string)
{
#if _APPLE_M1_ == 0
    struct timespec start, stop;
#else
    long long start_cc, stop_cc;
#endif
    long long fft, ifft;
    unsigned ntests = NTESTS;
    if (logn < 10)
    {
        ntests = NTESTS * 1000;
    }
    /* =================================== */
    TIME(start);
    for (unsigned i = 0; i < ntests; i++)
    {
        ZfN(poly_mul_autoadj_fft)(c, a, b, logn);
    }
    TIME(stop);
    fft = CALC(start, stop, ntests);

    printf("%s, %u, %lld\n", string, logn, fft);
}

void test_poly_LDL_fft(fpr *c, fpr *a, fpr *b, unsigned logn, char *string)
{
#if _APPLE_M1_ == 0
    struct timespec start, stop;
#else
    long long start_cc, stop_cc;
#endif
    long long fft, ifft;
    unsigned ntests = NTESTS;
    if (logn < 10)
    {
        ntests = NTESTS * 1000;
    }
    /* =================================== */
    TIME(start);
    for (unsigned i = 0; i < ntests; i++)
    {
        ZfN(poly_LDL_fft)(c, a, b, logn);
    }
    TIME(stop);
    fft = CALC(start, stop, ntests);

    printf("%s, %u, %lld\n", string, logn, fft);
}

void test_poly_LDLmv_fft(fpr *d11, fpr *l01, const fpr *c,const  fpr *a, const  fpr *b, unsigned logn, char *string)
{
#if _APPLE_M1_ == 0
    struct timespec start, stop;
#else
    long long start_cc, stop_cc;
#endif
    long long fft, ifft;
    unsigned ntests = NTESTS;
    if (logn < 10)
    {
        ntests = NTESTS * 1000;
    }
    /* =================================== */
    TIME(start);
    for (unsigned i = 0; i < ntests; i++)
    {
        ZfN(poly_LDLmv_fft)(d11, l01, c, a, b, logn);
    }
    TIME(stop);
    fft = CALC(start, stop, ntests);

    printf("%s, %u, %lld\n", string, logn, fft);
    printf("=======\n");
}

int main()
{
    fpr f[FALCON_N], fa[FALCON_N], fb[FALCON_N], fc[FALCON_N], tmp[FALCON_N]={0};
    uint16_t a[FALCON_N];
    for (int i = 0; i < FALCON_N; i++)
    {
        double_t t;
        t = (double)i;
        f[i] = t;
        fa[i] = t;
        fb[i] = t;
        fc[i] = t;
        a[i] = rand() % FALCON_Q;
    }

#if _APPLE_M1_ == 1
    setup_rdtsc();
#endif

    for (unsigned i = 0; i <= FALCON_LOGN; i++)
    {
        test_FFT(f, i);
    }

    test_NTT(a, FALCON_LOGN);

    for (unsigned i = 0; i <= FALCON_LOGN; i++)
    {
        test_poly_add(fc, fa, fb, i, "poly_add");
        test_poly_sub(fc, fa, fb, i, "poly_sub");
        test_poly_neg(fc, fa, fb, i, "poly_neg");
        test_poly_adj_fft(fc, fa, fb, i, "poly_adj_fft");
        test_poly_mul_fft(fc, fa, fb, i, "poly_mul_fft");
        test_poly_invnorm2_fft(fc, fa, fb, i, "poly_invnorm2_fft");
        test_poly_mul_autoadj_fft(fc, fa, fb, i, "poly_mul_autoadj_fft");
        test_poly_LDL_fft(fc, fa, fb, i, "poly_LDL_fft");
        test_poly_LDLmv_fft(f, tmp, fc, fa, fb, i, "poly_LDLmv_fft");
    }

    return 0;
}

/*
 * Result in nanosection
 * FFT          - iFFT
 *  9: 4609     - 4333
 * 10: 11389    - 10790
 *
 * NTT          - iNTT
 *  9: 2560     - 2646
 * 10: 5426     - 5721
 */
