/*
 * Speed benchmark code for Falcon implementation.
 *
 * ==========================(LICENSE BEGIN)============================
 *
 * Copyright (c) 2017-2019  Falcon Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * ===========================(LICENSE END)=============================
 *
 * @author   Thomas Pornin <thomas.pornin@nccgroup.com>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <papi.h>
#include "config.h"

#define TIME(s) s = PAPI_get_virt_cyc();
#define CALC(start, end, ntests) (end - start)/ntests;

int handle_error(int retval)
{
    printf("Error: ");
    PAPI_perror(PAPI_strerror(retval));
    exit(1);
}

/*
 * This code uses only the external API.
 */

#include "falcon.h"

static void *
xmalloc(size_t len)
{
    void *buf;

    if (len == 0)
    {
        return NULL;
    }
    buf = malloc(len);
    if (buf == NULL)
    {
        fprintf(stderr, "memory allocation error\n");
        exit(EXIT_FAILURE);
    }
    return buf;
}

static void
xfree(void *buf)
{
    if (buf != NULL)
    {
        free(buf);
    }
}

/*
 * Benchmark function takes an opaque context and an iteration count;
 * it returns 0 on success, a negative error code on error.
 */
typedef int (*bench_fun)(void *ctx, unsigned long num);

static long long
do_bench_cycles(bench_fun bf, void *ctx, int iteration)
{
    long long start, end;

    /*
     * Always do a few blank runs to "train" the caches and branch
     * prediction.
     */
    bf(ctx, 10);

    TIME(start);

    // Benchmark cycles
    bf(ctx, iteration);

    TIME(end);

    return CALC(start, end, iteration);
}

/*
 * Returned value is the time per iteration in nanoseconds. If the
 * benchmark function reports an error, 0.0 is returned.
 */
static double
do_bench(bench_fun bf, void *ctx, double threshold)
{
    unsigned long num;
    int r;

    /*
     * Always do a few blank runs to "train" the caches and branch
     * prediction.
     */
    r = bf(ctx, 5);
    if (r != 0)
    {
        fprintf(stderr, "ERR: %d\n", r);
        return 0.0;
    }

    num = 1;
    for (;;)
    {
        clock_t begin, end;
        double tt;

        begin = clock();
        r = bf(ctx, num);
        end = clock();
        if (r != 0)
        {
            fprintf(stderr, "ERR: %d\n", r);
            return 0.0;
        }
        tt = (double)(end - begin) / (double)CLOCKS_PER_SEC;
        if (tt >= threshold)
        {
            return tt * 1000000000.0 / (double)num;
        }

        /*
         * If the function ran for less than 0.1 seconds then
         * we simply double the iteration number; otherwise, we
         * use the run time to try to get a "correct" number of
         * iterations quickly.
         */
        if (tt < 0.1)
        {
            num <<= 1;
        }
        else
        {
            unsigned long num2;

            num2 = (unsigned long)((double)num * (threshold * 1.1) / tt);
            if (num2 <= num)
            {
                num2 = num + 1;
            }
            num = num2;
        }
    }
}

typedef struct
{
    unsigned logn;
    shake256_context rng;
    uint8_t *tmp;
    size_t tmp_len;
    uint8_t *pk;
    uint8_t *sk;
    uint8_t *esk;
    uint8_t *sig;
    size_t sig_len;
    uint8_t *sigct;
    size_t sigct_len;
} bench_context;

static inline size_t
maxsz(size_t a, size_t b)
{
    return a > b ? a : b;
}

#define CC(x)           \
    do                  \
    {                   \
        int ccr = (x);  \
        if (ccr != 0)   \
        {               \
            return ccr; \
        }               \
    } while (0)

static int
bench_keygen(void *ctx, unsigned long num)
{
    bench_context *bc;

    bc = ctx;
    while (num-- > 0)
    {
        CC(falcon_keygen_make(&bc->rng, bc->logn,
                              bc->sk, FALCON_PRIVKEY_SIZE(bc->logn),
                              bc->pk, FALCON_PUBKEY_SIZE(bc->logn),
                              bc->tmp, bc->tmp_len));
    }
    return 0;
}

static int
bench_sign_dyn(void *ctx, unsigned long num)
{
    bench_context *bc;

    bc = ctx;
    while (num-- > 0)
    {
        bc->sig_len = FALCON_SIG_COMPRESSED_MAXSIZE(bc->logn);
        CC(falcon_sign_dyn(&bc->rng,
                           bc->sig, &bc->sig_len, FALCON_SIG_COMPRESSED,
                           bc->sk, FALCON_PRIVKEY_SIZE(bc->logn),
                           "data", 4, bc->tmp, bc->tmp_len));
    }
    return 0;
}

static int
bench_sign_dyn_ct(void *ctx, unsigned long num)
{
    bench_context *bc;

    bc = ctx;
    while (num-- > 0)
    {
        bc->sigct_len = FALCON_SIG_CT_SIZE(bc->logn);
        CC(falcon_sign_dyn(&bc->rng,
                           bc->sigct, &bc->sigct_len, FALCON_SIG_CT,
                           bc->sk, FALCON_PRIVKEY_SIZE(bc->logn),
                           "data", 4, bc->tmp, bc->tmp_len));
    }
    return 0;
}

static int
bench_expand_privkey(void *ctx, unsigned long num)
{
    bench_context *bc;

    bc = ctx;
    while (num-- > 0)
    {
        CC(falcon_expand_privkey(
            bc->esk, FALCON_EXPANDEDKEY_SIZE(bc->logn),
            bc->sk, FALCON_PRIVKEY_SIZE(bc->logn),
            bc->tmp, bc->tmp_len));
    }
    return 0;
}

static int
bench_sign_tree(void *ctx, unsigned long num)
{
    bench_context *bc;

    bc = ctx;
    while (num-- > 0)
    {
        bc->sig_len = FALCON_SIG_COMPRESSED_MAXSIZE(bc->logn);
        CC(falcon_sign_tree(&bc->rng,
                            bc->sig, &bc->sig_len, FALCON_SIG_COMPRESSED,
                            bc->esk,
                            "data", 4, bc->tmp, bc->tmp_len));
    }
    return 0;
}

static int
bench_sign_tree_ct(void *ctx, unsigned long num)
{
    bench_context *bc;

    bc = ctx;
    while (num-- > 0)
    {
        bc->sigct_len = FALCON_SIG_CT_SIZE(bc->logn);
        CC(falcon_sign_tree(&bc->rng,
                            bc->sigct, &bc->sigct_len, FALCON_SIG_CT,
                            bc->esk,
                            "data", 4, bc->tmp, bc->tmp_len));
    }
    return 0;
}

static int
bench_verify(void *ctx, unsigned long num)
{
    bench_context *bc;
    size_t pk_len;

    bc = ctx;
    pk_len = FALCON_PUBKEY_SIZE(bc->logn);
    while (num-- > 0)
    {
        CC(falcon_verify(
            bc->sig, bc->sig_len, FALCON_SIG_COMPRESSED,
            bc->pk, pk_len,
            "data", 4, bc->tmp, bc->tmp_len));
    }
    return 0;
}

static int
bench_verify_ct(void *ctx, unsigned long num)
{
    bench_context *bc;
    size_t pk_len;

    bc = ctx;
    pk_len = FALCON_PUBKEY_SIZE(bc->logn);
    while (num-- > 0)
    {
        CC(falcon_verify(
            bc->sigct, bc->sigct_len, FALCON_SIG_CT,
            bc->pk, pk_len,
            "data", 4, bc->tmp, bc->tmp_len));
    }
    return 0;
}

static void
test_speed_falcon_cycles(unsigned logn, int iteration)
{
    printf("All numbers are in cycles\n\n");
    printf("degree  kg(c) \t  ek(c) \t  sd(c) \t sdc(c) \t  vv(c) \t vvc(c)\n");

    bench_context bc;

    size_t len;
    printf("%4u:", 1u << logn);
    fflush(stdout);

    bc.logn = logn;
    if (shake256_init_prng_from_system(&bc.rng) != 0)
    {
        fprintf(stderr, "random seeding failed\n");
        exit(EXIT_FAILURE);
    }
    len = FALCON_TMPSIZE_KEYGEN(logn);
    len = maxsz(len, FALCON_TMPSIZE_SIGNDYN(logn));
    len = maxsz(len, FALCON_TMPSIZE_SIGNTREE(logn));
    len = maxsz(len, FALCON_TMPSIZE_EXPANDPRIV(logn));
    len = maxsz(len, FALCON_TMPSIZE_VERIFY(logn));
    bc.tmp = xmalloc(len);
    bc.tmp_len = len;
    bc.pk = xmalloc(FALCON_PUBKEY_SIZE(logn));
    bc.sk = xmalloc(FALCON_PRIVKEY_SIZE(logn));
    bc.esk = xmalloc(FALCON_EXPANDEDKEY_SIZE(logn));
    bc.sig = xmalloc(FALCON_SIG_COMPRESSED_MAXSIZE(logn));
    bc.sig_len = 0;
    bc.sigct = xmalloc(FALCON_SIG_CT_SIZE(logn));
    bc.sigct_len = 0;

    printf(" %8lld",
           do_bench_cycles(&bench_keygen, &bc, iteration / 10));
    fflush(stdout);
    printf(" %8lld",
           do_bench_cycles(&bench_expand_privkey, &bc, iteration));
    fflush(stdout);
    printf(" %8lld",
           do_bench_cycles(&bench_sign_dyn, &bc, iteration));
    fflush(stdout);
    printf(" %8lld",
           do_bench_cycles(&bench_sign_dyn_ct, &bc, iteration));
    fflush(stdout);
    printf(" %8lld",
            do_bench_cycles(&bench_sign_tree, &bc, iteration));
    fflush(stdout);
    printf(" %8lld",
            do_bench_cycles(&bench_sign_tree_ct, &bc, iteration));
    fflush(stdout);
    printf(" %8lld",
           do_bench_cycles(&bench_verify, &bc, iteration));
    fflush(stdout);
    printf(" %8lld",
           do_bench_cycles(&bench_verify_ct, &bc, iteration));
    fflush(stdout);

    printf("\n\n");
    fflush(stdout);

    xfree(bc.tmp);
    xfree(bc.pk);
    xfree(bc.sk);
    xfree(bc.esk);
    xfree(bc.sig);
    xfree(bc.sigct);
}

static void
test_speed_falcon(unsigned logn, double threshold)
{
    printf("keygen in milliseconds, other values in microseconds\n");
    printf("\n");
    // printf("degree  kg(ms)   ek(us)   sd(us)  sdc(us)   vv(us)  vvc(us)\n");
    printf("degree  kg(ms)   ek(us)   sd(us)  sdc(us)   st(us)  stc(us)   vv(us)  vvc(us)\n");

    bench_context bc;
    size_t len;

    printf("%4u:", 1u << logn);
    fflush(stdout);

    bc.logn = logn;
    if (shake256_init_prng_from_system(&bc.rng) != 0)
    {
        fprintf(stderr, "random seeding failed\n");
        exit(EXIT_FAILURE);
    }
    len = FALCON_TMPSIZE_KEYGEN(logn);
    len = maxsz(len, FALCON_TMPSIZE_SIGNDYN(logn));
    len = maxsz(len, FALCON_TMPSIZE_SIGNTREE(logn));
    len = maxsz(len, FALCON_TMPSIZE_EXPANDPRIV(logn));
    len = maxsz(len, FALCON_TMPSIZE_VERIFY(logn));
    bc.tmp = xmalloc(len);
    bc.tmp_len = len;
    bc.pk = xmalloc(FALCON_PUBKEY_SIZE(logn));
    bc.sk = xmalloc(FALCON_PRIVKEY_SIZE(logn));
    bc.esk = xmalloc(FALCON_EXPANDEDKEY_SIZE(logn));
    bc.sig = xmalloc(FALCON_SIG_COMPRESSED_MAXSIZE(logn));
    bc.sig_len = 0;
    bc.sigct = xmalloc(FALCON_SIG_CT_SIZE(logn));
    bc.sigct_len = 0;

    printf(" %8.2f",
           do_bench(&bench_keygen, &bc, threshold) / 1000000.0);
    fflush(stdout);
    printf(" %8.2f",
           do_bench(&bench_expand_privkey, &bc, threshold) / 1000.0);
    fflush(stdout);
    printf(" %8.2f",
           do_bench(&bench_sign_dyn, &bc, threshold) / 1000.0);
    fflush(stdout);
    printf(" %8.2f",
           do_bench(&bench_sign_dyn_ct, &bc, threshold) / 1000.0);
    fflush(stdout);
    printf(" %8.2f",
            do_bench(&bench_sign_tree, &bc, threshold) / 1000.0);
    fflush(stdout);
    printf(" %8.2f",
            do_bench(&bench_sign_tree_ct, &bc, threshold) / 1000.0);
    fflush(stdout);
    printf(" %8.2f",
           do_bench(&bench_verify, &bc, threshold) / 1000.0);
    fflush(stdout);
    printf(" %8.2f",
           do_bench(&bench_verify_ct, &bc, threshold) / 1000.0);
    fflush(stdout);

    printf("\n\n");
    fflush(stdout);

    xfree(bc.tmp);
    xfree(bc.pk);
    xfree(bc.sk);
    xfree(bc.esk);
    xfree(bc.sig);
    xfree(bc.sigct);
}

int main()
{

    // PAPI stuffs https://bitbucket.org/icl/papi/wiki/PAPI-LL
    /* Initialize the PAPI library */
    int retval;

    retval = PAPI_library_init(PAPI_VER_CURRENT);

    if (retval != PAPI_VER_CURRENT && retval > 0)
    {
        fprintf(stderr, "PAPI library version mismatch!\n");
        exit(1);
    }

    if (retval < 0)
    {
        fprintf(stderr, "Initialization error!\n");
        exit(1);
    }

    //
    double threshold;
    int iteration;

    threshold = 2.0;
    iteration = 10000;

    printf("time threshold = %.4f s\n", threshold);
    printf("kg = keygen, ek = expand private key, sd = sign (without expanded key)\n");
    printf("st = sign (with expanded key), vv = verify\n");
    printf("sdc, stc, vvc: like sd, st and vv, but with constant-time hash-to-point\n");

    test_speed_falcon(FALCON_LOGN, threshold);
    test_speed_falcon_cycles(FALCON_LOGN, iteration);
    return 0;
}