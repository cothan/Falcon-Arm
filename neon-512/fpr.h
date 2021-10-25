#ifndef PQCLEAN_FALCON512_NEON_FPR_H
#define PQCLEAN_FALCON512_NEON_FPR_H

#include <math.h>

/*
 * We wrap the native 'double' type into a structure so that the C compiler
 * complains if we inadvertently use raw arithmetic operators on the 'fpr'
 * type instead of using the inline functions below. This should have no
 * extra runtime cost, since all the functions below are 'inline'.
 */
typedef double fpr;

static inline fpr
fpr_of(int64_t i)
{
	return (double)i;
}

static const fpr fpr_q = 12289.0 ;
static const fpr fpr_inverse_of_q = 1.0 / 12289.0 ;
static const fpr fpr_inv_2sqrsigma0 = .150865048875372721532312163019 ;
static const fpr fpr_inv_sigma = .005819826392951607426919370871 ;
static const fpr fpr_sigma_min_9 = 1.291500756233514568549480827642 ;
static const fpr fpr_sigma_min_10 = 1.311734375905083682667395805765 ;
static const fpr fpr_log2 = 0.69314718055994530941723212146 ;
static const fpr fpr_inv_log2 = 1.4426950408889634073599246810 ;
static const fpr fpr_bnorm_max = 16822.4121 ;
static const fpr fpr_zero = 0.0 ;
static const fpr fpr_one = 1.0 ;
static const fpr fpr_two = 2.0 ;
static const fpr fpr_onehalf = 0.5 ;
static const fpr fpr_invsqrt2 = 0.707106781186547524400844362105 ;
static const fpr fpr_invsqrt8 = 0.353553390593273762200422181052 ;
static const fpr fpr_ptwo31 = 2147483648.0 ;
static const fpr fpr_ptwo31m1 = 2147483647.0 ;
static const fpr fpr_mtwo31m1 = -2147483647.0 ;
static const fpr fpr_ptwo63m1 = 9223372036854775807.0 ;
static const fpr fpr_mtwo63m1 = -9223372036854775807.0 ;
static const fpr fpr_ptwo63 = 9223372036854775808.0 ;

static inline int64_t
fpr_rint(fpr x)
{
	/*
	 * We do not want to use llrint() since it might be not
	 * constant-time.
	 *
	 * Suppose that x >= 0. If x >= 2^52, then it is already an
	 * integer. Otherwise, if x < 2^52, then computing x+2^52 will
	 * yield a value that will be rounded to the nearest integer
	 * with exactly the right rules (round-to-nearest-even).
	 *
	 * In order to have constant-time processing, we must do the
	 * computation for both x >= 0 and x < 0 cases, and use a
	 * cast to an integer to access the sign and select the proper
	 * value. Such casts also allow us to find out if |x| < 2^52.
	 */
	int64_t sx, tx, rp, rn, m;
	uint32_t ub;

	sx = (int64_t)(x - 1.0);
	tx = (int64_t)x;
	rp = (int64_t)(x + 4503599627370496.0) - 4503599627370496;
	rn = (int64_t)(x - 4503599627370496.0) + 4503599627370496;

	/*
	 * If tx >= 2^52 or tx < -2^52, then result is tx.
	 * Otherwise, if sx >= 0, then result is rp.
	 * Otherwise, result is rn. We use the fact that when x is
	 * close to 0 (|x| <= 0.25) then both rp and rn are correct;
	 * and if x is not close to 0, then trunc(x-1.0) yields the
	 * appropriate sign.
	 */

	/*
	 * Clamp rp to zero if tx < 0.
	 * Clamp rn to zero if tx >= 0.
	 */
	m = sx >> 63;
	rn &= m;
	rp &= ~m;

	/*
	 * Get the 12 upper bits of tx; if they are not all zeros or
	 * all ones, then tx >= 2^52 or tx < -2^52, and we clamp both
	 * rp and rn to zero. Otherwise, we clamp tx to zero.
	 */
	ub = (uint32_t)((uint64_t)tx >> 52);
	m = -(int64_t)((((ub + 1) & 0xFFF) - 2) >> 31);
	rp &= m;
	rn &= m;
	tx &= ~m;

	/*
	 * Only one of tx, rn or rp (at most) can be non-zero at this
	 * point.
	 */
	return tx | rn | rp;
}

static inline int64_t
fpr_floor(fpr x)
{
	int64_t r;

	/*
	 * The cast performs a trunc() (rounding toward 0) and thus is
	 * wrong by 1 for most negative values. The correction below is
	 * constant-time as long as the compiler turns the
	 * floating-point conversion result into a 0/1 integer without a
	 * conditional branch or another non-constant-time construction.
	 * This should hold on all modern architectures with an FPU (and
	 * if it is false on a given arch, then chances are that the FPU
	 * itself is not constant-time, making the point moot).
	 */
	r = (int64_t)x;
	return r - (x < (double)r);
}

static inline int64_t
fpr_trunc(fpr x)
{
	return (int64_t)x;
}

static inline fpr
fpr_add(fpr x, fpr y)
{
	return (x + y);
}

static inline fpr
fpr_sub(fpr x, fpr y)
{
	return (x - y);
}

static inline fpr
fpr_neg(fpr x)
{
	return (-x);
}

static inline fpr
fpr_half(fpr x)
{
	return (x * 0.5);
}

static inline fpr
fpr_double(fpr x)
{
	return (x + x);
}

static inline fpr
fpr_mul(fpr x, fpr y)
{
	return (x * y);
}

static inline fpr
fpr_sqr(fpr x)
{
	return (x * x);
}

static inline fpr
fpr_inv(fpr x)
{
	return (1.0 / x);
}

static inline fpr
fpr_div(fpr x, fpr y)
{
	return (x / y);
}

static inline fpr
fpr_sqrt(fpr x)
{

#if defined __aarch64__ && __aarch64__
	__asm__("fsqrt   %d0, %d0"
			: "+w"(x)
			:
			:);
#else
	__asm__("fsqrtd  %P0, %P0"
			: "+w"(x)
			:
			:);
#endif
	return x;
}

static inline int
fpr_lt(fpr x, fpr y)
{
	return x < y;
}

static inline uint64_t
fpr_expm_p63(fpr x, fpr ccs)
{
	/*
	 * Polynomial approximation of exp(-x) is taken from FACCT:
	 *   https://eprint.iacr.org/2018/1234
	 * Specifically, values are extracted from the implementation
	 * referenced from the FACCT article, and available at:
	 *   https://github.com/raykzhao/gaussian
	 * Tests over more than 24 billions of random inputs in the
	 * 0..log(2) range have never shown a deviation larger than
	 * 2^(-50) from the true mathematical value.
	 */
	/*
	 * Normal implementation uses Horner's method, which minimizes
	 * the number of operations.
	 */

	double d, y;

	d = x;
	y = 0.000000002073772366009083061987;
	y = 0.000000025299506379442070029551 - y * d;
	y = 0.000000275607356160477811864927 - y * d;
	y = 0.000002755586350219122514855659 - y * d;
	y = 0.000024801566833585381209939524 - y * d;
	y = 0.000198412739277311890541063977 - y * d;
	y = 0.001388888894063186997887560103 - y * d;
	y = 0.008333333327800835146903501993 - y * d;
	y = 0.041666666666110491190622155955 - y * d;
	y = 0.166666666666984014666397229121 - y * d;
	y = 0.500000000000019206858326015208 - y * d;
	y = 0.999999999999994892974086724280 - y * d;
	y = 1.000000000000000000000000000000 - y * d;
	y *= ccs;
	return (uint64_t)(y * fpr_ptwo63);
}


#endif