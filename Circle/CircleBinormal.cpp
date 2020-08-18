#include "CircleBinormal.h"
#include "CircleIntersect.h"
//#include "../Binomial.h"

#define ROOT_EPS				1e-10		// [ x ] must satisfy |f(x)| <= ROOT_EPS to be root  
#define TOL						1e-5		// Epsilon for numerical stability : Search domain [ -1, 1 ] is extended by this
#define DOMAIN_EPS				1e-3		// Maximum domain width to find a single root ( If multiple roots coexist in the domain narrower than this, only one of them would be found )
#define NR_MAX_ITER				20			// Maximum step to take when doing numerical search
#define PROJECTION_EPS			1e-1		// Tolerance when projecting point from a circular arc to the other circular arc ( for stability )
#define BINORMAL_LUDCMP_EPS		1.0e-20		// EPS for LU Decomposition in binormal routine
#define PROXIMITY_EPS			1e-10		// Two points are considered to be same if distance between them is lower than this value
#define PROXIMITY_EPS2			1e-5		// Two points are close enough to make binormal objective function value unstable

namespace MN {
	static const Binomial oBin = Binomial::create(8);
	static const Real binomial[9][9] = {
		{ oBin.at(0, 0) },
		{ oBin.at(1, 0), oBin.at(1, 1) },
		{ oBin.at(2, 0), oBin.at(2, 1), oBin.at(2, 2) },
		{ oBin.at(3, 0), oBin.at(3, 1), oBin.at(3, 2), oBin.at(3, 3) },
		{ oBin.at(4, 0), oBin.at(4, 1), oBin.at(4, 2), oBin.at(4, 3), oBin.at(4, 4) },
		{ oBin.at(5, 0), oBin.at(5, 1), oBin.at(5, 2), oBin.at(5, 3), oBin.at(5, 4), oBin.at(5, 5) },
		{ oBin.at(6, 0), oBin.at(6, 1), oBin.at(6, 2), oBin.at(6, 3), oBin.at(6, 4), oBin.at(6, 5), oBin.at(6, 6) },
		{ oBin.at(7, 0), oBin.at(7, 1), oBin.at(7, 2), oBin.at(7, 3), oBin.at(7, 4), oBin.at(7, 5), oBin.at(7, 6), oBin.at(7, 7) },
		{ oBin.at(8, 0), oBin.at(8, 1), oBin.at(8, 2), oBin.at(8, 3), oBin.at(8, 4), oBin.at(8, 5), oBin.at(8, 6), oBin.at(8, 7), oBin.at(8, 8) }
	};

	// BP
	inline static void subdivideBP(const CircleBinormal::BP& M, Real t, CircleBinormal::BP& L, CircleBinormal::BP& R) {
		Real iCoefs[9];
		memcpy(iCoefs, M.coefs, sizeof(Real) * (M.degree + 1));
		L.degree = M.degree;
		R.degree = M.degree;

		int cnt = 0;
		int degree = M.degree;
		Real t_1 = 1.0 - t;
		do {
			L.coefs[cnt] = iCoefs[0];
			R.coefs[degree - cnt] = iCoefs[degree - cnt];

			for (int i = 0; i < degree - cnt; i++)
				iCoefs[i] = t_1 * iCoefs[i] + t * iCoefs[i + 1];
		} while (cnt++ < degree);
	}
	inline static bool purgeBP(const CircleBinormal::BP& M, Real validDomain[2]) {
		// Check if M's domain meets [ validDomain ]
		// If not, purge it
		if (M.domain[1] < validDomain[0] || M.domain[0] > validDomain[1])
			return true;
		// Assume none of coefficients is zero
		bool positive = M.coefs[0] > 0, tmpPos;
		for (size_t i = 0; i < M.degree + 1; i++) {
			tmpPos = M.coefs[i] > 0;
			if (tmpPos != positive)
				return false;
		}
		return true;
	}
	inline static bool purgeBP(const CircleBinormal::BP& M, const std::vector<Domain>& validDomains) {
		// Check if M's domain meets [ validDomains ]
		// If not, purge it
		bool inDomain = false;
		for (auto& domain : validDomains) {
			if (M.domain[1] >= domain.beg() && M.domain[0] <= domain.end()) {
				inDomain = true;
				break;
			}
		}
		if (!inDomain)
			return true;
		// Assume none of coefficients is zero
		bool positive = M.coefs[0] > 0, tmpPos;
		for (size_t i = 0; i < M.degree + 1; i++) {
			tmpPos = M.coefs[i] > 0;
			if (tmpPos != positive)
				return false;
		}
		return true;
	}
	inline static void factorBP(CircleBinormal::BP& M, bool atZero) {
		for (int i = 0; i < M.degree; i++) {
			M.coefs[i] = atZero ?
				M.coefs[i + 1] * (M.degree / (Real)(i + 1)) :
				M.coefs[i] * (M.degree / (Real)(M.degree - i));
		}
		M.degree--;
	}
	inline static void evaluateBP(const CircleBinormal::BP& M, Real t, Real& value, Real& deriv) {
		Real t_1 = 1.0 - t;
		Real memA[9];
		Real memB[9];
		memA[0] = 1.0;
		memB[0] = 1.0;
		for (int i = 1; i <= M.degree; i++) {
			memA[i] = memA[i - 1] * t_1;	// (1-t)^0, (1-t)^1, (1-t)^2, ... , (1-t)^degree
			memB[i] = memB[i - 1] * t;		// t^0, t^1, t^2, ... , t^degree
		}
		value = 0.0;
		deriv = 0.0;
		for (int i = 0; i <= M.degree; i++) {
			value += M.coefs[i] * binomial[M.degree][i] * memA[M.degree - i] * memB[i];
			if (i < M.degree)
				deriv += M.dCoefs[i] * binomial[M.degree - 1][i] * memA[M.degree - 1 - i] * memB[i];
		}
		deriv *= M.degree;
	}
	inline static bool solveNR(const CircleBinormal::BP& M, Real& t) {
		Real value, deriv, dt;
		for (int i = 0; i < NR_MAX_ITER; i++) {
			evaluateBP(M, t, value, deriv);
			if (fabs(value) < ROOT_EPS) return true;		// Found root
			if (deriv == 0.0) return false;		// Divergent step
			dt = value / deriv;
			t -= dt;
			if (t < 0 || t > 1) return false;	// Out of domain
		}
		return false;							// Max iteration
	}
	inline static Real estimateBP(const CircleBinormal::BP& M) {
		Real length = (1.0 / M.degree);
		Real t = 0;
		for (int i = 0; i < M.degree; i++) {
			if (M.coefs[i] * M.coefs[i + 1] < 0) {
				Real add = -M.coefs[i] * (length / (M.coefs[i + 1] - M.coefs[i]));
				return t + add;
			}
			t += length;
		}
		return 0.5;
	}
	CircleBinormal::BP CircleBinormal::initBP(const Real monoCoefs[9], const std::vector<Domain>& domains) {
		BP bp;
		// 1. Reparametrize monomial to [0, 1] 
		// Since cosine value spans [-1, 1], we have to shrink it into [0, 1]
		const static Real P0 = -1 - TOL;			// For numerical stability, we extend original span [-1, 1] with [ TOL ] value in both sides
		const static Real Q0 = 2 + TOL * 2.0;		// Therefore, it actually maps [-1-TOL, 1+TOL] in original domain to [0, 1]
		const static Real Ps[] = { 1, pow(P0, 1), pow(P0, 2), pow(P0, 3), pow(P0, 4), pow(P0, 5), pow(P0, 6), pow(P0, 7), pow(P0, 8) };
		const static Real Qs[] = { 1, pow(Q0, 1), pow(Q0, 2), pow(Q0, 3), pow(Q0, 4), pow(Q0, 5), pow(Q0, 6), pow(Q0, 7), pow(Q0, 8) };
		Real iCoefs[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		Real c;

		for (int i = 0; i <= 8; i++) {
			c = monoCoefs[i]; // c * t^i ----> c * (-1 + 2t)^i
			for (int j = 0; j <= i; j++)
				iCoefs[j] += c * binomial[i][j] * Ps[i - j] * Qs[j];
		}

		// 2. Change monomial to intermediate form
		for (int i = 0; i <= 8; i++)
			iCoefs[i] /= binomial[8][i];

		// 3. Change intermediate form to Bezier polynomial
		for (int j = 0; j <= 8; j++) {
			bp.coefs[j] = 0.0;
			for (int i = 0; i <= j; i++)
				bp.coefs[j] += binomial[j][i] * iCoefs[i];
		}

		// 4. Set valid domain
		// If we want to find root only in [ a, b ], in the following process, we can only search for [ 0.5a + 0.5, 0.5b + 0.5 ]
		const static Real Q1 = 1.0 / Q0;	// 0.5	
		const static Real P1 = P0 * Q1;		// -0.5
		validDomains = domains;
		for (auto& domain : validDomains) {
			Real nbeg = Q1 * (domain.beg() - TOL - P0);
			Real nend = Q1 * (domain.end() + TOL - P0);
			domain.set(nbeg, nend);
		}

		bp.degree = 8;
		bp.domain[0] = 0;
		bp.domain[1] = 1;
		return bp;
	}
	void CircleBinormal::solveBP(const Real monoCoefs[9], const std::vector<Domain>& domains, Real roots[], int& rootNum) {
		// Allocate storage for [ data ]
		static const int dataSize = (int)(1.0 / DOMAIN_EPS) + 1;
		if (data.size() < dataSize)
			data.resize(dataSize);

		// Before iteration, factor as much as we can
		BP m = initBP(monoCoefs, domains);
		rootNum = 0;
		if (m.degree > 0 && fabs(m.coefs[0]) < ROOT_EPS) {
			roots[rootNum++] = 0.0;
			while (m.degree > 0 && fabs(m.coefs[0]) < ROOT_EPS)
				factorBP(m, true);
		}
		if (m.degree > 0 && fabs(m.coefs[m.degree]) < ROOT_EPS) {
			roots[rootNum++] = 1.0;
			while (m.degree > 0 && fabs(m.coefs[m.degree]) < ROOT_EPS)
				factorBP(m, false);
		}

		if (m.degree == 0)
			return;

		// Iterative search
		int idx = 0;
		data[idx] = m;
		while (idx >= 0) {
			if (purgeBP(data[idx], validDomains)) {
				idx--;
				continue;
			}
			Real nrRoot = estimateBP(data[idx]);
			Real nrRootCopy = nrRoot;
			Real domWidth = data[idx].domain[1] - data[idx].domain[0];

			data[idx].setDcoefs();	// Have to get ready derivative coefficients before NR
			if (solveNR(data[idx], nrRoot)) {
				subdivideBP(data[idx], nrRoot, data[idx + 1], data[idx + 2]);
				do {
					factorBP(data[idx + 1], false);
					factorBP(data[idx + 2], true);
				} while (fabs(data[idx + 2].coefs[0]) < ROOT_EPS && data[idx + 2].degree > 0);

				Real root = data[idx].domain[0] + domWidth * nrRoot;
				roots[rootNum++] = root;
				if (data[idx + 1].degree > 0) {
					data[idx + 1].domain[0] = data[idx].domain[0];
					data[idx + 1].domain[1] = root;
					data[idx + 2].domain[0] = root;
					data[idx + 2].domain[1] = data[idx].domain[1];

					data[idx] = data[idx + 1];
					data[idx + 1] = data[idx + 2];
					idx++;
				}
				else
					idx--;
				continue;
			}
			else {
				if (domWidth < DOMAIN_EPS) {
					roots[rootNum++] = data[idx].domain[0] + domWidth * 0.5; // Since we do NR later, just push it
					idx--;
					continue;
				}
				if (nrRootCopy > 1 - DOMAIN_EPS)
					nrRootCopy = 1 - DOMAIN_EPS;
				else if (nrRootCopy < DOMAIN_EPS)
					nrRootCopy = DOMAIN_EPS;

				subdivideBP(data[idx], nrRootCopy, data[idx + 1], data[idx + 2]);
				Real domMid = data[idx].domain[0] + domWidth * nrRootCopy;

				data[idx + 1].domain[0] = data[idx].domain[0];
				data[idx + 1].domain[1] = domMid;
				data[idx + 2].domain[0] = domMid;
				data[idx + 2].domain[1] = data[idx].domain[1];

				data[idx] = data[idx + 1];
				data[idx + 1] = data[idx + 2];
				idx++;
				continue;
			}
		}

		// Change roots to original domain
		for (int i = 0; i < rootNum; i++)
			roots[i] = roots[i] * 2.0 * (1 + TOL) - (1 + TOL);
	}

	// NR
	inline static void binormalLUdcmp(Real a[3][3], int* idx, Real* d) {
		int i, imax, j, k;
		Real big, dum, sum, temp;
		Real vv[3];

		*d = 1.0;
		for (i = 1; i <= 2; i++) {
			big = 0.0;
			for (j = 1; j <= 2; j++)
				if ((temp = fabs(a[i][j])) > big) big = temp;
			if (big == 0.0)
				throw(std::runtime_error("Singular matrix in routine ludcmp"));
			vv[i] = 1.0 / big;
		}
		for (j = 1; j <= 2; j++) {
			for (i = 1; i < j; i++) {
				sum = a[i][j];
				for (k = 1; k < i; k++) sum -= a[i][k] * a[k][j];
				a[i][j] = sum;
			}
			big = 0.0;
			for (i = j; i <= 2; i++) {
				sum = a[i][j];
				for (k = 1; k < j; k++)
					sum -= a[i][k] * a[k][j];
				a[i][j] = sum;
				if ((dum = vv[i] * fabs(sum)) >= big) {
					big = dum;
					imax = i;
				}
			}
			if (j != imax) {
				for (k = 1; k <= 2; k++) {
					dum = a[imax][k];
					a[imax][k] = a[j][k];
					a[j][k] = dum;
				}
				*d = -(*d);
				vv[imax] = vv[j];
			}
			idx[j] = imax;
			if (a[j][j] == 0.0)
				a[j][j] = BINORMAL_LUDCMP_EPS;
			if (j != 2) {
				dum = 1.0 / (a[j][j]);
				for (i = j + 1; i <= 2; i++) a[i][j] *= dum;
			}
		}
	}
	inline static void binormalLUbksb(Real a[3][3], int* idx, Real b[]) {
		int i, ii = 0, ip, j;
		Real sum;

		for (i = 1; i <= 2; i++) {
			ip = idx[i];
			sum = b[ip];
			b[ip] = b[i];
			if (ii)
				for (j = ii; j <= i - 1; j++) sum -= a[i][j] * b[j];
			else if (sum) ii = i;
			b[i] = sum;
		}
		for (i = 2; i >= 1; i--) {
			sum = b[i];
			for (j = i + 1; j <= 2; j++) sum -= a[i][j] * b[j];
			b[i] = sum / a[i][i];
		}
	}
	// For given parameters, compute objective function
	// @dvec : Objective function value
	// @djac : Jacobian matrix of objective function
	// @ret :	If two points are close enough, |f(x) - g(y)| converges to zero and objective function becomes unstable. 
	//			Therefore, in those cases, we return TRUE to indicate that those points are close
	/*inline static bool binormalNRFunc(Real param[3], Real dvec[3], Real djac[3][3], const Circle& a, const Circle& b, const Transform& btoa) {
		Real
			rx = a.radius,
			ry = b.radius,
			x_ = param[1],				// Parameter A
			y_ = param[2];				// Parameter B
		Vec3
			afx = a.evaluate(x_),		// Point A in A's local coordinates		= f(x)
			bgy = b.evaluate(y_),		// Point B in B's local coordinates		
			agy = btoa.apply(bgy);		// Point B in A's local coordinates		= g(y)

		Vec3
			afx_gy = afx - agy;			// Difference vector between point A and point B in A's local coordinates	= f(x) - g(y)
		Real afx_gy_len = afx_gy.len();	// |f(x) - g(y)|

		if (afx_gy_len < PROXIMITY_EPS)	// If two points are close enough, just return true
			return true;
		
		Vec3 afx_1, agy_1, bgy_1;
		afx_1 = a.differentiate(x_);	// Tangent vector at point A in A's local coordinates	= f'(x)
		bgy_1 = b.differentiate(y_);	// Tangent vector at point B in B's local coordinates
		agy_1 = btoa.applyR(bgy_1);		// Tangent vector at point B in A's local coordinates	= g'(y)

		Vec3 afx_2, agy_2, bgy_2;
		afx_2 = { -afx[0], -afx[1], 0 };	// Second deriv vector at point A in A's local coordinates	= f''(x)
		bgy_2 = { -bgy[0], -bgy[1], 0 };	// Second deriv vector at point B in B's local coordinates	
		agy_2 = btoa.applyR(bgy_2);			// Second deriv vector at point B in A's local coordinates	= g''(y)

		Real afx_1_len, bgy_1_len;
		afx_1_len = afx_1.len();		// |f'(x)|
		bgy_1_len = bgy_1.len();		// |g'(y)|

		Real tmp0, tmp1, tmp2, tmp3;
		tmp0 = afx_1.dot(afx_gy);			// (f'(x) * (f(x) - g(y)))
		tmp1 = afx_1_len * afx_gy_len;		// (|f'(x)| * |f(x) - g(y)|)
		tmp2 = agy_1.dot(afx_gy);			// (g'(y) * (f(x) - g(y)))
		tmp3 = bgy_1_len * afx_gy_len;		// (|g'(y)| * |f(x) - g(y)|)

		// 1. dvec
		dvec[1] = tmp0 / tmp1;
		dvec[2] = tmp2 / tmp3;

		Real tmp4, tmp5, tmp6, tmp7, tmp8;
		tmp4 = afx_1.dot(afx_2) / afx_1_len;	// (d/dx)(|f'(x)|)
		tmp5 = afx_1.dot(afx_gy) / afx_gy_len;	// (d/dx)(|f(x) - g(y)|)
		tmp6 = agy_1.dot(agy_2) / bgy_1_len;	// (d/dy)(|g'(y)|)
		tmp7 = -agy_1.dot(afx_gy) / afx_gy_len;	// (d/dy)(|f(x) - g(y)|)
		tmp8 = afx_1.dot(agy_1);				// (f'(x) * g'(y))

		// 2. djac
		
		// ((f''(x) * (f(x) - g(y)) + f'(x) * f'(x)) * tmp1 - tmp0 * (tmp4 * |f(x) - g(y)| + |f'(x)| * tmp5)) / tmp1^2
		djac[1][1] = (((afx_2.dot(afx_gy) + SQUARE(afx_1_len)) * tmp1) - (tmp0 * (tmp4 * afx_gy_len + afx_1_len * tmp5))) / SQUARE(tmp1);			
		// ((-tmp8 * tmp1) - (tmp0 * (|f'(x)| * tmp7))) / tmp1^2
		djac[1][2] = ((-tmp8 * tmp1) - (tmp0 * (afx_1_len * tmp7))) / SQUARE(tmp1);
		// ((tmp8 * tmp3) - (tmp2 * (|g'(y)| * tmp5))) / tmp3^2
		djac[2][1] = ((tmp8 * tmp3) - (tmp2 * (bgy_1_len * tmp5))) / SQUARE(tmp3);
		// ((g''(y) * (f(x) - g(y)) - g'(y) * g'(y)) * tmp3 - tmp2 * (tmp6 * |f(x) - g(y)| + |g'(y)| * tmp7)) / tmp3^2
		djac[2][2] = (((agy_2.dot(afx_gy) - SQUARE(bgy_1_len)) * tmp3) - (tmp2 * (tmp6 * afx_gy_len + bgy_1_len * tmp7))) / SQUARE(tmp3);

		return false;
	}*/

	// For given parameters, compute objective function in inaccurate but fast way
	// @fvec : (Real) Objective function value
	// @dvec : (Square of) Objective function value
	// @djac : Jacobian matrix of objective function ( dvec )
	// @ret :  Return distance between two points on the circle, since equation becomes unstable as the distance decreases
	inline static Real binormalNRFunc(Real param[3], Real fvec[3], Real dvec[3], Real djac[3][3], const Circle& a, const Circle& b, const Transform& btoa) {
		Real
			rx = a.radius,
			ry = b.radius,
			x_ = param[1],
			y_ = param[2];
		Vec3
			afx = a.evaluate(x_),
			bgy = b.evaluate(y_),
			agy = btoa.apply(bgy);

		Vec3
			afx_gy = afx - agy;
		Real d = afx_gy.len();
		if (d < PROXIMITY_EPS)
			return d;

		Vec3 afx_1, agy_1, bgy_1;
		afx_1 = a.differentiate(x_, 1);
		bgy_1 = b.differentiate(y_, 1);
		agy_1 = btoa.applyR(bgy_1);

		Vec3 afx_2, bgy_2, agy_2;
		afx_2 = { -afx[0], -afx[1], 0 };
		bgy_2 = { -bgy[0], -bgy[1], 0 };
		agy_2 = btoa.applyR(bgy_2);

		// 1. dvec
		dvec[1] = afx_1.dot(afx_gy);
		dvec[2] = agy_1.dot(afx_gy);
		fvec[1] = dvec[1] / afx_1.len() / d;
		fvec[2] = dvec[2] / bgy_1.len() / d;

		// 2. djac
		djac[1][1] = afx_2.dot(afx_gy) + afx_1.dot(afx_1);
		djac[1][2] = -afx_1.dot(agy_1);
		djac[2][1] = -djac[1][2];
		djac[2][2] = agy_2.dot(afx_gy) - bgy_1.dot(bgy_1);
		return d;
	}
	// For given parameters, compute objective function in accurate but slow way
	// @dvec : Objective function value
	// @djac : Jacobian matrix of objective function ( dvec )
	// @ret :	If two points are close enough, |f(x) - g(y)| converges to zero and objective function becomes unstable. 
	//			Therefore, in those cases, we return TRUE to indicate that those points are close
	/*inline static bool binormalNRFunc(Real param[3], Real dvec[3], Real djac[3][3], const Circle& a, const Circle& b, const Transform& btoa) {
		Real u = param[1];
		Real v = param[2];

		Vec3 A = a.evaluate(u);
		Vec3 B = btoa.apply(b.evaluate(v));
		Vec3 D = A - B;
		Real d = D.len();
		if (d < PROXIMITY_EPS)
			return true;

		Vec3 Au = a.differentiate(u, 1);
		Vec3 Bv = btoa.applyR(b.differentiate(v, 1));
		Real AuLen = Au.len();
		Real BvLen = Bv.len();

		Vec3 Auu = a.differentiate(u, 2);
		Vec3 Bvv = btoa.applyR(b.differentiate(v, 2));

		Vec3 Du = D.cross(Au.cross(D)) / (d * d * d);
		Vec3 Dv = D.cross(Bv.cross(D)) / (d * d * d) * -1;

		dvec[1] = Au.dot(D) / AuLen / d;
		dvec[2] = Bv.dot(D) / BvLen / d;

		djac[1][1] = (Au.cross(Auu.cross(Au)) / (AuLen * AuLen * AuLen)).dot(D / d) + (Au / AuLen).dot(Du);
		djac[1][2] = (Au / AuLen).dot(Dv);
		djac[2][1] = (Bv / BvLen).dot(Du);
		djac[2][2] = (Bv.cross(Bvv.cross(Bv)) / (BvLen * BvLen * BvLen)).dot(D / d) + (Bv / BvLen).dot(Dv);
		return false;
	}*/
	inline static bool binormalNR(const Circle& a, const Circle& b, const Transform& btoa, Real param[3], Real eps) {
		int k, i, indx[3];
		Real errx, errf, d, dvec[3], fvec[3], djac[3][3], p[3], dist;
		bool largerEps = false;

		for (k = 0; k < NR_MAX_ITER; k++) {
			// User function supplies function values at [x] in [dvec] and Jacobian matrix in [djac].
			dist = binormalNRFunc(param, fvec, dvec, djac, a, b, btoa);

			if (dist < PROXIMITY_EPS)
				return true;											// If two points meet, just regard them as binormal
			else if (dist < PROXIMITY_EPS2 && !largerEps) {
				eps *= 100;												// If two points are close enough so that fvec values become unstable,
				largerEps = true;										// we mitigate the [ eps ] to find binormal easily
			}

			errf = 0.0;
			for (i = 1; i <= 2; i++)
				errf += fabs(fvec[i]);			// Check function convergence.
			if (errf <= eps)
				return true;

			for (i = 1; i <= 2; i++)
				p[i] = -dvec[i];										// Right-hand side of linear equations.
			binormalLUdcmp(djac, indx, &d);								// Solve linear equations using LU decomposition.
			binormalLUbksb(djac, indx, p);
			errx = 0.0;
			for (i = 1; i <= 2; i++) {									// Check root convergence.
				errx += fabs(p[i]);										// Update solution.
				param[i] += p[i];
			}
		}
		return false;
	}

	/* Subroutine */

	// Subroutine of [ exceptionB ] : We have to solve quadratic equation to find points on circle [ b ] that goes through axis of circle [ a ]
	inline static void solveQuadraticEquation(const Real c[3], Real roots[2], int& num) {
		// Solve second degree equation : c[0] * x^2 + 2 * c[1] * x + c[0] = 0
		num = 0;
		Real det = c[1] * c[1] - c[0] * c[2];
		if (det < 0) {
			static const Real eps = 1e-10;	// For numerical stability
			if (det > -eps) det = 0;
			else return;
		}
		det = sqrt(det);
		roots[num++] = (-c[1] + det) / c[0];
		if (det > 0)
			roots[num++] = (-c[1] - det) / c[0];
	}
	// Add new binormal to binormal vector without duplication
	inline static void safeAddBinormal(const CircularArc& a, const CircularArc& b, const Transform& btoa, Real paramA, Real paramB, std::vector<CircleBinormal::Binormal>& bins) {
		Real x, y, diffX, diffY;
		bool duplicate, sameX, sameY;
		
		duplicate = false;
		x = piDomain::regularize(paramA);
		y = piDomain::regularize(paramB);
		for (auto& bin : bins) {
			diffX = fabs(x - bin.paramA);
			diffY = fabs(y - bin.paramB);
			sameX = (diffX < DOMAIN_EPS) || (diffX > PI20 - DOMAIN_EPS);
			sameY = (diffY < DOMAIN_EPS) || (diffY > PI20 - DOMAIN_EPS);
			if (sameX && sameY) {
				duplicate = true;
				break;
			}
		}
		if (!duplicate) {
			CircleBinormal::Binormal nbin;
			nbin.paramA = x;
			nbin.paramB = y;
			nbin.pointA = a.evaluate(x);
			nbin.pointB = b.evaluate(y);
			nbin.distance = btoa.apply(nbin.pointB).dist(nbin.pointA);

			// Type check
			if (nbin.pointA.dist(btoa.T) < PROXIMITY_EPS) {
				Vec3 atan = a.differentiate(x, 1);
				atan.normalize();
				Real det = atan.dot(Vec3{ btoa.R[0][2], btoa.R[1][2], btoa.R[2][2] });
				if (fabs(det) > 1 - 1e-10) nbin.type = 3;
				else nbin.type = 0;
			}
			else nbin.type = 0;
			bins.push_back(nbin);
		}
	}
	// Find potential binormals for given [ cosB ] parameter
	inline static void findPotBinormal(const CircularArc& a, const CircularArc& b, const Transform& btoa, Real cosB, CircleBinormal::Binormal bin[4], int& num) {
		Real t1, t2, d;
		Vec3 bpt, bptA;
		num = 0;

		// We can use 2 different schemes to find potential binormal lines here :
		// 1. Project [ bpt ] on arc [ a ] 
		// 2. Intersect plane created by axis of arc [ b ] and vector [ bpt ] with arc [ a ]
		// Since scheme 1 is cheaper than 2, we first use 1 and use 2 only if scheme 1 is unstable :
		// When [ bptA ] is near to the origin, so that projection becomes unstable
		const static Real thresd = 0.2;	// [ bptA ] must be far away than this to recognised as stable

		t2 = acos(cosB);
		if (b.domain.has(t2)) {
			bpt = b.evaluate(t2);
			bptA = btoa.apply(bpt);

			d = sqrt(bptA[0] * bptA[0] + bptA[1] * bptA[1]);

			if (d / a.radius > thresd) {
				// Scheme 1
				a.Circle::findMinDistParam(bptA, t1);
				if (a.domain.has(t1)) {
					bin[num].paramA = t1;
					bin[num].paramB = t2;
					num++;
				}

				t1 = piDomain::regularize(t1 + PI);
				if (a.domain.has(t1)) {
					bin[num].paramA = t1;
					bin[num].paramB = t2;
					num++;
				}
			}
			else {
				// Scheme 2
				bpt = b.differentiate(t2, 1);
				bptA = btoa.applyR(bpt);
				Real intParam[2];
				int intNum;
				intersect(a, bptA, btoa.T, intParam, intNum);

				for (int i = 0; i < intNum; i++) {
					bin[num].paramA = intParam[i];
					bin[num].paramB = t2;
					num++;
				}
			}
		}
		
		t2 = piDomain::regularize(PI20 - t2);
		if (b.domain.has(t2)) {
			bpt = b.evaluate(t2);
			bptA = btoa.apply(bpt);
			
			d = sqrt(bptA[0] * bptA[0] + bptA[1] * bptA[1]);

			if (d / a.radius > thresd) {
				// Scheme 1
				a.Circle::findMinDistParam(bptA, t1);
				if (a.domain.has(t1)) {
					bin[num].paramA = t1;
					bin[num].paramB = t2;
					num++;
				}

				t1 = piDomain::regularize(t1 + PI);
				if (a.domain.has(t1)) {
					bin[num].paramA = t1;
					bin[num].paramB = t2;
					num++;
				}
			}
			else {
				// Scheme 2
				bpt = b.differentiate(t2, 1);
				bptA = btoa.applyR(bpt);
				Real intParam[2];
				int intNum;
				intersect(a, bptA, btoa.T, intParam, intNum);

				for (int i = 0; i < intNum; i++) {
					bin[num].paramA = intParam[i];
					bin[num].paramB = t2;
					num++;
				}
			}
		}
	}

	// Exceptions
	bool CircleBinormal::exceptionA(const CircularArc& a, const CircularArc& b, const Transform& btoa) {
		Real d = sqrt(btoa.T[0] * btoa.T[0] + btoa.T[1] * btoa.T[1]);
		if (d > PROXIMITY_EPS) 
			return false;

		Real axisD = fabs(btoa.R[2][2]);
		if (axisD < 1 - PROXIMITY_EPS) 
			return false;

		return true;
	}
	inline static void arcPlaneIntersection(const CircularArc& arc, Real u, const Vec3& normal, Real p, Vec3 pt[4], Real param[4], int& num) {
		u = piDomain::regularize(u);
		if (!arc.domain.has(u)) return;

		pt[num] = arc.evaluate(u);
		param[num] = u;
		num++;
	}
	inline static void arcPlaneIntersection(const CircularArc& arc, const Vec3& normal, Real p, Vec3 pt[4], Real param[4], int& num) {
		// Intersection between arc and plane defined by [ normal ] vector and [ p ]
		// At most 4 intersection points are created
		num = 0;
		if (normal[0] == 0) {
			Real sinu = p / (normal[1] * arc.radius);
			if (fabs(sinu) > 1) {
				if (fabs(sinu) < 1 + 1e-10) {
					if (sinu > 0) sinu = 1;
					else sinu = -1;
				}
				else
					return;
			}
			Real u = asin(sinu);
			arcPlaneIntersection(arc, u, normal, p, pt, param, num);
			if (fabs(sinu) != 1)
				arcPlaneIntersection(arc, PI - u, normal, p, pt, param, num);
		}
		else if (normal[1] == 0) {
			Real cosu = p / (normal[0] * arc.radius);
			if (fabs(cosu) > 1) {
				if (fabs(cosu) < 1 + 1e-10) {
					if (cosu > 0) cosu = 1;
					else cosu = -1;
				}
				else
					return;
			}
			Real u = acos(cosu);
			arcPlaneIntersection(arc, u, normal, p, pt, param, num);
			if (fabs(cosu) != 1)
				arcPlaneIntersection(arc, PI20 - u, normal, p, pt, param, num);
		}
		else {
			Real c[3];
			c[0] = SQ(arc.radius) * (SQ(normal[0]) + SQ(normal[1]));
			c[1] = -p * normal[0] * arc.radius;
			c[2] = SQ(p) - SQ(normal[1] * arc.radius);

			Real cosu[2];
			int cosnum;
			solveQuadraticEquation(c, cosu, cosnum);

			for (int i = 0; i < cosnum; i++) {
				if (fabs(cosu[i]) > 1) {
					if (fabs(cosu[i]) < 1 + 1e-10) {
						if (cosu > 0) cosu[i] = 1;
						else cosu[i] = -1;
					}
					else
						continue;
				}
				Real u = acos(cosu[i]);
				arcPlaneIntersection(arc, u, normal, p, pt, param, num);
				if (fabs(cosu[i]) != 1)
					arcPlaneIntersection(arc, PI20 - u, normal, p, pt, param, num);
			}
		}
	}
	// For given value [ bu ], check if there are points on circle [ a ] that form binormal line with [ B(bu) ]
	inline static void subroutineExceptionB2(const CircularArc& a, const CircularArc& b, const Transform& btoa, Real bu, std::vector<CircleBinormal::Binormal>& bins) {
		Vec3 bpt, bptA, btanA;
		bpt = b.evaluate(bu);
		bptA = btoa.apply(bpt);

		// [ bptA ] must be on circle A's axis
		Real d = sqrt(bptA[0] * bptA[0] + bptA[1] * bptA[1]);
		if (d > 1e-10) return;

		btanA = btoa.applyR(b.differentiate(bu, 1));
		btanA.normalize();
		if (fabs(btanA[2]) > 1 - 1e-10) {
			// When [ btanA ] is parallel to circle A's axis
			if (fabs(bptA[2]) < 1e-10) {
				// If [ bptA ] equals to the center of circle B
				CircleBinormal::Binormal bin;
				bin.type = 2;
				bin.pointB = bpt;
				bin.paramB = bu;
				bin.distance = (bptA - (a.evaluate(0))).len();
				bins.push_back(bin);
				return;
			}
			else 
				return;
		}
		else {
			//Vec3 apt[4];
			//arcPlaneIntersection(a, btanA, btanA.dot(btoa.T), apt, aparam, anum);
			Real aparam[2];
			int anum;
			intersect(a, btanA, btoa.T, aparam, anum);

			for (int i = 0; i < anum; i++) {
				Real param[3] = { 0, aparam[i], bu };
				bool success = binormalNR(a, b, btoa, param, 1e-10);
				// Since refinement is done, just use original domain, not extended one
				if (success && a.domain.has(param[1]) && b.domain.has(param[2]))
					safeAddBinormal(a, b, btoa, param[1], param[2], bins);
			}
		}
	}
	// For given value [ bTriU ], compute parameter [ u ] for circle [ b ] that goes through axis of circle [ a ]
	// @isCos : Trus if [ bTriU ] is cosine value of parameter [ u ]. False if it is sine value
	inline static void subroutineExceptionB(const CircularArc& a, const CircularArc& b, const Transform& btoa, Real bTriU, std::vector<CircleBinormal::Binormal>& bins, bool isCos) {
		Real abs = fabs(bTriU);
		if (abs > 1) {
			if (abs < 1 + 1e-10) {
				if (bTriU > 0) bTriU = 1;
				else bTriU = -1;
			}
			else
				return;
		}
		if (isCos) {
			Real u0 = acos(bTriU);
			if(b.domain.has(u0))
				subroutineExceptionB2(a, b, btoa, u0, bins);
			if (u0 != PI && u0 != 0 && b.domain.has(PI20 - u0)) 
				subroutineExceptionB2(a, b, btoa, PI20 - u0, bins);
		}
		else {
			Real u0 = asin(bTriU);
			if(b.domain.has(u0))
				subroutineExceptionB2(a, b, btoa, u0, bins);
			if (u0 != PI05 && u0 != -PI05 && b.domain.has(PI - u0)) 
				subroutineExceptionB2(a, b, btoa, PI - u0, bins);
		}
	}
	void CircleBinormal::exceptionB(const CircularArc& a, const CircularArc& b, const Transform& btoa, std::vector<Binormal>& bins) {
		Vec3 uVec, vVec;
		uVec = { btoa.R[0][0], btoa.R[1][0], btoa.R[2][0] };
		vVec = { btoa.R[0][1], btoa.R[1][1], btoa.R[2][1] };

		// Find point on [ b ] that goes through axis of [ a ]
		// We have to find [ t ] that satisfies following equations :
		//		btoa.T[0] + (rcos(t) * uVec[0]) + (rsin(t) * vVec[0]) = 0 ( r = radius of [ b ] )
		//		btoa.T[1] + (rcos(t) * uVec[1]) + (rsin(t) * vVec[1]) = 0 ( r = radius of [ b ] )
		if (uVec[1] == 0 && vVec[1] == 0) {
			if (btoa.T[1] != 0) 
				// Since second equation does not satisfy...
				return;

			if (vVec[0] == 0) {
				Real cos = -btoa.T[0] / uVec[0] / b.radius;
				subroutineExceptionB(a, b, btoa, cos, bins, true);
			}
			else {
				Real c[3];
				c[0] = (uVec[0] * uVec[0]) / (vVec[0] * vVec[0]) + 1;
				c[1] = (btoa.T[0] * uVec[0]) / (b.radius * vVec[0] * vVec[0]);
				c[2] = (btoa.T[0] * btoa.T[0]) / (b.radius * b.radius * vVec[0] * vVec[0]) - 1;

				Real cosu[2];
				int num;
				// Solve quadratic equation
				solveQuadraticEquation(c, cosu, num);

				for (int i = 0; i < num; i++) 
					subroutineExceptionB(a, b, btoa, cosu[i], bins, true);
			}
		}
		else if (uVec[1] == 0) {
			Real mult = (vVec[0] / vVec[1]);
			Real cosu = -(btoa.T[0] - btoa.T[1] * mult) / (b.radius * uVec[0]);
			subroutineExceptionB(a, b, btoa, cosu, bins, true);
		}
		else {
			Real mult = (uVec[0] / uVec[1]);
			Real sinu = -(btoa.T[0] - btoa.T[1] * mult) / (b.radius * (vVec[0] - vVec[1] * mult));
			subroutineExceptionB(a, b, btoa, sinu, bins, false);
		}
	}

	// Solve
	void CircleBinormal::subroutine(const CircularArc& a, const CircularArc& b, const Transform& btoa, const std::vector<Domain>& bCosDomains, std::vector<Binormal>& bins, bool refine, Real precision) {
		// Assume [ a ] is located on XY plane.
		Vec3 C, U, V;	// Center, orthonormal direction of [ b ] in [ a ]'s local coordinates.
		C = btoa.T;
		U = { btoa.R[0][0], btoa.R[1][0], btoa.R[2][0] };
		V = { btoa.R[0][1], btoa.R[1][1], btoa.R[2][1] };

		Real coefA[6];
		coefA[0] = C.dot(U);
		coefA[1] = C.dot(V);
		coefA[2] = C[2];
		coefA[3] = U[2];
		coefA[4] = V[2];
		coefA[5] = C.dot(C);

		Real r0, r1, r0_2, r0_3, r0_4, r1_2_m;
		r0 = b.radius;
		r1 = a.radius;
		r0_2 = SQ(r0);
		r0_3 = r0_2 * r0;
		r0_4 = r0_3 * r0;
		r1_2_m = -r1 * r1;

		Real s[9];
		s[3] = SQ(coefA[0]);
		s[4] = SQ(coefA[1]);
		s[5] = SQ(coefA[3]);
		s[6] = SQ(coefA[4]);
		s[7] = coefA[0] * coefA[1];
		s[8] = coefA[3] * coefA[4];
		s[0] = coefA[5] - SQ(coefA[2]) + r0_2;
		s[1] = coefA[1] - coefA[2] * coefA[4];
		s[2] = coefA[0] - coefA[2] * coefA[3];

		Real coefB[4];
		coefB[0] = s[3] - s[4];
		coefB[1] = s[5] - s[6];
		coefB[2] = 4.0 * s[7] * s[8];
		coefB[3] = 2.0 * s[1] * s[7];

		Real c[7];
		c[0] = s[5] - s[6];
		c[1] = s[2] * s[8];
		c[2] = s[8] * c[0];
		c[3] = s[1] * s[8];
		c[4] = 4.0 * SQ(s[8]) - SQ(c[0]);
		c[5] = s[2] * c[0];
		c[6] = SQ(s[2]);

		Real m[9];
		m[0] = 2.0 * r0_3 * s[1] * s[3];
		m[1] = -2.0 * r0_2 * (s[0] * s[7] + r0_2 * (s[3] * s[8] - s[6] * s[7]));
		m[2] = -2.0 * r0_3 * (s[1] * coefB[0] + 2.0 * s[2] * s[7]);
		m[3] = 2.0 * r0_4 * (s[8] * coefB[0] + s[7] * coefB[1]);
		m[4] = r0_2 * s[3] * (s[0] - r0_2 * s[6]);
		m[5] = 2.0 * r0_3 * (s[2] * s[3] - coefB[3]);
		m[6] = r0_2 * (-s[0] * coefB[0] + r0_2 * (-s[3] * coefB[1] + s[6] * coefB[0] + coefB[2]));
		m[7] = -2.0 * r0_3 * (s[2] * coefB[0] - coefB[3]);
		m[8] = r0_4 * (coefB[0] * coefB[1] - coefB[2]);

		Real n[9];
		n[0] = -2.0 * r0_3 * c[1];
		n[1] = 2.0 * r0_2 * (-s[1] * s[2] + r0_2 * c[2]);
		n[2] = 2.0 * r0_3 * (2.0 * c[1] + s[1] * c[0]);
		n[3] = -4.0 * r0_4 * c[2];
		n[4] = r0_2 * (c[6] + r0_2 * SQ(s[8]));
		n[5] = 2.0 * r0_3 * (-c[5] + c[3]);
		n[6] = r0_2 * (-c[6] + SQ(s[1]) - r0_2 * c[4]);
		n[7] = 2.0 * r0_3 * (c[5] - 2.0 * c[3]);
		n[8] = r0_4 * c[4];

		for (int i = 0; i < 9; i++) n[i] *= r1_2_m;

		Real x[9];
		for (int i = 0; i < 9; i++) x[i] = m[i] + n[i];

		Real z[6];
		z[0] = 2.0 * x[0] * x[1];
		z[1] = 2.0 * x[0] * x[2];
		z[2] = 2.0 * x[0] * x[3];
		z[3] = 2.0 * x[1] * x[2];
		z[4] = 2.0 * x[1] * x[3];
		z[5] = 2.0 * x[2] * x[3];

		Real coef[9];
		coef[0] = SQ(x[4]) - SQ(x[0]);
		coef[1] = 2.0 * x[4] * x[5] - z[0];
		coef[2] = 2.0 * x[4] * x[6] + SQ(x[5]) - SQ(x[1]) - z[1] + SQ(x[0]);
		coef[3] = 2.0 * (x[4] * x[7] + x[5] * x[6]) - z[2] - z[3] + z[0];
		coef[4] = 2.0 * (x[4] * x[8] + x[5] * x[7]) + SQ(x[6]) - SQ(x[2]) - z[4] + SQ(x[1]) + z[1];
		coef[5] = 2.0 * (x[5] * x[8] + x[6] * x[7]) - z[5] + z[2] + z[3];
		coef[6] = 2.0 * x[6] * x[8] + SQ(x[7]) - SQ(x[3]) + SQ(x[2]) + z[4];
		coef[7] = 2.0 * x[7] * x[8] + z[5];
		coef[8] = SQ(x[8]) + SQ(x[3]);

		// @BUGFIX : Make coefficients stable
		Real coefMin = maxDouble, coefMax = 0, coefAvg;
		for (int i = 0; i < 9; i++) {
			if (fabs(coef[i]) > coefMax)
				coefMax = fabs(coef[i]);
			if (fabs(coef[i]) < coefMin)
				coefMin = fabs(coef[i]);
		}
		coefAvg = (coefMin + coefMax) * 0.5;
		if (coefAvg > 0)
			for (int i = 0; i < 9; i++)
				coef[i] /= coefAvg;

		Real roots[16];
		int rootNum;

		// @ Since BezierPolynomial automatically culls out duplicate roots, we do not have to test it
		solveBP(coef, bCosDomains, roots, rootNum);

		// Slightly extend domain for stability
		CircularArc aCopy = a, bCopy = b;
		if (aCopy.domain.width() < PI20 - PROJECTION_EPS * 2) aCopy.domain.set(aCopy.domain.beg() - PROJECTION_EPS, aCopy.domain.end() + PROJECTION_EPS);
		else aCopy.domain.set(0, PI20);
		if (bCopy.domain.width() < PI20 - PROJECTION_EPS * 2) bCopy.domain.set(bCopy.domain.beg() - PROJECTION_EPS, bCopy.domain.end() + PROJECTION_EPS);
		else bCopy.domain.set(0, PI20);

		Binormal tmpBins[4];
		for (int i = 0; i < rootNum; i++) {
			if (fabs(roots[i]) > 1 + TOL) continue;
			else if (roots[i] > 1)	roots[i] = 1;
			else if (roots[i] < -1) roots[i] = -1;

			int tnum;
			findPotBinormal(aCopy, bCopy, btoa, roots[i], tmpBins, tnum);
			for (int j = 0; j < tnum; j++) {
				if (refine) {
					Real param[3] = { 0, tmpBins[j].paramA, tmpBins[j].paramB };
					bool success = binormalNR(a, b, btoa, param, precision);
					// Since refinement is done, just use original domain, not extended one
					if (success && a.domain.has(param[1]) && b.domain.has(param[2]))
						safeAddBinormal(a, b, btoa, param[1], param[2], bins);
				}
				else
					safeAddBinormal(a, b, btoa, tmpBins[j].paramA, tmpBins[j].paramB, bins);
			}
		}
	}
	void CircleBinormal::solve(const Circle& a, const Circle& b, const Transform& tA, const Transform& tB, std::vector<Binormal>& bins, bool refine, Real precision) {
		Transform btoa = Transform::connect(tB, tA);
		solve(a, b, btoa, bins, refine, precision);
	}
	void CircleBinormal::solve(const CircularArc& a, const CircularArc& b, const Transform& tA, const Transform& tB, std::vector<Binormal>& bins, bool refine, Real precision) {
		Transform btoa = Transform::connect(tB, tA);
		solve(a, b, btoa, bins, refine, precision);
	}
	void CircleBinormal::solve(const Circle& a, const Circle& b, const Transform& btoa, std::vector<Binormal>& bins, bool refine, Real precision) {
		CircularArc arcA, arcB;
		arcA.radius = a.radius;
		arcB.radius = b.radius;
		arcA.domain = piDomain::create(0, PI20);
		arcB.domain = piDomain::create(0, PI20);

		solve(arcA, arcB, btoa, bins, refine, precision);
	}
	void CircleBinormal::solve(const CircularArc& a, const CircularArc& b, const Transform& btoa, std::vector<Binormal>& bins, bool refine, Real precision) {
		// For numerical stability, scale circles to make average radius to be 1.0
		Real avgRadius = (a.radius + b.radius) * 0.5;

		Transform nbtoa = btoa;
		nbtoa.T /= avgRadius;

		CircularArc arcA = a, arcB = b;
		arcA.radius = a.radius / avgRadius;
		arcB.radius = b.radius / avgRadius;

		bins.clear();
		/* Check for exceptional cases */
		// Exception 1
		if (exceptionA(arcA, arcB, nbtoa)) {
			Binormal bin;
			bin.type = 1;
			bins.push_back(bin);
			return;
		}

		// Exception 2
		exceptionB(arcA, arcB, nbtoa, bins);

		// Set domains to solve 8-th degree polynomial
		Domain aCosDom, bCosDom;
		Real aBegCos = cos(a.domain.beg()), aEndCos = cos(a.domain.end());
		Real bBegCos = cos(b.domain.beg()), bEndCos = cos(b.domain.end());

		Real beg, end;
		if (a.domain.has(PI)) beg = -1.0;
		else beg = (aBegCos < aEndCos) ? aBegCos : aEndCos;
		if (a.domain.has(0.0)) end = 1.0;
		else end = (aBegCos > aEndCos) ? aBegCos : aEndCos;
		aCosDom.set(beg, end);

		if (b.domain.has(PI)) beg = -1.0;
		else beg = (bBegCos < bEndCos) ? bBegCos : bEndCos;
		if (b.domain.has(0.0)) end = 1.0;
		else end = (bBegCos > bEndCos) ? bBegCos : bEndCos;
		bCosDom.set(beg, end);

		// Solve 8-th degree polynomial of the circle with smaller domain
		if (aCosDom.width() >= bCosDom.width())
			subroutine(arcA, arcB, nbtoa, { bCosDom }, bins, refine, precision);
		else {
			subroutine(arcB, arcA, nbtoa.inverse(), { aCosDom }, bins, refine, precision);
			for (auto& bin : bins) {
				std::swap(bin.paramA, bin.paramB);
				std::swap(bin.pointA, bin.pointB);
			}
		}
		
		// Recover real radius and distance
		for (auto& bin : bins) {
			bin.pointA *= avgRadius;
			bin.pointB *= avgRadius;
			bin.distance *= avgRadius;
		}
	}

	void CircleBinormal::solve(const Circle& a, const Circle& b, const Transform& btoa, const std::vector<piDomain>& bDomain, std::vector<Binormal>& bins, bool refine, Real precision) {
		// For numerical stability, scale circles to make average radius to be 1.0
		Real avgRadius = (a.radius + b.radius) * 0.5;

		Transform nbtoa = btoa;
		nbtoa.T /= avgRadius;

		CircularArc arcA, arcB;
		arcA.radius = a.radius / avgRadius;
		arcB.radius = b.radius / avgRadius;
		arcA.domain.set(0, PI20);
		arcB.domain.set(0, PI20);	// We use [ bDomain ] instead

		bins.clear();
		// Check for exceptional case 
		// Exception 1
		if (exceptionA(arcA, arcB, nbtoa)) {
			Binormal bin;
			bin.type = 1;
			bins.push_back(bin);
			return;
		}

		// Exception 2
		exceptionB(arcA, arcB, nbtoa, bins);	// @TODO : It is not problem only for torus binormal with gaussmap...

		std::vector<Domain> bCosDomains;
		bCosDomains.reserve(bDomain.size());

		Domain cosDomain;
		for (auto& domain : bDomain) {
			Real beg, end;
			Real bBegCos = cos(domain.beg()), bEndCos = cos(domain.end());
			if (domain.has(PI)) beg = -1;
			else beg = (bBegCos < bEndCos) ? bBegCos : bEndCos;
			if (domain.has(0) || domain.has(PI20)) end = 1;
			else end = (bBegCos > bEndCos) ? bBegCos : bEndCos;
			cosDomain.set(beg, end);
			bCosDomains.push_back(cosDomain);
		}
		subroutine(arcA, arcB, nbtoa, bCosDomains, bins, refine, precision);

		for (auto& bin : bins) {
			bin.pointA *= avgRadius;
			bin.pointB *= avgRadius;
			bin.distance *= avgRadius;
		}
	}
	void CircleBinormal::solve(const CircularArc& a, const CircularArc& b, const Transform& btoa, const std::vector<piDomain>& bDomain, std::vector<Binormal>& bins, bool refine, Real precision) {
		// For numerical stability, scale circles to make average radius to be 1.0
		Real avgRadius = (a.radius + b.radius) * 0.5;

		Transform nbtoa = btoa;
		nbtoa.T /= avgRadius;

		CircularArc arcA = a, arcB = b;
		arcA.radius = a.radius / avgRadius;
		arcB.radius = b.radius / avgRadius;

		bins.clear();
		// Check for exceptional case 
		// Exception 1
		if (exceptionA(arcA, arcB, nbtoa)) {
			Binormal bin;
			bin.type = 1;
			bins.push_back(bin);
			return;
		}

		// Exception 2
		exceptionB(arcA, arcB, nbtoa, bins);	// @TODO : It is not problem only for torus binormal with gaussmap...

		std::vector<Domain> bCosDomains;
		bCosDomains.reserve(bDomain.size());

		Domain cosDomain;
		for (auto& domain : bDomain) {
			Real beg, end;
			Real bBegCos = cos(domain.beg()), bEndCos = cos(domain.end());
			if (domain.has(PI)) beg = -1;
			else beg = (bBegCos < bEndCos) ? bBegCos : bEndCos;
			if (domain.has(0) || domain.has(PI20)) end = 1;
			else end = (bBegCos > bEndCos) ? bBegCos : bEndCos;
			cosDomain.set(beg, end);
			bCosDomains.push_back(cosDomain);
		}
		subroutine(arcA, arcB, nbtoa, bCosDomains, bins, refine, precision);

		std::vector<Binormal> nbins;
		nbins.reserve(bins.size());
		for (auto& bin : bins) {
			if (a.domain.has(bin.paramA) && b.domain.has(bin.paramB)) {
				bin.pointA *= avgRadius;
				bin.pointB *= avgRadius;
				bin.distance *= avgRadius;
				nbins.push_back(bin);
			}
		}
		bins = nbins;
	}
	// Brute Solve
	static void bruteSolveAtGivenParam(const CircularArc& a, const CircularArc& b, Real aparam, Real bparam, const Transform& btoa, std::vector<CircleBinormal::Binormal>& bins, Real precision);
	static void bruteSolveAtGivenBParam(const CircularArc& a, const CircularArc& b, Real bparam, const Transform& btoa, std::vector<CircleBinormal::Binormal>& bins, Real precision);
	void CircleBinormal::bruteSolve(const CircularArc& a, const CircularArc& b, const Transform& tA, const Transform& tB, std::vector<Binormal>& bins, Real precision, int sampleNum) {
		Transform btoa = Transform::connect(tB, tA);
		const auto& bDomain = b.domain;
		//Binormal bin;

		for (int i = 0; i <= sampleNum; i++) {
			Real bt = bDomain.beg() + bDomain.width() * ((Real)i / (Real)sampleNum);
			bruteSolveAtGivenBParam(a, b, bt, btoa, bins, precision);
		}
	}

	static void bruteSolveAtGivenParam(const CircularArc& a, const CircularArc& b, Real aparam, Real bparam, const Transform& btoa, std::vector<CircleBinormal::Binormal>& bins, Real precision) {
		Real param[3], fvec[3], dvec[3], djac[3][3];
		param[1] = aparam;
		param[2] = bparam;
		Real dist = binormalNRFunc(param, fvec, dvec, djac, a, b, btoa);
		if (dist < PROXIMITY_EPS)
			safeAddBinormal(a, b, btoa, aparam, bparam, bins);
		else {
			Real det = fabs(fvec[1]) + fabs(fvec[1]);
			if (det < 0.2) {	// Only if current parameters give approximate binormals
				bool success = binormalNR(a, b, btoa, param, precision);
				if (success && a.domain.has(param[1]) && b.domain.has(param[2]))
					safeAddBinormal(a, b, btoa, param[1], param[2], bins);
			}
		}
	}
	static void bruteSolveAtGivenBParam(const CircularArc& a, const CircularArc& b, Real bparam, const Transform& btoa, std::vector<CircleBinormal::Binormal>& bins, Real precision) {
		auto bpt = b.evaluate(bparam);
		bpt = btoa.apply(bpt);

		Real minp, maxp;
		a.Circle::findExtDistParam(bpt, minp, maxp);

		bruteSolveAtGivenParam(a, b, minp, bparam, btoa, bins, precision);
		bruteSolveAtGivenParam(a, b, maxp, bparam, btoa, bins, precision);
	}
}