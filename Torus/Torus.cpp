#include "Torus.h"

#define NR_MAX_ITER		20
#define LUDCMP_EPS		1.0e-20

namespace MN {
	inline static void ludcmp(Real a[3][3], int* idx, Real* d) {
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
				a[j][j] = LUDCMP_EPS;
			if (j != 2) {
				dum = 1.0 / (a[j][j]);
				for (i = j + 1; i <= 2; i++) a[i][j] *= dum;
			}
		}
	}
	inline static void lubksb(Real a[3][3], int* idx, Real b[]) {
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
	void Torus::minDistParamRefine(const Vec3& pt, Real& u, Real& v) const {
		// Use KKT condition & Simple gradient descent
		const static Real eps = 1e-11;
		Vec3 T, Tu, Tv, Tuu, Tuv, Tvv, diff;
		Real D;

		T = evaluate(u, v);
		diff = T - pt;
		D = sqrt(diff.dot(diff));

		Real du, dv, pu, pv;
		Real addu, addv;

		int cnt = 0;
		while (true) {
			cnt++;
			Tu = differentiate(u, v, 1, 0);
			Tv = differentiate(u, v, 0, 1);
			Tuu = differentiate(u, v, 2, 0);
			Tuv = differentiate(u, v, 1, 1);
			Tvv = differentiate(u, v, 0, 2);

			du = Tu.dot(diff) / D;
			dv = Tv.dot(diff) / D;

			// KKT
			pu = u;
			pv = v;

			Real dvec[3], djac[3][3];
			dvec[1] = du * D;
			dvec[2] = dv * D;

			djac[1][1] = Tuu.dot(diff) + Tu.dot(Tu);
			djac[1][2] = Tuv.dot(diff) + Tu.dot(Tv);
			djac[2][1] = djac[1][2];
			djac[2][2] = Tvv.dot(diff) + Tv.dot(Tv);

			Real p[3], d;
			int indx[3];
			for (int i = 1; i <= 2; i++)
				p[i] = -dvec[i];								// Right-hand side of linear equations.
			ludcmp(djac, indx, &d);								// Solve linear equations using LU decomposition.
			lubksb(djac, indx, p);
			u = pu + p[1];
			v = pv + p[2];

			// Test KKT result
			T = evaluate(u, v);
			diff = T - pt;
			Real tmpD = diff.len();
			if (tmpD <= D) {
				if (tmpD > D - eps)
					return;
				D = tmpD;
				continue;
			}

			// GR
			Real maxAlpha = 1.0;

			Real startAlpha = maxAlpha * 0.1;
			while (true) {
				cnt++;
				Real tmpAlpha = startAlpha;
				u = pu - tmpAlpha * du;
				v = pv - tmpAlpha * dv;

				T = evaluate(u, v);
				diff = T - pt;
				Real tmpD = sqrt(diff.dot(diff));
				if (tmpD <= D) {
					if (tmpD > D - eps)
						return;
					D = tmpD;
					break;
				}
				startAlpha *= 0.1;
			}
		}
	}
	void TorusPatch::minDistParamRefine(const Vec3& pt, Real& u, Real& v) const {
		// Use KKT condition & Simple gradient descent
		const static Real eps = 1e-11;
		Vec3 T, Tu, Tv, Tuu, Tuv, Tvv, diff;
		Real D;

		T = evaluate(u, v);
		diff = T - pt;
		D = sqrt(diff.dot(diff));

		Real du, dv, pu, pv;
		Real addu, addv;

		int cnt = 0;
		while (true) {
			cnt++;
			Tu = differentiate(u, v, 1, 0);
			Tv = differentiate(u, v, 0, 1);
			Tuu = differentiate(u, v, 2, 0);
			Tuv = differentiate(u, v, 1, 1);
			Tvv = differentiate(u, v, 0, 2);

			du = Tu.dot(diff) / D;
			dv = Tv.dot(diff) / D;
			
			pu = u;
			pv = v;

			// KKT
			bool kktU, kktV;		// KKT along boundaries
			kktU = false;
			kktV = false;

			if (u == uDomain.beg()) {
				if (v == vDomain.beg()) {
					if (du >= 0 && dv >= 0) {
						return;
					}
					else if (du >= 0 && dv < 0)
						kktV = true;
					else if (du < 0 && dv >= 0)
						kktU = true;
				}
				else if (v == vDomain.end()) {
					if (du >= 0 && dv <= 0) {
						return;
					}
					else if (du >= 0 && dv > 0)
						kktV = true;
					else if (du < 0 && dv <= 0)
						kktU = true;
				}
				else {
					if (du >= 0)
						kktV = true;
				}
			}
			else if (u == uDomain.end()) {
				if (v == vDomain.beg()) {
					if (du <= 0 && dv >= 0) {
						return;
					}
					else if (du <= 0 && dv < 0)
						kktV = true;
					else if (du > 0 && dv >= 0)
						kktU = true;
				}
				else if (v == vDomain.end()) {
					if (du <= 0 && dv <= 0) {
						return;
					}
					else if (du <= 0 && dv > 0)
						kktV = true;
					else if (du > 0 && dv <= 0)
						kktU = true;
				}
				else {
					if (du <= 0)
						kktV = true;
				}
			}
			else {
				if (v == vDomain.beg()) {
					if (dv >= 0)
						kktU = true;
				}
				else if (v == vDomain.end()) {
					if (dv <= 0)
						kktU = true;
				}
			}

			if (kktU) {
				Real X = du * D;
				Real Xu = Tuu.dot(diff) + Tu.dot(Tu);
				if (Xu != 0) {
					Real deltaU = -X / Xu;

					if (deltaU * du < 0) {
						u = pu + deltaU;
						u = piDomain::regularize(u);

						if (!uDomain.has(u)) {
							if (deltaU > 0)
								u = piDomain::regularize(uDomain.end());
							else
								u = piDomain::regularize(uDomain.beg());
						}

						T = evaluate(u, v);
						diff = T - pt;
						Real tmpD = sqrt(diff.dot(diff));
						if (tmpD <= D) {
							if (tmpD > D - eps) {
								return;
							}
							D = tmpD;
							continue;
						}
					}
				}
			}
			else if (kktV) {
				Real X = dv * D;
				Real Xv = Tvv.dot(diff) + Tv.dot(Tv);
				if (Xv != 0) {
					Real deltaV = -X / Xv;

					if (deltaV * dv < 0) {
						v = pv + deltaV;
						v = piDomain::regularize(v);

						if (!vDomain.has(v)) {
							if (deltaV > 0)
								v = piDomain::regularize(vDomain.end());
							else
								v = piDomain::regularize(vDomain.beg());
						}

						T = evaluate(u, v);
						diff = T - pt;
						Real tmpD = sqrt(diff.dot(diff));
						if (tmpD <= D) {
							if (tmpD > D - eps) {
								return;
							}
							D = tmpD;
							continue;
						}
					}
				}
			}
			else {
				Real dvec[3], djac[3][3];
				dvec[1] = du * D;
				dvec[2] = dv * D;

				djac[1][1] = Tuu.dot(diff) + Tu.dot(Tu);
				djac[1][2] = Tuv.dot(diff) + Tu.dot(Tv);
				djac[2][1] = djac[1][2];
				djac[2][2] = Tvv.dot(diff) + Tv.dot(Tv);

				Real p[3], d;
				int indx[3];
				for (int i = 1; i <= 2; i++)
					p[i] = -dvec[i];								// Right-hand side of linear equations.
				ludcmp(djac, indx, &d);								// Solve linear equations using LU decomposition.
				lubksb(djac, indx, p);
				u = pu + p[1];
				v = pv + p[2];
				u = piDomain::regularize(u);
				v = piDomain::regularize(v);

				if (!uDomain.has(u)) {
					if (p[1] > 0)
						u = piDomain::regularize(uDomain.end());
					else
						u = piDomain::regularize(uDomain.beg());
				}

				if (!vDomain.has(v)) {
					if (p[2] > 0)
						v = piDomain::regularize(vDomain.end());
					else
						v = piDomain::regularize(vDomain.beg());
				}

				// Test KKT result
				T = evaluate(u, v);
				diff = T - pt;
				Real tmpD = diff.len();
				if (tmpD <= D) {
					if (tmpD > D - eps)
						return;
					D = tmpD;
					continue;
				}
			}

			// GR
			Real maxAlpha = 1.0;

			Real startAlpha = maxAlpha * 0.1;
			while (true) {
				cnt++;
				Real tmpAlpha = startAlpha;
				Real deltaU = -tmpAlpha * du;
				Real deltaV = -tmpAlpha * dv;
				u = pu + deltaU;
				v = pv + deltaV;
				u = piDomain::regularize(u);
				v = piDomain::regularize(v);

				if (!uDomain.has(u)) {
					if (deltaU > 0)
						u = piDomain::regularize(uDomain.end());
					else
						u = piDomain::regularize(uDomain.beg());
				}

				if (!vDomain.has(v)) {
					if (deltaV > 0)
						v = piDomain::regularize(vDomain.end());
					else
						v = piDomain::regularize(vDomain.beg());
				}

				T = evaluate(u, v);
				diff = T - pt;
				Real tmpD = sqrt(diff.dot(diff));
				if (tmpD <= D) {
					if (tmpD > D - eps)
						return;
					D = tmpD;
					break;
				}
				startAlpha *= 0.1;
			}
		}
	}
}