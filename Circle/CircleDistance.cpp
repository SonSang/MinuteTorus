/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#include "CircleDistance.h"

#define ITMAX 100
#define EPS 1.0e-10
#define ZEPS 1.0e-10
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

namespace MN {
	struct Vranek {
		Real a[10];
		Real t, cost, sint, cost2, sint2, costsint, root;
	};

	inline static void updateT(Real t, Vranek& vv) {
		if (vv.t != t) {
			vv.t = t;
			vv.cost = cos(t);
			vv.sint = sin(t);
			vv.cost2 = SQ(vv.cost);
			vv.sint2 = SQ(vv.sint);
			vv.costsint = vv.cost * vv.sint;
			vv.root = sqrt(vv.a[5] * vv.cost2 + vv.a[4] * vv.sint2 + vv.a[3] * vv.cost + vv.a[2] * vv.sint + vv.a[1] * vv.costsint + vv.a[0]);
		}
	}

	// Squared distance function f(t) that is to be minimized in vranek's scheme
	// [t] refers to the radian angle of [circle_b]
	inline static Real vranekEvaluate(Real t, Vranek& vv) {
		updateT(t, vv);
		return vv.a[9] * vv.cost + vv.a[8] * vv.sint + vv.a[7] + vv.a[6] * vv.root;
	}

	// First derivative version of [vranekEvaluate].
	inline static Real vranekDifferentiate(Real t, Vranek& vv) {
		updateT(t, vv);
		return -vv.a[9] * vv.sint + vv.a[8] * vv.cost + vv.a[6] * 0.5 * (1.0 / vv.root) * (2.0 * vv.costsint * (vv.a[4] - vv.a[5]) - vv.a[3] * vv.sint + vv.a[2] * vv.cost + vv.a[1] * (vv.cost2 - vv.sint2));
	}

	inline static Real vranekZbrent(Real x1, Real x2, Real tol, Vranek& vv) {
		int iter;
		Real a = x1, b = x2, c = x2, d, e, min1, min2;
		Real(*func)(Real, Vranek&) = vranekDifferentiate;
		Real fa = (func)(a, vv), fb = (func)(b, vv), fc, p, q, r, s, tol1, xm;
		if ((fa > 0. && fb > 0.) || (fa < 0. && fb < 0.)) 
			throw(std::runtime_error("Root must be bracketed in zbrent"));
		fc = fb;
		for (iter = 1; iter <= ITMAX; iter++) {
			if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
				c = a;
				fc = fa;
				e = d = b - a;
			}
			if (fabs(fc) < fabs(fb)) {
				a = b;
				b = c;
				c = a;
				fa = fb;
				fb = fc;
				fc = fa;
			}
			tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
			xm = 0.5 * (c - b);
			if (fabs(xm) <= tol1 || fb == 0.0) return b;
			if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
				s = fb / fa;
				if (a == c) {
					p = 2.0 * xm * s;
					q = 1.0 - s;
				}
				else {
					q = fa / fc;
					r = fb / fc;
					p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
					q = (q - 1.0) * (r - 1.0) * (s - 1.0);
				}
				if (p > 0.0) q = -q;
				p = fabs(p);
				min1 = 3.0 * xm * q - fabs(tol1 * q);
				min2 = fabs(e * q);
				if (2.0 * p < (min1 < min2 ? min1 : min2)) {
					e = d;
					d = p / q;
				}
				else {
					d = xm;
					e = d;
				}
			}
			else {
				d = xm;
				e = d;
			}
			a = b;
			fa = fb;
			if (fabs(d) > tol1)
				b += d;
			else
				b += SIGN(tol1, xm);
			fb = (*func)(b, vv);
		}
		throw(std::runtime_error("Maximum number of iterations exceeded in zbrent"));
	}
	inline static Real vranekDbrent(Real ax, Real bx, Real cx, Vranek& vv, Real tol, Real* xmin) {
		int iter, ok1, ok2;
		Real a, b, d, d1, d2, du, dv, dw, dx, e = 0.0;
		Real fu, fv, fw, fx, olde, tol1, tol2, u, u1, u2, v, w, x, xm;
		Real(*f)(Real, Vranek&) = vranekEvaluate;
		Real(*df)(Real, Vranek&) = vranekDifferentiate;

		a = (ax < cx ? ax : cx);
		b = (ax > cx ? ax : cx);
		x = w = v = bx;
		fw = fv = fx = (*f)(x, vv);
		dw = dv = dx = (*df)(x, vv);
		for (iter = 1; iter <= ITMAX; iter++) {
			xm = 0.5 * (a + b);
			tol1 = tol * fabs(x) + ZEPS;
			tol2 = 2.0 * tol1;
			if (fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
				*xmin = x;
				return fx;
			}
			if (fabs(e) > tol1) {
				d1 = 2.0 * (b - a);
				d2 = d1;
				if (dw != dx) d1 = (w - x) * dx / (dx - dw);
				if (dv != dx) d2 = (v - x) * dx / (dx - dv);
				u1 = x + d1;
				u2 = x + d2;
				ok1 = (a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0;
				ok2 = (a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0;
				olde = e;
				e = d;
				if (ok1 || ok2) {
					if (ok1 && ok2)
						d = (fabs(d1) < fabs(d2) ? d1 : d2);
					else if (ok1)
						d = d1;
					else
						d = d2;
					if (fabs(d) <= fabs(0.5 * olde)) {
						u = x + d;
						if (u - a < tol2 || b - u < tol2)
							d = SIGN(tol1, xm - x);
					}
					else {
						d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
					}
				}
				else {
					d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
				}
			}
			else {
				d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
			}
			if (fabs(d) >= tol1) {
				u = x + d;
				fu = (*f)(u, vv);
			}
			else {
				u = x + SIGN(tol1, d);
				fu = (*f)(u, vv);
				if (fu > fx) {
					*xmin = x;
					return fx;
				}
			}
			du = (*df)(u, vv);
			if (fu <= fx) {
				if (u >= x) a = x; else b = x;
				MOV3(v, fv, dv, w, fw, dw)
					MOV3(w, fw, dw, x, fx, dx)
					MOV3(x, fx, dx, u, fu, du)
			}
			else {
				if (u < x) a = u; else b = u;
				if (fu <= fw || w == x) {
					MOV3(v, fv, dv, w, fw, dw)
						MOV3(w, fw, dw, u, fu, du)
				}
				else if (fu < fv || v == x || v == w) {
					MOV3(v, fv, dv, u, fu, du)
				}
			}
		}
		throw(std::runtime_error("Too many iterations in Vranek Dbrent"));
	}

	inline static void cosclamp(Real& cosval) {
		if (cosval > 1 && cosval < 1 + 1e-5)
			cosval = 1;
		else if (cosval < -1 && cosval > -1 - 1e-5)
			cosval = -1;
	}

	Distance circleDistance(const Circle& a, const Circle& b, const Transform& tA, const Transform& tB) {
		// @cen_a, cen_b	: center of circle [ca], [cb].
		// @norm_a, norm_b	: axis of circle [ca], [cb].
		// @U, V			: local X, Y axis of [cb].
		Vec3 cen_a{ tA.T }, cen_b{ tB.T },
			norm_a{ tA.R[0][2], tA.R[1][2], tA.R[2][2] }, norm_b{ tB.R[0][2], tB.R[1][2], tB.R[2][2] },
			U{ tB.R[0][0], tB.R[1][0], tB.R[2][0] }, V{ tB.R[0][1], tB.R[1][1], tB.R[2][1] };
		Transform btoa = Transform::connect(tB, tA);
		
		Vranek vv;
		// calculate a[10].
		{
			Real tmp[11];
			tmp[0] = cen_a.dot(cen_a);
			tmp[1] = cen_a.dot(cen_b);
			tmp[2] = cen_b.dot(cen_b);
			tmp[3] = norm_a.dot(cen_a);
			tmp[4] = norm_a.dot(cen_b);
			tmp[5] = norm_a.dot(U);
			tmp[6] = norm_a.dot(V);
			tmp[7] = cen_a.dot(U);
			tmp[8] = cen_a.dot(V);
			tmp[9] = cen_b.dot(U);
			tmp[10] = cen_b.dot(V);
			vv.a[0] = b.radius * b.radius + tmp[2] + tmp[0] - 2.0 * tmp[1] - pow(tmp[4] - tmp[3], 2.0);
			vv.a[1] = -2.0 * b.radius * b.radius * tmp[5] * tmp[6];
			vv.a[2] = 2.0 * b.radius * tmp[10] - 2.0 * b.radius * tmp[8] - 2.0 * b.radius * (tmp[4] - tmp[3]) * tmp[6];
			vv.a[3] = 2.0 * b.radius * tmp[9] - 2.0 * b.radius * tmp[7] - 2.0 * b.radius * (tmp[4] - tmp[3]) * tmp[5];
			vv.a[4] = -pow(b.radius * tmp[6], 2.0);
			vv.a[5] = -pow(b.radius * tmp[5], 2.0);
			vv.a[6] = -2.0 * a.radius;
			vv.a[7] = a.radius * a.radius + b.radius * b.radius + tmp[0] + tmp[2] - 2.0 * tmp[1];
			vv.a[8] = 2.0 * b.radius * tmp[10] - 2.0 * b.radius * tmp[8];
			vv.a[9] = 2.0 * b.radius * tmp[9] - 2.0 * b.radius * tmp[7];
		}

		// 1. set random triplet and find single local minimum.
		Real
			loc_min,
			loc_min_t;
		{
			Real
				a = (-2.0 / 3.0) * PI,
				b = 0,
				c = (2.0 / 3.0) * PI * 1.1;
			Real
				fa = vranekEvaluate(a, vv),
				fb = vranekEvaluate(b, vv),
				fc = vranekEvaluate(c, vv);
			if (fb < fa) {
				if (fb > fc) {
					Real
						tmp = a;
					a = b;
					b = c;
					c = tmp + PI20;
				}
				else if (fb == fc) 
					throw(std::runtime_error("Vranek error"));
			}
			else if (fb > fa) {
				if (fb > fc) {
					if (fa > fc) {
						Real
							tmp = a;
						a = b;
						b = c;
						c = tmp + PI20;
					}
					else if (fa < fc) {
						Real
							tmp = a;
						a = c - PI20;
						c = b;
						b = tmp;
					}
					else 
						throw(std::runtime_error("Vranek error"));
				}
				else if (fb < fc) {
					Real
						tmp = a;
					a = c - PI20;
					c = b;
					b = tmp;
				}
				else 
					throw(std::runtime_error("Vranek error"));
			}
			else 
				throw(std::runtime_error("Vranek error"));

			loc_min = vranekDbrent(a, b, c, vv, 1e-10, &loc_min_t);
			if (loc_min <0 && loc_min > -1e-10)
				loc_min = 0;
		}

		// 2. formulate g(t), which is 4th degree polynomial.
		Real
			coefB[5];	// coefficients of g(t).
		{
			Real
				d[3] = { vv.a[7] - loc_min, vv.a[6] * vv.a[6], vv.a[8] * vv.a[8] },
				c1 = 2.0 * vv.a[8] * d[0] - d[1] * vv.a[2],
				c2 = d[1] * vv.a[4],
				c3 = d[0] * d[0] - d[1] * vv.a[0] + d[2] - d[1] * vv.a[4],
				c4 = 2.0 * vv.a[9] * d[0] - d[1] * vv.a[3],
				c5 = 2.0 * vv.a[9] * vv.a[8] - d[1] * vv.a[1],
				c6 = -d[2] + d[1] * vv.a[4] + vv.a[9] * vv.a[9] - d[1] * vv.a[5];
			//b[0] = c3 * c3 - c1 * c1;
			//b[1] = -2.0 * c1 * c5 + 2.0 * c3 * c4;
			coefB[2] = c1 * c1 - c5 * c5 + 2.0 * c3 * c6 + c4 * c4;
			coefB[3] = 2.0 * c1 * c5 + 2.0 * c4 * c6;
			coefB[4] = c6 * c6 + c5 * c5;
		}

		// 3. reduce g(t) into 2nd degree polynomial by [loc_min_t].
		Real
			eq[3];
		{
			Real
				tmp[4];
			Real
				z = cos(loc_min_t),
				z2 = z * z,
				z3 = z2 * z;
			//tmp[0] = b[1] + b[2] * z + b[3] * z2 + b[4] * z3;
			tmp[1] = coefB[2] + coefB[3] * z + coefB[4] * z2;
			tmp[2] = coefB[3] + coefB[4] * z;
			tmp[3] = coefB[4];

			eq[0] = tmp[1] + tmp[2] * z + tmp[3] * z2;
			eq[1] = tmp[2] + tmp[3] * z;
			eq[2] = tmp[3];
		}

		// 4. find real global minimum.
		{
			Real
				b_4ac = eq[1] * eq[1] - 4.0 * eq[0] * eq[2];
			Real
				gl_min = loc_min,
				gl_min_t = loc_min_t;
			// 4-1. global minimum already found.
			/*if (b_4ac <= 1e-7)
			{
				gl_min = loc_min;
				gl_min_t = loc_min_t;
			}*/
			// 4-2. global minimum has to be found.
			if (b_4ac >= 0.0)
			{
				Real
					root_b_4ac = sqrt(b_4ac),
					cosa = (-eq[1] + root_b_4ac) / (2.0 * eq[2]),
					cosb = (-eq[1] - root_b_4ac) / (2.0 * eq[2]);
				cosclamp(cosa);
				cosclamp(cosb);

				/*tmp[0] = (SQ(vv.a[9]) - SQ(vv.a[8]) - vv.a[5] * SQ(vv.a[6]) + vv.a[4] * SQ(vv.a[6]));
				tmp[1] = (-2.0 * loc_min * vv.a[9] + 2.0 * vv.a[7] * vv.a[9] - vv.a[3] * SQ(vv.a[6]));
				tmp[2] = (SQ(loc_min) + SQ(vv.a[8]) + SQ(vv.a[7]) - 2.0 * loc_min * vv.a[7] - vv.a[4] * SQ(vv.a[6]) - vv.a[0] * SQ(vv.a[6]));
				tmp[3] = (vv.a[1] * SQ(vv.a[6]) - 2.0 * vv.a[8] * vv.a[9]);
				tmp[4] = vv.a[2] * SQ(vv.a[6]) + 2.0 * loc_min * vv.a[8] - 2.0 * vv.a[7] * vv.a[8];
				{
					denom = SQ(cosa) * tmp[0] + cosa * tmp[1] + tmp[2];
					nom = cosa * tmp[3] + tmp[4];
					if (denom * nom < 0)
						a = pi20 - a;
				}
				{
					denom = SQ(cosb) * tmp[0] + cosb * tmp[1] + tmp[2];
					nom = cosb * tmp[3] + tmp[4];
					if (denom * nom < 0)
						b = pi20 - b;
				}*/
				{
					Real a = acos(cosa), b = acos(cosb);
					if (vranekDifferentiate(a, vv) * vranekDifferentiate(b, vv) < 0.0)
					{
						loc_min_t = vranekZbrent(a, b, 1e-10, vv);
						loc_min = vranekEvaluate(loc_min_t, vv);
						if (loc_min < gl_min) {
							gl_min = loc_min;
							gl_min_t = loc_min_t;
						}
					}
				}
				{
					Real a = acos(cosa), b = PI20 - acos(cosb);
					if (vranekDifferentiate(a, vv) * vranekDifferentiate(b, vv) < 0.0)
					{
						loc_min_t = vranekZbrent(a, b, 1e-10, vv);
						loc_min = vranekEvaluate(loc_min_t, vv);
						if (loc_min < gl_min) {
							gl_min = loc_min;
							gl_min_t = loc_min_t;
						}
					}
				}
				{
					Real a = PI20 - acos(cosa), b = acos(cosb);
					if (vranekDifferentiate(a, vv) * vranekDifferentiate(b, vv) < 0.0)
					{
						loc_min_t = vranekZbrent(a, b, 1e-10, vv);
						loc_min = vranekEvaluate(loc_min_t, vv);
						if (loc_min < gl_min) {
							gl_min = loc_min;
							gl_min_t = loc_min_t;
						}
					}
				}
				{
					Real a = PI20 - acos(cosa), b = PI20 - acos(cosb);
					if (vranekDifferentiate(a, vv) * vranekDifferentiate(b, vv) < 0.0)
					{
						loc_min_t = vranekZbrent(a, b, 1e-10, vv);
						loc_min = vranekEvaluate(loc_min_t, vv);
						if (loc_min < gl_min) {
							gl_min = loc_min;
							gl_min_t = loc_min_t;
						}
					}
				}
			}
			Distance distance;
			distance.length = sqrt(gl_min);
			distance.paramB[0] = gl_min_t;
			distance.pointB = b.evaluate(gl_min_t);
			a.findMinDistParam(btoa.apply(distance.pointB), distance.paramA[0]);
			distance.pointA = a.evaluate(distance.paramA[0]);
			return distance;
		}
	}
}