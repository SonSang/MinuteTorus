/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#include "CircleIntersect.h"

namespace MN {
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
	// Insert [ t ] into [ intParam ] if [ circle(t) ] point is on the plane defined by [ planeNormal, p ]
	inline static void addIntersect(const Circle& circle, Real t, const Vec3& planeNormal, Real p, Real intParam[2], int& intNum) {
		const static Real eps = 1e-10;
		auto pt = circle.evaluate(t);
		if (fabs(pt.dot(planeNormal) - p) < eps) {
			if (intNum == 2)
				throw("This should not happen, may be numerical error");
			intParam[intNum++] = t;
		}
			
	}
	int intersect(const Circle& circle, const Vec3& planeNormal, const Vec3& planePoint, Real intParam[2], int& intNum) {
		auto pn = planeNormal;
		pn.normalize();
		Real p = pn.dot(planePoint);

		intNum = 0;
		if (p == 0 && fabs(pn[2]) == 1) {
			// [ circle ] is fully included in the given plane
			return 0;
		}

		if (pn[0] == 0) {
			Real sinu = p / (pn[1] * circle.radius);
			if (fabs(sinu) > 1) {
				if (fabs(sinu) < 1 + 1e-10) {
					if (sinu > 0) sinu = 1;
					else sinu = -1;
				}
				else
					return 1;
			}
			Real u = asin(sinu);
			addIntersect(circle, u, pn, p, intParam, intNum);
			if (fabs(sinu) != 1)
				addIntersect(circle, PI - u, pn, p, intParam, intNum);
		}
		else if (pn[1] == 0) {
			Real cosu = p / (pn[0] * circle.radius);
			if (fabs(cosu) > 1) {
				if (fabs(cosu) < 1 + 1e-10) {
					if (cosu > 0) cosu = 1;
					else cosu = -1;
				}
				else
					return 1;
			}
			Real u = acos(cosu);
			addIntersect(circle, u, pn, p, intParam, intNum);
			if (fabs(cosu) != 1)
				addIntersect(circle, PI20 - u, pn, p, intParam, intNum);
		}
		else {
			Real c[3];
			c[0] = SQ(circle.radius) * (SQ(pn[0]) + SQ(pn[1]));
			c[1] = -p * pn[0] * circle.radius;
			c[2] = SQ(p) - SQ(pn[1] * circle.radius);

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
				addIntersect(circle, u, pn, p, intParam, intNum);
				if (fabs(cosu[i]) != 1)
					addIntersect(circle, PI20 - u, pn, p, intParam, intNum);
			}
		}
		return 1;
	}
	int intersect(const CircularArc& arc, const Vec3& planeNormal, const Vec3& planePoint, Real intParam[2], int& intNum) {
		auto pn = planeNormal;
		pn.normalize();
		Real p = pn.dot(planePoint);

		intNum = 0;
		if (p == 0 && fabs(pn[2]) == 1) {
			// [ circle ] is fully included in the given plane
			return 0;
		}

		if (pn[0] == 0) {
			Real sinu = p / (pn[1] * arc.radius);
			if (fabs(sinu) > 1) {
				if (fabs(sinu) < 1 + 1e-10) {
					if (sinu > 0) sinu = 1;
					else sinu = -1;
				}
				else
					return 1;
			}
			Real u = asin(sinu);
			if(arc.domain.has(u))
				addIntersect(arc, u, pn, p, intParam, intNum);
			if (fabs(sinu) != 1 && arc.domain.has(PI - u))
				addIntersect(arc, PI - u, pn, p, intParam, intNum);
		}
		else if (pn[1] == 0) {
			Real cosu = p / (pn[0] * arc.radius);
			if (fabs(cosu) > 1) {
				if (fabs(cosu) < 1 + 1e-10) {
					if (cosu > 0) cosu = 1;
					else cosu = -1;
				}
				else
					return 1;
			}
			Real u = acos(cosu);
			if(arc.domain.has(u))
				addIntersect(arc, u, pn, p, intParam, intNum);
			if (fabs(cosu) != 1 && arc.domain.has(PI20 - u))
				addIntersect(arc, PI20 - u, pn, p, intParam, intNum);
		}
		else {
			Real c[3];
			c[0] = SQ(arc.radius) * (SQ(pn[0]) + SQ(pn[1]));
			c[1] = -p * pn[0] * arc.radius;
			c[2] = SQ(p) - SQ(pn[1] * arc.radius);

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
				if(arc.domain.has(u))
					addIntersect(arc, u, pn, p, intParam, intNum);
				if (fabs(cosu[i]) != 1 && arc.domain.has(PI20 - u))
					addIntersect(arc, PI20 - u, pn, p, intParam, intNum);
			}
		}
		return 1;
	}
}