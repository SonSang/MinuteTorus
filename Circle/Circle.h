#ifndef __MN_CIRCLE_H__
#define __MN_CIRCLE_H__

#ifdef _MSC_VER
#pragma once
#endif

#include <vector>
#include "../MinuteUtils/Utils.h"

namespace MN {
	// Circle on XY plane
	class Circle {
	public:
		Real radius = (Real)0.0;

		// [ rcos(t), rsin(t), 0 ]
		inline Vec3 evaluate(Real t) const noexcept {
			return {
				radius * cos(t),
				radius * sin(t),
				0.0
			};
		}

		// @order : Differentiation order
		inline Vec3 differentiate(Real t, int order) const noexcept {
			int det = order % 4;
			if (det == 0)
				return evaluate(t);
			else if (det == 1)
				return { radius * -sin(t), radius * cos(t), 0 };
			else if (det == 2)
				return { radius * -cos(t), radius * -sin(t), 0 };
			else
				return { radius * sin(t), radius * -cos(t), 0 };
		}

		// Find parameter [ t0 ] such that [ C(t0) ] is the closest point on circle [ C(t) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this circle, there is no unique [ param ]
		// @ 1 = Unique [ param ]
		inline virtual int findMinDistParam(const Vec3& pt, Real& param) const {
			Real x = pt[0], y = pt[1], len;
			if (y == (Real)0.0) {
				if (x == (Real)0.0) {
					param = (Real)0.0;
					return 0;
				}
				else {
					if (x > 0) param = (Real)0.0;
					else if (x < 0) param = PI;
					return 1;
				}
			}
			else {
				len = sqrt(SQ(x) + SQ(y));
				x /= len;
				if (y > 0.0) param = acos(x);
				else param = (PI20 - acos(x));
				return 1;
			}
		}
		// Find parameter [ t0 ] such that [ C(t0) ] is the farthest point on circle [ C(t) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this circle, there is no unique [ param ]
		// @ 1 = Unique [ param ]
		inline virtual int findMaxDistParam(const Vec3& pt, Real& param) const {
			if (Circle::findMinDistParam(pt, param)) {
				param += PI;
				return 1;
			}
			else
				return 0;
		}
		// Find parameter [ t0, t1 ] such that [ C(t0), C(t1) ] are the closest & farthest points on circle [ C(t) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this circle, there is no unique [ param ]
		// @ 1 = Unique [ param ]
		inline virtual int findExtDistParam(const Vec3& pt, Real& minParam, Real& maxParam) const {
			auto ret = Circle::findMinDistParam(pt, minParam);
			maxParam = piDomain::regularize(minParam + PI);
			return ret;
		}

		// Find point [ fpt ] on this circle that is the closest point to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this circle, there is no unique [ param ]
		// @ 1 = Unique [ param ]
		inline virtual int findMinDistPoint(const Vec3& pt, Vec3& fpt) const {
			Real param;
			auto ret = Circle::findMinDistParam(pt, param);
			fpt = evaluate(param);
			return ret;
		}
		// Find point [ fpt ] on this circle that is the farthest point to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this circle, there is no unique [ param ]
		// @ 1 = Unique [ param ]
		inline virtual int findMaxDistPoint(const Vec3& pt, Vec3& fpt) const {
			Real param;
			auto ret = Circle::findMaxDistParam(pt, param);
			fpt = evaluate(param);
			return ret;
		}
		// Find point [ fpt ] on this circle that is the closest & farthest point to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this circle, there is no unique [ param ]
		// @ 1 = Unique [ param ]
		inline virtual int findExtDistPoint(const Vec3& pt, Vec3& minFpt, Vec3& maxFpt) const {
			Real minp, maxp;
			auto ret = Circle::findExtDistParam(pt, minp, maxp);
			minFpt = evaluate(minp);
			maxFpt = evaluate(maxp);
			return ret;
		}
	};

	// Circular Arc on XY plane
	class CircularArc : public Circle {
	public:
		piDomain domain;	// Needed to define "arc", which is subset of circle

		// Find parameter [ t0 ] such that [ C(t0) ] is the closest point on circlular arc [ C(t) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this circle, there is no unique [ param ]
		// @ 1 = [ param ] falls to the boundary of this circular arc's domain, since the closest point on the full circle is out of domain
		// @ 2 = [ param ] falls to the inner part of this circular arc's domain
		inline virtual int findMinDistParam(const Vec3& pt, Real& param) const {
			if (Circle::findMinDistParam(pt, param)) {
				if (!domain.has(param)) {
					// [ param ] is out of domain
					auto begDist = pt.distsq(evaluate(domain.beg()));
					auto endDist = pt.distsq(evaluate(domain.end()));
					param = (begDist < endDist) ? domain.beg() : domain.end();
					return 1;
				}
				else
					return 2;
			}
			else {
				param = domain.middle();
				return 0;
			}
		}
		// Find parameter [ t0 ] such that [ C(t0) ] is the farthest point on circlular arc [ C(t) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this circle, there is no unique [ param ]
		// @ 1 = [ param ] falls to the boundary of this circular arc's domain, since the farthest point on the full circle is out of domain
		// @ 2 = [ param ] falls to the inner part of this circular arc's domain
		inline virtual int findMaxDistParam(const Vec3& pt, Real& param) const {
			if (Circle::findMaxDistParam(pt, param)) {
				if (!domain.has(param)) {
					// [ param ] is out of domain
					auto begDist = pt.distsq(evaluate(domain.beg()));
					auto endDist = pt.distsq(evaluate(domain.end()));
					param = (begDist > endDist) ? domain.beg() : domain.end();
					return 1;
				}
				else
					return 2;
			}
			else {
				param = domain.middle();
				return 0;
			}
		}
		// Find parameter [ t0 ] such that [ C(t0) ] is the closest & farthest point on circlular arc [ C(t) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this circle, there is no unique [ param ]
		// @ 1 = [ minParam ] is on boundary, [ maxParam ] is on boundary ( of domain )
		// @ 2 = [ minParam ] is on boundary, [ maxParam ] is on inner ( of domain )
		// @ 3 = [ minParam ] is on inner, [ maxParam ] is on boundary ( of domain )
		// @ 4 = [ minParam ] is on inner, [ maxParam ] is on inner ( of domain )
		inline virtual int findExtDistParam(const Vec3& pt, Real& minParam, Real& maxParam) const {
			if (Circle::findExtDistParam(pt, minParam, maxParam)) {
				bool hasMin = domain.has(minParam);
				bool hasMax = domain.has(maxParam);
				if (hasMin && hasMax)
					return 4;
				else if (hasMin) {
					// [ maxParam ] is out of domain
					auto begDist = pt.distsq(evaluate(domain.beg()));
					auto endDist = pt.distsq(evaluate(domain.end()));
					if (begDist < endDist)
						maxParam = piDomain::regularize(domain.end());
					else
						maxParam = piDomain::regularize(domain.beg());
					return 3;
				}
				else if (hasMax) {
					// [ maxParam ] is out of domain
					auto begDist = pt.distsq(evaluate(domain.beg()));
					auto endDist = pt.distsq(evaluate(domain.end()));
					if (begDist < endDist)
						minParam = piDomain::regularize(domain.beg());
					else
						minParam = piDomain::regularize(domain.end());
					return 2;
				}
				else {
					// [ minParam ], [ maxParam ] are out of domain
					auto begDist = pt.distsq(evaluate(domain.beg()));
					auto endDist = pt.distsq(evaluate(domain.end()));
					if (begDist < endDist) {
						minParam = piDomain::regularize(domain.beg());
						maxParam = piDomain::regularize(domain.end());
					}
					else {
						minParam = piDomain::regularize(domain.end());
						maxParam = piDomain::regularize(domain.beg());
					}
					return 1;
				}
			}
			else {
				minParam = maxParam = piDomain::regularize(domain.middle());
				return 0;
			}
		}
	};
}

#endif