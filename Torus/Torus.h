/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_TORUS_H__
#define __MN_TORUS_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "../Circle/Circle.h"

namespace MN {
	class Torus {
	public:
		Real majorRadius;
		Real minorRadius;
		
		// [(R + rcosv)cosu, (R + rcosv)sinu, rsinv] : Outward normal
		inline Vec3 evaluate(Real u, Real v) const {
			Real tmp = majorRadius + minorRadius * cos(v);
			return {
				tmp * cos(u),
				tmp * sin(u),
				minorRadius * sin(v)
			};
		}
		// Up to 3rd differentiation is allowed
		inline Vec3 differentiate(Real u, Real v, int uOrder, int vOrder) const {
			Real
				tmp;
			if (uOrder == 0 && vOrder == 0)
				return evaluate(u, v);
			else if (uOrder == 1 && vOrder == 0) {
				tmp = majorRadius + minorRadius * cos(v);
				return { -tmp * sin(u), tmp * cos(u), 0.0 };
			}
			else if (uOrder == 0 && vOrder == 1) {
				tmp = -minorRadius * sin(v);
				return { tmp * cos(u), tmp * sin(u), minorRadius * cos(v) };
			}
			else if (uOrder == 2 && vOrder == 0) {
				tmp = majorRadius + minorRadius * cos(v);
				return { -tmp * cos(u), -tmp * sin(u), 0.0 };
			}
			else if (uOrder == 1 && vOrder == 1) {
				tmp = minorRadius * sin(v);
				return { tmp * sin(u), -tmp * cos(u), 0.0 };
			}
			else if (uOrder == 0 && vOrder == 2) {
				tmp = -minorRadius * cos(v);
				return { tmp * cos(u), tmp * sin(u), -minorRadius * sin(v) };
			}
			else 
				throw(std::runtime_error("Invalid torus differentiation order"));
		}
		inline Vec3 normal(Real u, Real v) const {
			Vec3
				Su = differentiate(u, v, 1, 0),
				Sv = differentiate(u, v, 0, 1),
				N = Su.cross(Sv);
			N.normalize();
			return N;
		}

		// For given u, give transform that takes torus coordinates to local coordiantes of [ u ] circular arc
		inline Transform uTransform(Real u) const {
			Transform transform;
			Vec3 pt;
			pt[0] = cos(u);
			pt[1] = sin(u);
			pt[2] = 0.0;

			transform.R[0][0] = pt[0]; transform.R[0][1] = pt[1]; transform.R[0][2] = 0;
			transform.R[1][0] = 0; transform.R[1][1] = 0; transform.R[1][2] = 1;
			transform.R[2][0] = -pt[1]; transform.R[2][1] = pt[0]; transform.R[2][2] = 0;

			pt[0] *= -majorRadius;
			pt[1] *= -majorRadius;
			transform.T = transform.R * pt;

			return transform;
		}

		inline Circle majorCircle() const noexcept {
			Circle mc;
			mc.radius = majorRadius;
			return mc;
		}
		inline Circle minorCircle() const noexcept {
			Circle mc;
			mc.radius = minorRadius;
			return mc;
		}

		/* Parameter */
		// Find parameter [ u0 ] such that [ MC(u0) ] is the closest point on major circle [ MC(u) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this torus, there is no unique [ u ]
		// @ 1 = Unique [ u0 ]
		inline virtual int findMinDistParamU(const Vec3& pt, Real& u) const {
			auto mc = majorCircle();
			return mc.findMinDistParam(pt, u);
		}

		// Find parameter [ u0 ] such that [ MC(u0) ] is the farthest point on major circle [ MC(u) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this torus, there is no unique [ u ]
		// @ 1 = Unique [ u0 ]
		inline virtual int findMaxDistParamU(const Vec3& pt, Real& u) const {
			auto mc = majorCircle();
			return mc.findMaxDistParam(pt, u);
		}

		// Find parameter [ u0 & u1 ] such that [ MC(u0), MC(u1) ] are the closest & farthest point on major circle [ MC(u) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this torus, there is no unique [ u0, u1 ]
		// @ 1 = Unique [ u0, u1 ]
		inline virtual int findExtDistParamU(const Vec3& pt, Real& minU, Real& maxU) const {
			auto mc = majorCircle();
			return mc.findExtDistParam(pt, minU, maxU);
		}

		// Find parameter [ v0 ] such that [ T(u, v0) ] is the closest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the major circle of this torus, there is no unique [ v ]
		// @ 1 = Unique [ v0 ]
		inline virtual int findMinDistParamV(const Vec3& pt, Real u, Real& v) const {
			auto mc = minorCircle();
			auto utransform = uTransform(u);
			Vec3 tmppt = utransform.apply(pt);
			return mc.findMinDistParam(pt, v);
		}

		// Find parameter [ v0 ] such that [ T(u, v0) ] is the farthest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the major circle of this torus, there is no unique [ v ]
		// @ 1 = Unique [ v0 ]
		inline virtual int findMaxDistParamV(const Vec3& pt, Real u, Real& v) const {
			auto mc = minorCircle();
			auto utransform = uTransform(u);
			Vec3 tmppt = utransform.apply(pt);
			return mc.findMaxDistParam(pt, v);
		}

		// Find parameter [ v0, v1 ] such that [ T(u, v0) & T(u, v1) ] are the closest & farthest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the major circle of this torus, there is no unique [ v0, v1 ]
		// @ 1 = Unique [ v0 ]
		inline virtual int findExtDistParamV(const Vec3& pt, Real u, Real& minV, Real& maxV) const {
			auto mc = minorCircle();
			auto utransform = uTransform(u);
			Vec3 tmppt = utransform.apply(pt);
			return mc.findExtDistParam(pt, minV, maxV);
		}

		// Find parameter [ u0, v0 ] such that [ T(u0, v0) ] is the closest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this torus and major circle, there is no unique [ u, v ]
		// @ 1 = Since [ pt ] is on the axis of this torus, there is no unique [ u ]
		// @ 2 = Since [ pt ] is on the major circle, there is no unique [ v ]
		// @ 3 = Unique [ u0, v0 ]
		inline virtual int findMinDistParam(const Vec3& pt, Real& u, Real& v) const {
			int uresult, vresult;
			uresult = findMinDistParamU(pt, u);
			if (uresult == 0) {
				u = 0;
				vresult = findMinDistParamV(pt, u, v);
				if (vresult == 0) {
					v = 0;
					return 0;
				}
				return 1;
			}
			else {
				vresult = findMinDistParamV(pt, u, v);
				if (vresult == 0) {
					v = 0;
					return 2;
				}
				return 3;
			}
		}

		// Find parameter [ u0, v0 ] such that [ T(u0, v0) ] is the farthest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this torus and major circle, there is no unique [ u, v ]
		// @ 1 = Since [ pt ] is on the axis of this torus, there is no unique [ u ]
		// @ 2 = Since [ pt ] is on the major circle, there is no unique [ v ]
		// @ 3 = Unique [ u0, v0 ]
		inline virtual int findMaxDistParam(const Vec3& pt, Real& u, Real& v) const {
			int uresult, vresult;
			uresult = findMaxDistParamU(pt, u);
			if (uresult == 0) {
				u = 0;
				vresult = findMaxDistParamV(pt, u, v);
				if (vresult == 0) {
					v = 0;
					return 0;
				}
				return 1;
			}
			else {
				vresult = findMaxDistParamV(pt, u, v);
				if (vresult == 0) {
					v = 0;
					return 2;
				}
				return 3;
			}
		}

		// Find parameter [ u0, v0 ] & [ u1, v1] such that [ T(u0, v0) ] is the closest point 
		// and [ T(u1, v1) ] is the farthest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = There is no unique [ minU, minV ], [ maxU, maxV ]
		// @ 1 = There is no unique [ minU, minV ]
		// @ 2 = There is no unique [ maxU, maxV ]
		// @ 3 = Unique [ minU, minV ], [ maxU, maxV ]
		inline virtual int findExtDistParam(const Vec3& pt, Real2& minParam, Real2& maxParam) const {
			int uresult, vresult0, vresult1;
			uresult = findExtDistParamU(pt, minParam.first, maxParam.first);
			if (uresult == 0) {
				minParam.first = 0;
				maxParam.first = PI;
				findExtDistParamV(pt, minParam.first, minParam.second, maxParam.second);
				return 0;
			}
			else {
				vresult0 = findMinDistParamV(pt, minParam.first, minParam.second);
				vresult1 = findMaxDistParamV(pt, maxParam.first, maxParam.second);
				if (vresult0 == 0 && vresult1 == 0)
					throw("It should not happen");
				else if (vresult0 == 0)
					return 1;
				else if (vresult1 == 0)
					return 2;
				else
					return 3;
			}
		}

		/* Point */
		// Find point [ fpt ] which is the closest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this torus or major circle, there is no unique [ fpt ]
		// @ 1 = Unique [ fpt ]
		inline virtual int findMinDistPoint(const Vec3& pt, Vec3& fpt) const {
			Real u, v;
			int result = findMinDistParam(pt, u, v);
			fpt = evaluate(u, v);
			return (result == 3);
		}

		// Find point [ fpt ] which is the farthest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this torus or major circle, there is no unique [ fpt ]
		// @ 1 = Unique [ fpt ]
		inline virtual int findMaxDistPoint(const Vec3& pt, Vec3& fpt) const {
			Real u, v;
			int result = findMaxDistParam(pt, u, v);
			fpt = evaluate(u, v);
			return (result == 3);
		}

		// Find point [ minfpt ], [ maxfpt ] which is the closest & farthest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = There is no unique [ minfpt ] and [ maxfpt ]
		// @ 1 = There is no unique [ minfpt ]
		// @ 2 = There is no unique [ maxfpt ]
		// @ 3 = Unique [ fpt ]s
		inline virtual int findExtDistPoint(const Vec3& pt, Vec3& minfpt, Vec3& maxfpt) const {
			Real2 minp, maxp;
			int result = findExtDistParam(pt, minp, maxp);
			minfpt = evaluate(minp.first, minp.second);
			maxfpt = evaluate(maxp.first, maxp.second);
			return result;
		}
		
		/*
		// For given [ pt ], project it on [ u ] or [ v ] parameter ( for only minimum distance )
		// @ret : 0 for too many projections case(center), 1 for success
		virtual inline int projectOnU(const Vec3& pt, Real& u) const {
			if (pt[0] == 0.0 && pt[1] == 0.0) {
				u = 0;
				return 0;
			}
			Real len = sqrt(pt[0] * pt[0] + pt[1] * pt[1]);
			Real cosu = pt[0] / len;
			Real sinu = pt[1] / len;
			u = acos(cosu);
			if (sinu < 0)
				u = PI20 - u;
			return 1;
		}
		virtual inline int projectOnV(const Vec3& pt, Real u, Real& v) const {
			Transform uTR;
			uTransform(u, uTR);
			Vec3 tmp = uTR.apply(pt);

			if (tmp[0] == 0.0 && tmp[1] == 0.0) {
				v = 0;
				return 0;
			}
			Real len = sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1]);
			Real cosv = tmp[0] / len;
			Real sinv = tmp[1] / len;
			v = acos(cosv);
			if (sinv < 0)
				v = PI20 - v;
			return 1;
		}
		
		// For given [ pt ], project it on [ u ] or [ v ] parameter ( for both min & max distance )
		virtual int projectAllOnU(const Vec3& pt, Real u[2]) const;
		virtual int projectAllOnV(const Vec3& pt, Real v[2]) const;

		// Projection not only for min distance, but also for max distance
		void projectAllOnV(const Vec3& pt, Real u, Real v[2], bool valid[2]) const;
		
		virtual bool project(const Vec3& pt, Vec3& fpt) const;
		*/
		
	};

	class TorusPatch : public Torus {
	public:
		piDomain uDomain;
		piDomain vDomain;

		inline CircularArc majorCircularArc() const noexcept {
			CircularArc mc;
			mc.radius = majorRadius;
			mc.domain = uDomain;
			return mc;
		}
		inline CircularArc minorCircularArc() const noexcept {
			CircularArc mc;
			mc.radius = minorRadius;
			mc.domain = vDomain;
			return mc;
		}

		/* Parameter */
		// Find parameter [ u0 ] such that [ MC(u0) ] is the closest point on major circular arc [ MC(u) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this torus, there is no unique [ u ]
		// @ 1 = [ u ] falls to the boundary of this major circular arc's domain, since the closest point on the full circle is out of domain
		// @ 2 = [ u ] falls to the inner part of this major circular arc's domain
		inline virtual int findMinDistParamU(const Vec3& pt, Real& u) const {
			auto mc = majorCircularArc();
			return mc.findMinDistParam(pt, u);
		}

		// Find parameter [ u0 ] such that [ MC(u0) ] is the farthest point on major circular arc [ MC(u) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this torus, there is no unique [ u ]
		// @ 1 = [ u ] falls to the boundary of this major circular arc's domain, since the farthest point on the full circle is out of domain
		// @ 2 = [ u ] falls to the inner part of this major circular arc's domain
		inline virtual int findMaxDistParamU(const Vec3& pt, Real& u) const {
			auto mc = majorCircularArc();
			return mc.findMaxDistParam(pt, u);
		}

		// Find parameter [ u0 & u1 ] such that [ MC(u0) & MC(u1) ] is the closest & farthest point on major circular arc [ MC(u) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this torus, there is no unique [ minU, maxU ]
		// @ 1 = [ minU ] is on boundary, [ maxU ] is on boundary ( of domain )
		// @ 2 = [ minU ] is on boundary, [ maxU ] is on inner ( of domain )
		// @ 3 = [ minU ] is on inner, [ maxU ] is on boundary ( of domain )
		// @ 4 = [ minU ] is on inner, [ maxU ] is on inner ( of domain )
		inline virtual int findExtDistParamU(const Vec3& pt, Real& minU, Real& maxU) const {
			auto mc = majorCircularArc();
			return mc.findExtDistParam(pt, minU, maxU);
		}

		// Find parameter [ v0 ] such that [ T(u, v0) ] is the closest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the major circle of this torus, there is no unique [ v ]
		// @ 1 = [ v ] falls to the boundary of this minor circular arc's domain, since the closest point on the full circle is out of domain
		// @ 2 = [ v ] falls to the inner part of this minor circular arc's domain
		inline virtual int findMinDistParamV(const Vec3& pt, Real u, Real& v) const {
			auto mc = minorCircularArc();
			auto utransform = uTransform(u);
			Vec3 tmppt = utransform.apply(pt);
			return mc.findMinDistParam(tmppt, v);
		}

		// Find parameter [ v0 ] such that [ T(u, v0) ] is the farthest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the major circle of this torus, there is no unique [ v ]
		// @ 1 = [ v ] falls to the boundary of this minor circular arc's domain, since the closest point on the full circle is out of domain
		// @ 2 = [ v ] falls to the inner part of this minor circular arc's domain
		inline virtual int findMaxDistParamV(const Vec3& pt, Real u, Real& v) const {
			auto mc = minorCircularArc();
			auto utransform = uTransform(u);
			Vec3 tmppt = utransform.apply(pt);
			return mc.findMaxDistParam(pt, v);
		}

		// Find parameter [ v0 & v1 ] such that [ T(u, v0) & T(u, v1) ] is the closest & farthest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the major circle of this torus, there is no unique [ minV, maxV ]
		// @ 1 = [ minV ] is on boundary, [ maxV ] is on boundary ( of domain )
		// @ 2 = [ minV ] is on boundary, [ maxV ] is on inner ( of domain )
		// @ 3 = [ minV ] is on inner, [ maxV ] is on boundary ( of domain )
		// @ 4 = [ minV ] is on inner, [ maxV ] is on inner ( of domain )
		inline virtual int findExtDistParamV(const Vec3& pt, Real u, Real& minV, Real& maxV) const {
			auto mc = minorCircularArc();
			auto utransform = uTransform(u);
			Vec3 tmppt = utransform.apply(pt);
			return mc.findExtDistParam(pt, minV, maxV);
		}

		// Find parameter [ u0, v0 ] such that [ T(u0, v0) ] is the closest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this torus and major circle, there is no unique [ u, v ]
		// @ 1 = Since [ pt ] is on the axis of this torus, there is no unique [ u ]
		// @ 2 = Since [ pt ] is on the major circle, there is no unique [ v ]
		// @ 3 = Unique [ u0, v0 ]
		inline virtual int findMinDistParam(const Vec3& pt, Real& u, Real& v) const {
			int uresult, vresult;
			uresult = findMinDistParamU(pt, u);
			if (uresult == 0) {
				vresult = findMinDistParamV(pt, u, v);
				if (vresult == 0)
					return 0;
				return 1;
			}
			else {
				vresult = findMinDistParamV(pt, u, v);
				if (vresult == 0)
					return 2;
				return 3;
			}
		}

		// Find parameter [ u0, v0 ] such that [ T(u0, v0) ] is the farthest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this torus and major circle, there is no unique [ u, v ]
		// @ 1 = Since [ pt ] is on the axis of this torus, there is no unique [ u ]
		// @ 2 = Since [ pt ] is on the major circle, there is no unique [ v ]
		// @ 3 = Unique [ u0, v0 ]
		inline virtual int findMaxDistParam(const Vec3& pt, Real& u, Real& v) const {
			int uresult, vresult;
			uresult = findMaxDistParamU(pt, u);
			if (uresult == 0) {
				vresult = findMaxDistParamV(pt, u, v);
				if (vresult == 0)
					return 0;
				return 1;
			}
			else {
				vresult = findMaxDistParamV(pt, u, v);
				if (vresult == 0)
					return 2;
				return 3;
			}
		}

		// Find parameter [ u0, v0 ] & [ u1, v1] such that [ T(u0, v0) ] is the closest point 
		// and [ T(u1, v1) ] is the farthest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = There is no unique [ minU, minV ], [ maxU, maxV ]
		// @ 1 = There is no unique [ minU, minV ]
		// @ 2 = There is no unique [ maxU, maxV ]
		// @ 3 = Unique [ minU, minV ], [ maxU, maxV ]
		inline virtual int findExtDistParam(const Vec3& pt, Real2& minParam, Real2& maxParam) const {
			int uresult, vresult0, vresult1;
			uresult = findExtDistParamU(pt, minParam.first, maxParam.first);
			if (uresult == 0) {
				findExtDistParamV(pt, minParam.first, minParam.second, maxParam.second);
				return 0;
			}
			else {
				vresult0 = findMinDistParamV(pt, minParam.first, minParam.second);
				vresult1 = findMaxDistParamV(pt, maxParam.first, maxParam.second);
				if (vresult0 == 0 && vresult1 == 0)
					throw("It should not happen");
				else if (vresult0 == 0)
					return 1;
				else if (vresult1 == 0)
					return 2;
				else
					return 3;
			}
		}

		/* Point */
		// Find point [ fpt ] which is the closest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this torus or major circle, there is no unique [ fpt ]
		// @ 1 = Unique [ fpt ]
		inline virtual int findMinDistPoint(const Vec3& pt, Vec3& fpt) const {
			Real u, v;
			int result = findMinDistParam(pt, u, v);
			fpt = evaluate(u, v);
			return (result == 3);
		}

		// Find point [ fpt ] which is the farthest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this torus or major circle, there is no unique [ fpt ]
		// @ 1 = Unique [ fpt ]
		inline virtual int findMaxDistPoint(const Vec3& pt, Vec3& fpt) const {
			Real u, v;
			int result = findMaxDistParam(pt, u, v);
			fpt = evaluate(u, v);
			return (result == 3);
		}

		// Find point [ minfpt ], [ maxfpt ] which is the closest & farthest point on torus [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = There is no unique [ minfpt ] and [ maxfpt ]
		// @ 1 = There is no unique [ minfpt ]
		// @ 2 = There is no unique [ maxfpt ]
		// @ 3 = Unique [ fpt ]s
		inline virtual int findExtDistPoint(const Vec3& pt, Vec3& minfpt, Vec3& maxfpt) const {
			Real2 minp, maxp;
			int result = findExtDistParam(pt, minp, maxp);
			minfpt = evaluate(minp.first, minp.second);
			maxfpt = evaluate(maxp.first, maxp.second);
			return result;
		}
	};
}

#endif