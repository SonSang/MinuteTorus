/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_CYLINDER_H__
#define __MN_CYLINDER_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "../Circle/Circle.h"

namespace MN {
	class Cylinder {
	public:
		Real radius;
		Domain vDomain;

		// [rcosu, rsinu, v] : Outward normal
		inline Vec3 evaluate(Real u, Real v) const {
			return {
				radius * cos(u),
				radius * sin(u),
				v
			};
		}
		// Up to 4th differentiation is allowed
		inline Vec3 differentiate(Real u, Real v, int uOrder, int vOrder) const {
			Real
				tmp;
			if (uOrder == 0 && vOrder == 0)
				return evaluate(u, v);
			else if (uOrder == 1 && vOrder == 0) {
				return { -radius * sin(u), radius * cos(u), 0.0 };
			}
			else if (uOrder == 0 && vOrder == 1) {
				return { 0, 0, 1 };
			}
			else if (uOrder == 2 && vOrder == 0) {
				return { -radius * cos(u), -radius * sin(u), 0.0 };
			}
			else if (uOrder == 1 && vOrder == 1) {
				return { 0, 0, 0.0 };
			}
			else if (uOrder == 0 && vOrder == 2) {
				return { 0, 0, 0 };
			}
			else if (uOrder == 3 && vOrder == 0)
			{
				return { radius * sin(u), -radius * cos(u), 0 };
			}
			else if (uOrder == 2 && vOrder == 1)
			{
				return { 0, 0, 0 };
			}
			else if (uOrder == 1 && vOrder == 2)
			{
				return { 0, 0, 0 };
			}
			else if (uOrder == 0 && vOrder == 3)
			{
				return { 0, 0, 0 };
			}
			else
				throw(std::runtime_error("Invalid cylinder differentiation order"));
		}
		inline Vec3 normal(Real u, Real v) const {
			Vec3
				Su = differentiate(u, v, 1, 0),
				Sv = differentiate(u, v, 0, 1),
				N = Su.cross(Sv);
			N.normalize();
			return N;
		}
		
		inline Circle majorCircle() const noexcept {
			Circle mc;
			mc.radius = radius;
			return mc;
		}

		/* Parameter */
		// Find parameter [ u0 ] such that [ MC(u0) ] is the closest point on major circle [ MC(u) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this cylinder, there is no unique [ u ]
		// @ 1 = Unique [ u0 ]
		inline virtual int findMinDistParamU(const Vec3& pt, Real& u) const {
			auto mc = majorCircle();
			return mc.findMinDistParam(pt, u);
		}

		// Find parameter [ u0 ] such that [ MC(u0) ] is the farthest point on major circle [ MC(u) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this cylinder, there is no unique [ u ]
		// @ 1 = Unique [ u0 ]
		inline virtual int findMaxDistParamU(const Vec3& pt, Real& u) const {
			auto mc = majorCircle();
			return mc.findMaxDistParam(pt, u);
		}

		// Find parameter [ u0 & u1 ] such that [ MC(u0), MC(u1) ] are the closest & farthest point on major circle [ MC(u) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this cylinder, there is no unique [ u0, u1 ]
		// @ 1 = Unique [ u0, u1 ]
		inline virtual int findExtDistParamU(const Vec3& pt, Real& minU, Real& maxU) const {
			auto mc = majorCircle();
			return mc.findExtDistParam(pt, minU, maxU);
		}

		// Find parameter [ v0 ] such that [ T(u, v0) ] is the closest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = [ v ] falls to the boundary of this cylinder
		// @ 1 = [ v ] falls to the inner part of this cylinder
		inline virtual int findMinDistParamV(const Vec3& pt, Real u, Real& v) const {
			if (pt[2] > vDomain.end()) {
				v = vDomain.end();
				return 0;
			}	
			else if (pt[2] < vDomain.beg()) {
				v = vDomain.beg();
				return 0;
			}	
			else {
				v = pt[2];
				return 1;
			}
		}

		// Find parameter [ v0 ] such that [ T(u, v0) ] is the farthest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = [ v ] falls to the boundary of this cylinder
		// @ 1 = [ v ] falls to the inner part of this cylinder
		inline virtual int findMaxDistParamV(const Vec3& pt, Real u, Real& v) const {
			Real d0 = fabs(pt[2] - vDomain.beg());
			Real d1 = fabs(pt[2] - vDomain.end());
			if (d0 < d1)
				v = vDomain.beg();
			else
				v = vDomain.end();
			return 0;
		}

		// Find parameter [ v0, v1 ] such that [ T(u, v0) & T(u, v1) ] are the closest & farthest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = [ minV ] is on boundary, [ maxV ] is on boundary ( of domain )
		// @ 1 = [ minV ] is on boundary, [ maxV ] is on inner ( of domain )
		// @ 2 = [ minV ] is on inner, [ maxV ] is on boundary ( of domain )
		// @ 3 = [ minV ] is on inner, [ maxV ] is on inner ( of domain )
		inline virtual int findExtDistParamV(const Vec3& pt, Real u, Real& minV, Real& maxV) const {
			int result0 = findMinDistParamV(pt, u, minV);
			int result1 = findMaxDistParamV(pt, u, maxV);
			if (result0 == 0 && result1 == 0)
				return 0;
			else if (result0 == 0 && result1 == 1)
				return 1;
			else if (result0 == 1 && result1 == 0)
				return 2;
			else
				return 3;
		}

		// Find parameter [ u0, v0 ] such that [ T(u0, v0) ] is the closest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this cylinder, there is no unique [ u, v ]
		// @ 1 = Unique [ u0, v0 ]
		inline virtual int findMinDistParam(const Vec3& pt, Real& u, Real& v) const {
			int uresult, vresult;
			uresult = findMinDistParamU(pt, u);
			vresult = findMinDistParamV(pt, u, v);
			return uresult;
		}

		// Find parameter [ u0, v0 ] such that [ T(u0, v0) ] is the farthest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this cylinder, there is no unique [ u, v ]
		// @ 1 = Unique [ u0, v0 ]
		inline virtual int findMaxDistParam(const Vec3& pt, Real& u, Real& v) const {
			int uresult, vresult;
			uresult = findMaxDistParamU(pt, u);
			vresult = findMaxDistParamV(pt, u, v);
			return uresult;
		}

		// Find parameter [ u0, v0 ] & [ u1, v1] such that [ T(u0, v0) ] is the closest point 
		// and [ T(u1, v1) ] is the farthest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = There is no unique [ minU, minV ], [ maxU, maxV ]
		// @ 1 = Unique [ minU, minV ], [ maxU, maxV ]
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
				return 1;
			}
		}

		/* Point */
		// Find point [ fpt ] which is the closest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this cylinder, there is no unique [ fpt ]
		// @ 1 = Unique [ fpt ]
		inline virtual int findMinDistPoint(const Vec3& pt, Vec3& fpt) const {
			Real u, v;
			int result = findMinDistParam(pt, u, v);
			fpt = evaluate(u, v);
			return result;
		}

		// Find point [ fpt ] which is the farthest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this cylinder, there is no unique [ fpt ]
		// @ 1 = Unique [ fpt ]
		inline virtual int findMaxDistPoint(const Vec3& pt, Vec3& fpt) const {
			Real u, v;
			int result = findMaxDistParam(pt, u, v);
			fpt = evaluate(u, v);
			return result;
		}

		// Find point [ minfpt ], [ maxfpt ] which is the closest & farthest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = There is no unique [ minfpt ] and [ maxfpt ]
		// @ 1 = Unique [ fpt ]s
		inline virtual int findExtDistPoint(const Vec3& pt, Vec3& minfpt, Vec3& maxfpt) const {
			Real2 minp, maxp;
			int result = findExtDistParam(pt, minp, maxp);
			minfpt = evaluate(minp.first, minp.second);
			maxfpt = evaluate(maxp.first, maxp.second);
			return result;
		}
	};

	class CylinderPatch : public Cylinder {
	public:
		piDomain uDomain;

		inline CircularArc majorCircularArc() const noexcept {
			CircularArc mc;
			mc.radius = radius;
			mc.domain = uDomain;
			return mc;
		}

		/* Parameter */
		// Find parameter [ u0 ] such that [ MC(u0) ] is the closest point on major circular arc [ MC(u) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this cylinder, there is no unique [ u ]
		// @ 1 = [ u ] falls to the boundary of this major circular arc's domain, since the closest point on the full circle is out of domain
		// @ 2 = [ u ] falls to the inner part of this major circular arc's domain
		inline virtual int findMinDistParamU(const Vec3& pt, Real& u) const {
			auto mc = majorCircularArc();
			return mc.findMinDistParam(pt, u);
		}

		// Find parameter [ u0 ] such that [ MC(u0) ] is the farthest point on major circular arc [ MC(u) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this cylinder, there is no unique [ u ]
		// @ 1 = [ u ] falls to the boundary of this major circular arc's domain, since the farthest point on the full circle is out of domain
		// @ 2 = [ u ] falls to the inner part of this major circular arc's domain
		inline virtual int findMaxDistParamU(const Vec3& pt, Real& u) const {
			auto mc = majorCircularArc();
			return mc.findMaxDistParam(pt, u);
		}

		// Find parameter [ u0 & u1 ] such that [ MC(u0) & MC(u1) ] is the closest & farthest point on major circular arc [ MC(u) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this cylinder, there is no unique [ minU, maxU ]
		// @ 1 = [ minU ] is on boundary, [ maxU ] is on boundary ( of domain )
		// @ 2 = [ minU ] is on boundary, [ maxU ] is on inner ( of domain )
		// @ 3 = [ minU ] is on inner, [ maxU ] is on boundary ( of domain )
		// @ 4 = [ minU ] is on inner, [ maxU ] is on inner ( of domain )
		inline virtual int findExtDistParamU(const Vec3& pt, Real& minU, Real& maxU) const {
			auto mc = majorCircularArc();
			return mc.findExtDistParam(pt, minU, maxU);
		}

		// Find parameter [ v0 ] such that [ T(u, v0) ] is the closest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = [ v ] falls to the boundary of this cylinder
		// @ 1 = [ v ] falls to the inner part of this cylinder
		inline virtual int findMinDistParamV(const Vec3& pt, Real u, Real& v) const {
			return Cylinder::findMinDistParamV(pt, u, v);
		}

		// Find parameter [ v0 ] such that [ T(u, v0) ] is the farthest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = [ v ] falls to the boundary of this cylinder
		// @ 1 = [ v ] falls to the inner part of this cylinder
		inline virtual int findMaxDistParamV(const Vec3& pt, Real u, Real& v) const {
			return Cylinder::findMaxDistParamV(pt, u, v);
		}

		// Find parameter [ v0 & v1 ] such that [ T(u, v0) & T(u, v1) ] is the closest & farthest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the major circle of this torus, there is no unique [ minV, maxV ]
		// @ 1 = [ minV ] is on boundary, [ maxV ] is on boundary ( of domain )
		// @ 2 = [ minV ] is on boundary, [ maxV ] is on inner ( of domain )
		// @ 3 = [ minV ] is on inner, [ maxV ] is on boundary ( of domain )
		// @ 4 = [ minV ] is on inner, [ maxV ] is on inner ( of domain )
		inline virtual int findExtDistParamV(const Vec3& pt, Real u, Real& minV, Real& maxV) const {
			return Cylinder::findExtDistParamV(pt, u, minV, maxV);
		}

		// Find parameter [ u0, v0 ] such that [ T(u0, v0) ] is the closest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this cylinder, there is no unique [ u, v ]
		// @ 1 = Unique [ u0, v0 ]
		inline virtual int findMinDistParam(const Vec3& pt, Real& u, Real& v) const {
			int uresult, vresult;
			uresult = findMinDistParamU(pt, u);
			vresult = findMinDistParamV(pt, u, v);
			return uresult;
		}

		// Find parameter [ u0, v0 ] such that [ T(u0, v0) ] is the farthest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this cylinder, there is no unique [ u, v ]
		// @ 1 = Unique [ u0, v0 ]
		inline virtual int findMaxDistParam(const Vec3& pt, Real& u, Real& v) const {
			int uresult, vresult;
			uresult = findMaxDistParamU(pt, u);
			vresult = findMaxDistParamV(pt, u, v);
			return uresult;
		}

		// Find parameter [ u0, v0 ] & [ u1, v1] such that [ T(u0, v0) ] is the closest point 
		// and [ T(u1, v1) ] is the farthest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = There is no unique [ minU, minV ], [ maxU, maxV ]
		// @ 1 = Unique [ minU, minV ], [ maxU, maxV ]
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
				return 1;
			}
		}

		/* Point */
		// Find point [ fpt ] which is the closest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this cylinder, there is no unique [ fpt ]
		// @ 1 = Unique [ fpt ]
		inline virtual int findMinDistPoint(const Vec3& pt, Vec3& fpt) const {
			Real u, v;
			int result = findMinDistParam(pt, u, v);
			fpt = evaluate(u, v);
			return result;
		}

		// Find point [ fpt ] which is the farthest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = Since [ pt ] is on the axis of this cylinder, there is no unique [ fpt ]
		// @ 1 = Unique [ fpt ]
		inline virtual int findMaxDistPoint(const Vec3& pt, Vec3& fpt) const {
			Real u, v;
			int result = findMaxDistParam(pt, u, v);
			fpt = evaluate(u, v);
			return result;
		}

		// Find point [ minfpt ], [ maxfpt ] which is the closest & farthest point on cylinder [ T(u, v) ] to the given point [ pt ]
		// @return : 
		// @ 0 = There is no unique [ minfpt ] and [ maxfpt ]
		// @ 1 = Unique [ fpt ]s
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