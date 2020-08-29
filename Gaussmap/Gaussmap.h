/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_GAUSSMAP_H__
#define __MN_GAUSSMAP_H__

#include "../MinuteUtils/Utils.h"

namespace MN {
	/* 
	 *	Gauss map in spherical coordinates : (r * cosu * sinv, r * sinu * sinv, r * cosv)
	 */
	class Gaussmap {
	private:
		Gaussmap() = default;
	public:
		piDomain uDomain;
		piDomain vDomain;		// Width must be smaller than PI

		inline static void validDomainV(const piDomain& vDomain) {
			if (vDomain.beg() > PI || vDomain.end() < 0 || vDomain.width() > PI)
				throw("Invalid V domain for gaussmap");
		}
		
		inline static Gaussmap create(const piDomain& uDomain = piDomain::create(0, PI20), const piDomain& vDomain = piDomain::create(0, PI)) {
			validDomainV(vDomain);

			Gaussmap gm;
			gm.uDomain = uDomain;
			gm.vDomain = vDomain;
			return gm;
		}

		static inline Vec3 evaluate(Real u, Real v) noexcept {
			return { cos(u) * sin(v), sin(u) * sin(v), cos(v) };
		}
		static inline Vec3 normal(Real u, Real v) noexcept {
			return evaluate(u, v);
		}

		inline Gaussmap inverse() const noexcept {
			return create(piDomain::create(uDomain.beg() + PI, uDomain.end() + PI), piDomain::create(PI - vDomain.end(), PI - vDomain.beg()));
		}
	};

	// Legacy
	/*class Gaussmap {
	private:
		Gaussmap() = default;
	public:
		//// V - Boundary circle of gaussmap :
		//// Gaussmap has two boundary circles at its min & max V values
		//class bCircle {
		//public:
		//	Vec3 C;		// Center
		//	Vec3 U;		// X
		//	Vec3 V;		// Y
		//};
	private:
		// Compare this gaussmap against [ circle ]
		// Assume this gaussmap's local coordinates as standard one
		// @~Z : Z value of upper and lower bounding circle
		// @domain : Domain of [ circle ] that intersects with this gaussmap
		//void intersection(const Gcircle& circle, Real upperZ, Real lowerZ, Domain domain[2], int& domainNum) const;
	public:
		piDomain uDomain;
		piDomain vDomain;		// Width must be smaller than PI

		inline static void validDomainV(const piDomain& vDomain) {
			if (vDomain.begin > PI || vDomain.end < 0 || vDomain.width() > PI)
				throw("Invalid V domain for gaussmap");
		}

		inline static Gaussmap create(const piDomain& uDomain = piDomain::create(0, PI20), const piDomain& vDomain = piDomain::create(0, PI)) {
			validDomainV(vDomain);

			Gaussmap gm;
			gm.uDomain = uDomain;
			gm.vDomain = vDomain;
			return gm;
		}

		inline Gaussmap inverse() const noexcept {
			return create(piDomain::create(uDomain.begin() + PI, uDomain.end() + PI), piDomain::create(PI - vDomain.end(), PI - vDomain.begin()));
		}

		//// Compare this gaussmap against [ other ] gaussmap and find intersection domain of [ other ] gaussmap
		//// @transform : Transform that takes [ other ] gaussmap to this local coordinates
		//void intersection(const Gaussmap& other, const Transform& transform, Domain domain[3], int& domainNum) const;

		//// Compare two gaussmap [ a ] , [ b ] and find domain of each gaussmap that intersect with each other
		//static void intersection(
		//	const Gaussmap& a, const Gaussmap& b, const Transform& atob, const Transform& btoa,
		//	Domain domainA[3], Domain domainB[3], int& aNum, int& bNum);
	};*/
}

#endif