#ifndef __MN_TORUS_BINORMAL_H__
#define __MN_TORUS_BINORMAL_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "Torus.h"
#include "TorusGaussmap.h"
#include "../Circle/CircleBinormal.h"

namespace MN {
	class TorusBinormal {
	public:
		struct Binormal {
			int type;	// 0 : Point - Point binormal
						// 1 : (uDomainA, vA) - (uDomainB, vB) binormal ( Major circles share same center on XY plane, same axis, but different radius or center's Z value )
						// 2 : (uA, vDomainA) - (uB, vDomainB) binormal ( Same major cirlce point and parallel minor axis, but different minor radius )
						// 3 : (uDomainA, vDomainA) - (uDomainB, vDomainB) binormal ( Perfectly same major circle )
						// 4 : (uDomainA, vA) - (uB, vDomainB) binormal
						// 5 : (uA, vDomainA) - (uDomainB, vB) binormal
			Real uA, vA, uB, vB;
			Vec3 pointA;
			Vec3 pointB;
			Real length;			// These information is always valid

			//piDomain uDomainA;
			//piDomain vDomainA;
			//piDomain uDomainB;
			//piDomain vDomainB;		// These information depends on type
		};
		CircleBinormal circleBinormal;

		void solve(const Torus& a, const Torus& b, const Transform& ta, const Transform& tb, std::vector<Binormal>& bins);
		void fSolve(const Torus& a, const Torus& b, const Transform& atob, const Transform& btoa, std::vector<Binormal>& bins);

		void solve(const TorusPatch& a, const TorusPatch& b, const Transform& ta, const Transform& tb, std::vector<Binormal>& bins);
		void fSolve(const TorusPatch& a, const TorusPatch& b, const Transform& atob, const Transform& btoa, std::vector<Binormal>& bins);

		/*
		 * Find torus binormal with gaussmap information
		 *  @aOption, bOption : option[0] = Use outward normals, option[1] = Use inward normals
		 */
		void solve(const TPatchGmap& a, const TPatchGmap& b, const Transform& ta, const Transform& tb, bool aOutward, bool aInward, bool bOutward, bool bInward, std::vector<Binormal>& bins);
		void fSolve(const TPatchGmap& a, const TPatchGmap& b, const Transform& atob, const Transform& btoa, bool aOutward, bool aInward, bool bOutward, bool bInward, std::vector<Binormal>& bins);

		//void solve(const TorusPatch& a, const TorusPatch& b, const Transform& atob, const Transform& btoa, std::vector<Binormal>& bins);
		//void solve(const VTorus& a, const VTorus& b, const Transform& atob, const Transform& btoa, std::vector<Binormal>& bins);

		// Use gaussmap to accelerate process : Just use intersecting subset of [ b ]'s u domain between [ ga ] and [ gb ]
		// @gbOption :  0 if gaussmap parameter coincides with [ u ] parameter of [ gb ]
		//				1 if gaussmap parameter is opposite from [ u ] parameter of [ gb ]
		//				2 for safe routine
		//void solve(const VTorus& a, const VTorus& b, const Gaussmap& ga, const Gaussmap& gb, int gbOption, const Transform& atob, const Transform& btoa, std::vector<Binormal>& bins);
	};
}

#endif