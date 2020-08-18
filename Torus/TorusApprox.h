#ifndef __MN_TORUS_APPROX_H__
#define __MN_TORUS_APPROX_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "Torus.h"

namespace MN {
	// Torus patch to approximate certain surface [F]
	// This class contains mapping class for connecting two different parameter spaces : [F]'s & [G]'s.('G' is torus surface)
	class TorusApprox {
	public:
		using Ptr = std::shared_ptr<TorusApprox>;
		struct SurfaceInfo {
			Vec3 F;
			Vec3 Fu;
			Vec3 Fv;
			Vec3 Fuu;
			Vec3 Fuv;
			Vec3 Fvv;

			Real u;
			Real v;				// F(u, v) is the approximating point

			Real M1;
			Real M2;
			Real M3;
			Real M4;
			Domain uDomain;
			Domain vDomain;		// Deriv and domain information is needed to determine error bound
		};
		class Mapping {
		public:
			Real m0;
			Real n0;		// Parameter of the approximation point
			std::array<Real, 6> mCoefs;	// Coefficients for [ m ] : m(u, v) = c[0] * u^2 + c[1] * v^2 + c[2] * u * v + c[3] * u + c[4] * v + c[5]
			std::array<Real, 6> nCoefs;	// Coefficients for [ n ] : n(u, v) = c[0] * u^2 + c[1] * v^2 + c[2] * u * v + c[3] * u + c[4] * v + c[5]
			std::array<Real, 6> uCoefs;	// Coefficients for [ u ] : u(m, n) = c[0] * m^2 + c[1] * n^2 + c[2] * m * n + c[3] * m + c[4] * n + c[5]
			std::array<Real, 6> vCoefs;	// Coefficients for [ v ] : v(m, n) = c[0] * m^2 + c[1] * n^2 + c[2] * m * n + c[3] * m + c[4] * n + c[5]

			// Set [mConstant], [nConstant], [uCoefs], [vCoefs] with information about [F] and [G]
			// @ uv : The very parameter of [F] that is approximated by [G]
			// @ mn : The matching parameter of [G] to the point F(u, v) : G(m, n) == F(u, v)
			void set(const SurfaceInfo& surface, const TorusApprox& torus, Real m, Real n);

			// For given parameter [uv] in [F]'s domain, evaluate matching parameter [mn] in [G]'s domain
			void calMN(Real u, Real v, Real& m, Real& n) const;
			// For given parameter [mn] in [G]'s domain, evaluate matching parameter [uv] in [F]'s domain
			void calUV(Real m, Real n, Real& u, Real& v) const;

			static Real evaluate(const std::array<Real, 6>& c, Real u, Real v);

			// Calculate M, N domains for given U, V domains
			// If the resulting domain is valid, return true. Else, return false
			bool calDomainM(const Domain& uDomain, const Domain& vDomain, piDomain& mDomain);
			bool calDomainN(const Domain& uDomain, const Domain& vDomain, piDomain& nDomain);
			
			void calPositionErrorUpperBound(const TorusPatch& patch, const Domain& uDomain, const Domain& vDomain, Real& N1, Real& N2, Real& N3, Real& N4) const;
		};
		Mapping mapping;
		TorusPatch patch;
		bool validPatchDomain;	// Whether [ patch ] domain is valid one

		Transform transform;	// Transform that takes local coordinates in torus to world coordinates
		Transform iTransform;	// Transform that takes world coordinates to local coordinates in torus
		
		Real positionError;		// Upper bound of position error to original surface
		Real normalError;		// Upper bound of normal error with original surface

		static TorusApprox create(const SurfaceInfo& surface);
		static Ptr createPtr(const SurfaceInfo& surface);
		void setTransform(const Vec3& center, const Vec3& axis);
	};
}

#endif