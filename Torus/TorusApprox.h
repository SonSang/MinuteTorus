/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_TORUS_APPROX_H__
#define __MN_TORUS_APPROX_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "Torus.h"
#include <memory>
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

			Real M1;			// Maximum | Fuuu | value in given domain
			Real M2;			// Maximum | Fuuv | value in given domain
			Real M3;			// Maximum | Fuvv | value in given domain
			Real M4;			// Maximum | Fvvv | value in given domain
			Domain uDomain;
			Domain vDomain;		// Deriv and domain information is needed to determine error bound
		};
		// This class connects parameter space between original surface [ F(u, v) ] and torus approximation [ G(m, n) ]
		class Mapping {
		public:
			Real m0;
			Real n0;		// Parameter of the approximation point
			std::array<Real, 6> mCoefs;	// Coefficients for [ m ] : m(u, v) = c[0] * u^2 + c[1] * v^2 + c[2] * u * v + c[3] * u + c[4] * v + c[5]
			std::array<Real, 6> nCoefs;	// Coefficients for [ n ] : n(u, v) = c[0] * u^2 + c[1] * v^2 + c[2] * u * v + c[3] * u + c[4] * v + c[5]
			Domain uDomain;
			Domain vDomain;			// Domain of original surface [ F ] that this torus patch approximates

			// Set [mConstant], [nConstant], [uCoefs], [vCoefs] with information about [F] and [G]
			// @ uv : The very parameter of [F] that is approximated by [G]
			// @ mn : The matching parameter of [G] to the point F(u, v) : G(m, n) == F(u, v)
			void set(const SurfaceInfo& surface, const TorusApprox& torus, Real m, Real n);

			// For given parameter [uv] in [F]'s domain, evaluate matching parameter [mn] in [G]'s domain
			void calMN(Real u, Real v, Real& m, Real& n) const;
			inline Real calMu(Real u, Real v) const {
				return (2 * mCoefs[0] * u) + (mCoefs[2] * v) + mCoefs[3];
			}
			inline Real calMv(Real u, Real v) const {
				return (2 * mCoefs[1] * v) + (mCoefs[2] * u) + mCoefs[4];
			}
			inline Real calNu(Real u, Real v) const {
				return (2 * nCoefs[0] * u) + (nCoefs[2] * v) + nCoefs[3];
			}
			inline Real calNv(Real u, Real v) const {
				return (2 * nCoefs[1] * v) + (nCoefs[2] * u) + nCoefs[4];
			}
			inline Real calMuu(Real u, Real v) const {
				return 2 * mCoefs[0];
			}
			inline Real calMuv(Real u, Real v) const {
				return mCoefs[2];
			}
			inline Real calMvv(Real u, Real v) const {
				return 2 * mCoefs[1];
			}
			inline Real calNuu(Real u, Real v) const {
				return 2 * nCoefs[0];
			}
			inline Real calNuv(Real u, Real v) const {
				return nCoefs[2];
			}
			inline Real calNvv(Real u, Real v) const {
				return 2 * nCoefs[1];
			}

			// For given parameter [mn] in [G]'s domain, evaluate matching parameter [uv] in [F]'s domain
			// We use numerical algorithm in here to find them
			// @ret : If there is no valid matching [uv], return false, with last valid [uv] 
			bool calUV(Real m, Real n, Real& u, Real& v, Real eps = 1e-10) const;

			static Real evaluate(const std::array<Real, 6>& c, Real u, Real v);

			// Calculate M, N domains for given U, V domains
			// If the resulting domain is valid, return true. Else, return false
			void calDomainM(const Domain& uDomain, const Domain& vDomain, piDomain& mDomain);
			void calDomainN(const Domain& uDomain, const Domain& vDomain, piDomain& nDomain);
			void calPositionErrorUpperBound(const TorusPatch& patch, const Domain& uDomain, const Domain& vDomain, Real& N1, Real& N2, Real& N3, Real& N4) const;

			/* LEGACY : Inverse mapping is not valid in this simple form!
			//std::array<Real, 6> uCoefs;	// Coefficients for [ u ] : u(m, n) = c[0] * m^2 + c[1] * n^2 + c[2] * m * n + c[3] * m + c[4] * n + c[5]
			//std::array<Real, 6> vCoefs;	// Coefficients for [ v ] : v(m, n) = c[0] * m^2 + c[1] * n^2 + c[2] * m * n + c[3] * m + c[4] * n + c[5]
			// For given parameter [mn] in [G]'s domain, evaluate matching parameter [uv] in [F]'s domain
			// Inaccurate
			//void calUV(Real m, Real n, Real& u, Real& v) const;
			*/
		};
		Mapping mapping;
		TorusPatch patch;
		//bool validPatchDomain;	// Whether [ patch ] domain is valid one

		Transform transform;	// Transform that takes local coordinates in torus to world coordinates
		Transform iTransform;	// Transform that takes world coordinates to local coordinates in torus

		Real pErrorT;			// Upper bound of position error to original surface ( Computed by theory, STRONG but LARGE )
		Real pErrorS;			// Upper bound of position error to original surface ( Computed by sampling, WEAK but SMALL )
		Real nError;			// Upper bound of normal error with original surface

		
		Real k1, k2;	// Principal curvatures
		Vec3 w1, w2;	// Principal directions
		
		static TorusApprox create(const SurfaceInfo& surface);
		static Ptr createPtr(const SurfaceInfo& surface);
		void setTransform(const Vec3& center, const Vec3& axis);

		struct SamplePoint {
			Real u;
			Real v;
			Vec3 point;
		};

		// Set [ pErrorS ] by given [ samplePoints ]
		// @multiplier : We can multiply it to computed upper bound for satisfaction
		void setPositionErrorBySampling(const std::vector<SamplePoint>& samplePoints, Real multiplier);
		void CreateFromSingleArc(Biarc2d::Circle2d singleArc);
	};
}

#endif