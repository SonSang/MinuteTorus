/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_CIRCLE_BINORMAL_H__
#define __MN_CIRCLE_BINORMAL_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "Circle.h"

namespace MN {
	class CircleBinormal {
	public:
		struct Binormal {
		public:
			int type;			// Type of this binormal 
								// 0 : Common binormal, which consists of two distinct points from each circle ( or arc )
								// 1 : When two circles share same axis and center ( on XY plane ), every pair of points can form binormal
								// 2 : When a point on circle B is on A's center ( in space ) and also has first deriv parallel to A's axis (RARE)
								//		Then that point can form binormal line with every point on the circle A
								// 3 : When a point on circle A is on B's center ( in space ) and also has first deriv parallel to B's axis (RARE)
								//		Then that point can form binormal line with every point on the circle B

			Real distance;		// Distance between [ pointA ] on arc A and [ pointB ] on arc B in global coordinates
			Vec3 pointA;		// [ pointA ] in its local coordinates
			Vec3 pointB;		// [ pointB ] in its local coordinates
			Real paramA;		// Parameter of [ pointA ] 
			Real paramB;		// Parameter of [ pointB ]
		};
		struct BP {				// Bezier Polynomial, used for solving 8-th degree polynomial 
			Real	coefs[9];
			Real	dCoefs[8];
			int		degree;
			Real	domain[2];

			inline void setCoefs(Real coefs[9], int degree) {
				memcpy(this->coefs, coefs, sizeof(Real) * 9);
				this->degree = degree;
			}
			inline void setDcoefs() {
				for (int i = 0; i < degree; i++)
					dCoefs[i] = coefs[i + 1] - coefs[i];
			}
			inline void operator=(const BP& bp) {
				memcpy(this->coefs, bp.coefs, sizeof(Real) * 9);
				memcpy(this->domain, bp.domain, sizeof(Real) * 2);
				this->degree = bp.degree;
			}
		};
	private:
		std::vector<BP>		data;
		std::vector<Domain> validDomains;

		// BP functions
		BP		initBP(const Real monoCoefs[9], const std::vector<Domain>& domains);
		void	solveBP(const Real monoCoefs[9], const std::vector<Domain>& domains, Real roots[], int& rootNum);

		// This is where actual search process runs. Assume [ a ] is on XY plane, and [ b ] has smaller domain than [ a ]
		void	subroutine(const CircularArc& a, const CircularArc& b, const Transform& btoa, const std::vector<Domain>& bCosDomains, std::vector<Binormal>& bins, bool refine = true, Real precision = 1e-10);
	
		// Exception 1 : If two circles share same axis and center ( on XY plane ), every point on each circle could form binormal
		//				If that is the case, return true
		bool exceptionA(const CircularArc& a, const CircularArc& b, const Transform& btoa);

		// Exception 2 : If arc B goes through axis of arc A, the rendeavue point could generate a binormal that is not detected by following procedure
		//				Therefore, detect those cases in advance
		void exceptionB(const CircularArc& a, const CircularArc& b, const Transform& btoa, std::vector<Binormal>& bins);
	public:
		// Find binormals between two circles ( or arcs )
		// @tA, tB :	Transformations that map circle [ a, b ] to global coordinates
		// @bins :		Result binormals. It contains point information that compose the binormal
		// @refine :	If true, refine binormals until it reaches given precision ( safe )
		//				If false, it returns all possible sets of point pairs that could be binormal without refinement ( unsafe )
		// @precision : Precision to be used for refinement process. 
		//				Precision of binormal line is determined by ===>>> norm(X'(x)) * norm( X(x) - Y(y) ) + norm(Y'(y)) * norm(( X(x) - Y(y) ))
		void solve(const Circle& a, const Circle& b, const Transform& tA, const Transform& tB, std::vector<Binormal>& bins, bool refine = true, Real precision = 1e-10);
		void solve(const CircularArc& a, const CircularArc& b, const Transform& tA, const Transform& tB, std::vector<Binormal>& bins, bool refine = true, Real precision = 1e-10);

		void solve(const Circle& a, const Circle& b, const Transform& btoa, std::vector<Binormal>& bins, bool refine = true, Real precision = 1e-10);
		void solve(const CircularArc& a, const CircularArc& b, const Transform& btoa, std::vector<Binormal>& bins, bool refine = true, Real precision = 1e-10);

		// @bDomain : Valid domains of circle B, needed to speedup VTorus binormal computation with gaussmap
		void solve(const Circle& a, const Circle& b, const Transform& btoa, const std::vector<piDomain>& bDomain, std::vector<Binormal>& bins, bool refine = true, Real precision = 1e-10);
		void solve(const CircularArc& a, const CircularArc& b, const Transform& btoa, const std::vector<piDomain>& bDomain, std::vector<Binormal>& bins, bool refine = true, Real precision = 1e-10);

		// Function to test validity of above functions
		// Just sample points from [ b ] and find binormals
		void bruteSolve(const CircularArc& a, const CircularArc& b, const Transform& tA, const Transform& tB, std::vector<Binormal>& bins, Real precision = 1e-10, int sampleNum = 1000);
	};
}
#endif