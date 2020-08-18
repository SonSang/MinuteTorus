#include "GaussmapIntersect.h"

namespace MN {
	// For given circle defined by [ X, Y ], ( C(t) = center + rcos(t) * X + rsin(t) * Y )
	// find [ t0 ] such that C(t0)[index] has the largest value of the circle
	inline static Real circleLargestValueT(const Vec3& X, const Vec3& Y, int index) {
		Real
			xval = X[index],
			yval = Y[index];
		if (xval == 0.0 && yval == 0.0)
			return 0;
		else if (xval == 0.0 && yval != 0.0) {
			if (yval > 0)
				return PI05;
			else
				return PI15;
		}
		else if (xval != 0 && yval == 0.0) {
			if (xval > 0)
				return 0;
			else
				return PI;
		}
		else {
			Real
				tant = yval / xval,
				t = atan(tant);
			Real
				det = cos(t) * xval + sin(t) * yval;
			if (det > 0)
				return t;
			else
				return (t + PI);
		}
	}
	inline static Real gaussmapLargestZT(const Transform& transform) {
		Vec3 X, Y;
		X = { transform.R[0][0], transform.R[1][0], transform.R[2][0] };
		Y = { transform.R[0][1], transform.R[1][1], transform.R[2][1] };
		return circleLargestValueT(X, Y, 2);
	}
	
	// Intersect [ gaussmap ] with the circle defined by [ center, xvec, yvec ]
	// @ gmUpperZ, gmLowerZ : Z value of upper and lower boundary circle of [ gaussmap ]
	inline static int intersect(const Gaussmap& gaussmap, Real gmUpperZ, Real gmLowerZ,
								const Vec3& center, const Vec3& xvec, const Vec3& yvec,
								Domain circleDomain[2], int& domainNum) {
		bool
			upperBdryInter = false,
			lowerBdryInter = false;
		if (xvec[2] == 0.0) {
			if (isbet(gmLowerZ, gmUpperZ, center[2])) {
				circleDomain[0].set(0, PI);
				circleDomain[1].set(PI, PI20);
				domainNum = 2;
			}
			else
				domainNum = 0;
			return 1;
		}

		Real
			costUpperInter = (gmUpperZ - center[2]) / xvec[2],
			costLowerInter = (gmLowerZ - center[2]) / xvec[2];

		if (fabs(costUpperInter) > 1) {
			if (fabs(costUpperInter) < 1 + 1e-10)
				costUpperInter = (costUpperInter > 0) ? 1 : -1;
		}
		if (fabs(costLowerInter) > 1) {
			if (fabs(costLowerInter) < 1 + 1e-10)
				costLowerInter = (costLowerInter > 0) ? 1 : -1;
		}

		upperBdryInter = (costUpperInter >= -1 && costUpperInter <= 1);
		lowerBdryInter = (costLowerInter >= -1 && costLowerInter <= 1);

		if (!upperBdryInter && !lowerBdryInter) {
			Real det = center[2] + xvec[2];
			if (det >= gmLowerZ && det <= gmUpperZ) {
				circleDomain[0].set(0, PI);
				circleDomain[1].set(PI, PI20);
				domainNum = 2;
			}
			else
				domainNum = 0;
		}
		else if (upperBdryInter && lowerBdryInter) {
			Real
				upperInterT0 = acos(costUpperInter),
				lowerInterT0 = acos(costLowerInter),
				upperInterT1 = PI20 - upperInterT0,
				lowerInterT1 = PI20 - lowerInterT0;
			circleDomain[0].set(upperInterT0, lowerInterT0);
			circleDomain[1].set(lowerInterT1, upperInterT1);
			domainNum = 2;
		}
		else if (upperBdryInter) {
			Real
				upperInterT0 = acos(costUpperInter),
				upperInterT1 = PI20 - upperInterT0;
			circleDomain[0].set(upperInterT0, PI);
			circleDomain[1].set(PI, upperInterT1);
			domainNum = 2;
		}
		else {
			Real
				lowerInterT0 = acos(costLowerInter),
				lowerInterT1 = PI20 - lowerInterT0;
			circleDomain[0].set(0, lowerInterT0);
			circleDomain[1].set(lowerInterT1, PI20);
			domainNum = 2;
		}
		return 1;
	}

	// Add intersecting part between new domain [ begin, end ] and [ validDomain ] to [ domain ]
	// [ begin, end ] must be in [ 0, PI20 ]
	void addDomain(const piDomain& validDomain, Real begin, Real end, piDomain domain[3], int& domainNum) {
		bool hasBegin, hasEnd;
		hasBegin = validDomain.has(begin);
		hasEnd = validDomain.has(end);

		if (hasBegin && hasEnd)
			domain[domainNum++].set(begin, end);
		else if (hasBegin)
			domain[domainNum++].set(begin, piDomain::regularize(validDomain.end()));
		else if (hasEnd)
			domain[domainNum++].set(piDomain::regularize(validDomain.beg()), end);
		else {
			Real vdbegin = piDomain::regularize(validDomain.beg());
			if (begin <= vdbegin && vdbegin <= end) 
				domain[domainNum++] = validDomain;
		}
	}
	int intersect(const Gaussmap& a, const Gaussmap& b, const Transform& btoa, piDomain uDomainB[3], int& domainNum) {
		Vec3
			bUpperCenter{ 0, 0, cos(b.vDomain.beg()) },
			bLowerCenter{ 0, 0, cos(b.vDomain.end()) };
		bUpperCenter = btoa.applyR(bUpperCenter);
		bLowerCenter = btoa.applyR(bLowerCenter);

		// sParam : [ s0 ] such that boundary circle of [ b ] has largest Z value in [ a ]'s local coordinates
		Real
			sParam = gaussmapLargestZT(btoa);

		Vec3
			bUpperSvec{ cos(sParam), sin(sParam), 0 },
			bUpperTvec{ -bUpperSvec[1], bUpperSvec[0], 0 },
			bLowerSvec = bUpperSvec,
			bLowerTvec = bUpperTvec;
		bUpperSvec *= sin(b.vDomain.beg());
		bUpperTvec *= sin(b.vDomain.beg());
		bLowerSvec *= sin(b.vDomain.end());
		bLowerTvec *= sin(b.vDomain.end());

		bUpperSvec = btoa.applyR(bUpperSvec);
		bUpperTvec = btoa.applyR(bUpperTvec);
		bLowerSvec = btoa.applyR(bLowerSvec);
		bLowerTvec = btoa.applyR(bLowerTvec);

		Real
			aUpperZ = cos(a.vDomain.beg()),
			aLowerZ = cos(a.vDomain.end());

		Domain
			bUpperDomain[2],
			bLowerDomain[2];
		int
			bUpperDomainNum,
			bLowerDomainNum;

		intersect(a, aUpperZ, aLowerZ, bUpperCenter, bUpperSvec, bUpperTvec, bUpperDomain, bUpperDomainNum);
		intersect(a, aUpperZ, aLowerZ, bLowerCenter, bLowerSvec, bLowerTvec, bLowerDomain, bLowerDomainNum);

		// Sum result domains.
		Domain domainSum[2];
		if (bUpperDomainNum == 2 && bLowerDomainNum == 2) {
			domainSum[0] = bUpperDomain[0] + bLowerDomain[0];
			domainSum[1] = bUpperDomain[1] + bLowerDomain[1];
		}
		else if (
			(bUpperDomainNum == 2 && bLowerDomainNum == 0) ||
			(bUpperDomainNum == 0 && bLowerDomainNum == 2)) {
			Real
				bUpperZeroZ = bUpperCenter[2] + bUpperSvec[2],
				bUpperPiZ = bUpperCenter[2] - bUpperSvec[2],
				bLowerZeroZ = bLowerCenter[2] + bLowerSvec[2],
				bLowerPiZ = bLowerCenter[2] - bLowerSvec[2];
			bool
				includeZero = (
					isbet(bUpperZeroZ, bLowerZeroZ, aUpperZ) ||
					isbet(bUpperZeroZ, bLowerZeroZ, aLowerZ)),
				includePi = (
					isbet(bUpperPiZ, bLowerPiZ, aUpperZ) ||
					isbet(bUpperPiZ, bLowerPiZ, aLowerZ));
			domainSum[0] = (bUpperDomainNum == 2) ? bUpperDomain[0] : bLowerDomain[0];
			domainSum[1] = (bUpperDomainNum == 2) ? bUpperDomain[1] : bLowerDomain[1];
			if (includeZero) {
				domainSum[0] = domainSum[0] + Domain::create(0, 0);
				domainSum[1] = domainSum[1] + Domain::create(PI20, PI20);
			}
			if (includePi) {
				domainSum[0] = domainSum[0] + Domain::create(PI, PI);
				domainSum[1] = domainSum[1] + Domain::create(PI, PI);
			}
		}
		else {
			Real
				bUpperZeroZ = bUpperCenter[2] + bUpperSvec[2],
				bLowerZeroZ = bLowerCenter[2] + bLowerSvec[2];
			bool
				include = (
					isbet(bUpperZeroZ, bLowerZeroZ, aUpperZ) ||
					isbet(bUpperZeroZ, bLowerZeroZ, aLowerZ));
			if (include) {
				uDomainB[0].set(0, PI20);
				domainNum = 1;
			}
			else
				domainNum = 0;
			return 1;
		}

		// Change domain from [ s ] to original domain.
		Real begin, end;
		domainNum = 0;
		if (sParam < PI) {
			begin = domainSum[0].beg() + sParam;
			end = domainSum[0].end() + sParam;
			addDomain(b.uDomain, begin, end, uDomainB, domainNum);
			if (domainSum[1].end() <= PI20 - sParam) {
				begin = domainSum[1].beg() + sParam;
				end = domainSum[1].end() + sParam;
				addDomain(b.uDomain, begin, end, uDomainB, domainNum);
			}
			else if (domainSum[1].beg() >= PI20 - sParam) {
				begin = domainSum[1].beg() + sParam - PI20;
				end = domainSum[1].end() + sParam - PI20;
				addDomain(b.uDomain, begin, end, uDomainB, domainNum);
			}
			else {
				begin = domainSum[1].beg() + sParam;
				end = PI20;
				addDomain(b.uDomain, begin, end, uDomainB, domainNum);

				begin = 0;
				end = domainSum[1].end() + sParam - PI20;
				addDomain(b.uDomain, begin, end, uDomainB, domainNum);
			}
		}
		else {
			begin = domainSum[1].beg() + sParam - PI20;
			end = domainSum[1].end() + sParam - PI20;
			addDomain(b.uDomain, begin, end, uDomainB, domainNum);
			if (domainSum[0].end() <= PI20 - sParam) {
				begin = domainSum[0].beg() + sParam;
				end = domainSum[0].end() + sParam;
				addDomain(b.uDomain, begin, end, uDomainB, domainNum);
			}
			else if (domainSum[0].beg() >= PI20 - sParam) {
				begin = domainSum[0].beg() + sParam - PI20;
				end = domainSum[0].end() + sParam - PI20;
				addDomain(b.uDomain, begin, end, uDomainB, domainNum);
			}
			else {
				begin = domainSum[0].beg() + sParam;
				end = PI20;
				addDomain(b.uDomain, begin, end, uDomainB, domainNum);
				
				begin = 0;
				end = domainSum[0].end() + sParam - PI20;
				addDomain(b.uDomain, begin, end, uDomainB, domainNum);
			}
		}
		return 1;
	}
	int intersect(const Gaussmap& a, const Gaussmap& b, const Transform& ta, const Transform& tb, piDomain uDomainA[3], piDomain uDomainB[3], int& aNum, int& bNum) {
		Transform atob, btoa;
		atob = Transform::connect(ta, tb);
		btoa = Transform::connect(tb, ta);
		return fIntersect(a, b, atob, btoa, uDomainA, uDomainB, aNum, bNum);
	}
	int fIntersect(const Gaussmap& a, const Gaussmap& b, const Transform& atob, const Transform& btoa, piDomain uDomainA[3], piDomain uDomainB[3], int& aNum, int& bNum) {
		intersect(a, b, btoa, uDomainB, bNum);
		intersect(b, a, atob, uDomainA, aNum);
		return 1;
	}
}