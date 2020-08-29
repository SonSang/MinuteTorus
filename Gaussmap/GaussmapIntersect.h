/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_GAUSSMAP_INTERSECT_H__
#define __MN_GAUSSMAP_INTERSECT_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "Gaussmap.h"

namespace MN {
	// Gaussmap - Gaussmap intersection
	// @ a, b : Gaussmaps to intersect
	// @ btoa : Transform that takes a local coordinate of [ b ] to that of [ a ]
	// @ uDomainB : Intersecting domains in [ b ]'s U parameter space ( maximum 3 distinct domains can occur )
	// @ domainNum : Number of [ uDomainB ]
	// @ ret : Always return 1
	int intersect(const Gaussmap& a, const Gaussmap& b, const Transform& btoa, piDomain uDomainB[3], int& domainNum);

	// Gaussmap - Gaussmap intersection
	// @ a, b : Gaussmaps to intersect
	// @ ta, tb : Transform that takes [ a, b ] to world coordinates
	// @ uDomainA, uDomainB : Intersecting domains in U parameter space ( maximum 3 distinct domains can occur )
	// @ aNum, bNum : Number of [ uDomainA, uDomainB ]
	// @ ret : Always return 1
	int intersect(const Gaussmap& a, const Gaussmap& b, const Transform& ta, const Transform& tb, piDomain uDomainA[3], piDomain uDomainB[3], int& aNum, int& bNum);

	// Same as above, except ...
	// @ atob, btoa : Transform that takes local coordinates in A to B and vice versa
	int fIntersect(const Gaussmap& a, const Gaussmap& b, const Transform& atob, const Transform& btoa, piDomain uDomainA[3], piDomain uDomainB[3], int& aNum, int& bNum);
}

#endif