/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_CIRCLE_INTERSECT_H__
#define __MN_CIRCLE_INTERSECT_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "Circle.h"

namespace MN {
	// Circle - Plane intersection
	// @ planeNormal : Normal of intersecting plane in [ circle ]'s local coordinates
	// @ planePoint : A point on intersecting plane in [ circle ]'s local coordinates
	// @ intParam : Intersecting parameters of [ circle ]
	// @ intNum : Number of intersecting parameters ( Normally at most 2 )
	// @ ret :	0 = Since [ circle ] is fully included in the given plane, all points on [ circle ] are intersecting points
	//			1 = [ circle ] is not included in the given plane, so at most 2 points of intersection occurred
	int intersect(const Circle& circle, const Vec3& planeNormal, const Vec3& planePoint, Real intParam[2], int& intNum);

	// Arc version of above
	int intersect(const CircularArc& arc, const Vec3& planeNormal, const Vec3& planePoint, Real intParam[2], int& intNum);
}

#endif