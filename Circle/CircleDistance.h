#ifndef __MN_CIRCLE_DISTANCE_H__
#define __MN_CIRCLE_DISTANCE_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "Circle.h"
#include "../Distance.h"

namespace MN {
	// Find minimum distance between two circles by Vranek's algorithm
	Distance distance(const Circle& a, const Circle& b, const Transform& tA, const Transform& tB);
	// Find minimum distance between two circular arcs by finding binormals
	Distance distance(const CircularArc& a, const CircularArc& b, const Transform& tA, const Transform& tB);
}

#endif