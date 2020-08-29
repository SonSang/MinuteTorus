/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_DISTANCE_H__
#define __MN_DISTANCE_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "../MinuteUtils/Utils.h"

namespace MN {
	class Distance {
	public:
		Real length;
		Vec3 pointA;
		Vec3 pointB;
		Real paramA[2];
		Real paramB[2];
	};
}

#endif