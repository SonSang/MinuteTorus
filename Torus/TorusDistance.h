/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#ifndef __MN_TORUS_DISTANCE_H__
#define __MN_TORUS_DISTANCE_H__

#ifdef _MSC_VER
#pragma once
#endif

#include "..//Distance.h"
#include "Torus.h"

namespace MN {
	class TorusDistance {
	public:
		static Distance minDistanceLocal(const TorusPatch& a, const TorusPatch& b, const Transform& ta, const Transform& tb, Real2 pa, Real2 pb);
		static Distance fMinDistanceLocal(const TorusPatch& a, const TorusPatch& b, const Transform& atob, const Transform& btoa, Real2 pa, Real2 pb);
	};
}

#endif