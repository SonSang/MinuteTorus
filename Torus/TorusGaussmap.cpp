/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#include "TorusGaussmap.h"

namespace MN {
	const piDomain TorusPatchGaussmap::h1 = piDomain::create(-PI05, PI05);
	const piDomain TorusPatchGaussmap::h2 = piDomain::create(PI05, PI15);
	const bool TorusPatchGaussmap::uInversion[4] = { false, true, true, false };
	const bool TorusPatchGaussmap::iuInversion[4] = { true, false, false, true };
}