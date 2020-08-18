#include "TorusDistance.h"

namespace MN {
	Distance TorusDistance::minDistanceLocal(const TorusPatch& a, const TorusPatch& b, const Transform& ta, const Transform& tb, Real2 pa, Real2 pb) {
		Transform atob, btoa;
		atob = Transform::connect(ta, tb);
		btoa = Transform::connect(tb, ta);
		return fMinDistanceLocal(a, b, atob, btoa, pa, pb);
	}
	Distance TorusDistance::fMinDistanceLocal(const TorusPatch& a, const TorusPatch& b, const Transform& atob, const Transform& btoa, Real2 pa, Real2 pb) {
		const static int itermax = 100;
		Distance distance;

		bool projectOnB = true;
		Vec3 tmppt, diff;
		Real mind = maxDouble, curd;
		int iter = 0;

		distance.pointA = a.evaluate(pa.first, pa.second);

		while (true) {
			if (projectOnB) {
				tmppt = atob.apply(distance.pointA);
				b.findMinDistPoint(tmppt, distance.pointB);
				diff = tmppt - distance.pointB;
			}
			else {
				tmppt = btoa.apply(distance.pointB);
				a.findMinDistPoint(tmppt, distance.pointA);
				diff = tmppt - distance.pointA;
			}
			curd = diff.lensq();
			if (curd < mind) {
				if (iter > 0)
					mind = curd;
			}
			else
				break;
			projectOnB = !projectOnB;
			if (++iter > itermax)
				break;
		}
		distance.length = sqrt(mind);
		return distance;
	}
}