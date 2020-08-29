/*
 *******************************************************************************************
 * Author	: Sang Hyun Son
 * Email	: shh1295@gmail.com
 * Github	: github.com/SonSang
 *******************************************************************************************
 */

#include "TorusBinormal.h"
#include "../Circle/CircleIntersect.h"
#include "../Gaussmap/GaussmapIntersect.h"

#define PROXIMITY_EPS	1e-10
#define SAME_ANGLE_EPS	(1-1e-10)	// If dot product between two vectors exceed this value, they are considered as same vectors 

namespace MN {
	// Full torus
	void TorusBinormal::solve(const Torus& a, const Torus& b, const Transform& ta, const Transform& tb, std::vector<Binormal>& bins) {
		Transform atob, btoa;
		atob = Transform::connect(ta, tb);
		btoa = Transform::connect(tb, ta);
		fSolve(a, b, atob, btoa, bins);
	}
	void TorusBinormal::fSolve(const Torus& a, const Torus& b, const Transform& atob, const Transform& btoa, std::vector<Binormal>& bins) {
		TorusPatch ta, tb;
		ta.majorRadius = a.majorRadius;
		ta.minorRadius = a.minorRadius;
		ta.uDomain.set(0, PI20);
		ta.vDomain.set(0, PI20);

		tb.majorRadius = b.majorRadius;
		tb.minorRadius = b.minorRadius;
		tb.uDomain.set(0, PI20);
		tb.vDomain.set(0, PI20);

		fSolve(ta, tb, atob, btoa, bins);
	}

	// Torus patch
	static void exceptionSameMajorCircle(const TorusPatch& a, const TorusPatch& b, const Transform& atob, const Transform& btoa, std::vector<TorusBinormal::Binormal>& bins) {
		bins.reserve(6);

		//Vec3 apt, aptB;
		Vec3 bpt, bptA;
		TorusBinormal::Binormal bin;
		//bins.reserve(16 + 4);	// (4 X 4) + 4
		piDomain shareDomain[4];
		bool validShareDomain[4];

		if (btoa.R[2][2] > 0) {
			// [ b ]'s axis is in same direction with that of [ a ]
			// Type 3 binormal
			bin.type = 3;

			// Check V domain
			Real2 vpair[4];
			bool validVpair[4] = { false, false, false, false };
			{
				a.vDomain.intersect(b.vDomain, shareDomain, validShareDomain);
				for (int i = 0; i < 4; i++) {
					if (validShareDomain[i]) {
						validVpair[0] = true;
						vpair[0].first = shareDomain[i].beg();
						vpair[0].second = shareDomain[i].beg();
						break;
					}
				}

				a.vDomain.intersect(piDomain::create(b.vDomain.beg() + PI, b.vDomain.width()), shareDomain, validShareDomain);
				for (int i = 0; i < 4; i++) {
					if (validShareDomain[i]) {
						validVpair[1] = true;
						vpair[1].first = shareDomain[i].beg();
						vpair[1].second = piDomain::regularize(shareDomain[i].beg() + PI);
						break;
					}
				}
			}

			// Check U domain
			Real buDomainBegin = b.uDomain.beg();
			bpt = b.majorCircle().evaluate(buDomainBegin);
			bptA = btoa.apply(bpt);

			a.majorCircle().findMinDistParam(bptA, buDomainBegin);
			piDomain buDomainA = piDomain::create(buDomainBegin, buDomainBegin + b.uDomain.width());

			a.uDomain.intersect(buDomainA, shareDomain, validShareDomain);
			for (int i = 0; i < 4; i++) {
				if (validShareDomain[i]) {
					bin.uA = shareDomain[i].beg();
					bin.uB = piDomain::regularize(bin.uA + (b.uDomain.beg() - buDomainA.beg()));

					for (int j = 0; j < 2; j++) {
						if (validVpair[j]) {
							bin.vA = vpair[j].first;
							bin.vB = vpair[j].second;
							bin.pointA = a.evaluate(bin.uA, bin.vA);
							bin.pointB = b.evaluate(bin.uB, bin.vB);
							bin.length = atob.apply(bin.pointA).dist(bin.pointB);
							bins.push_back(bin);
						}
					}
				}
				break;
			}

			// Type 1 binormal
			bin.type = 1;

			// Check V domain
			{
				vpair[0].first = 0; vpair[0].second = 0;
				vpair[1].first = 0; vpair[1].second = PI;
				vpair[2].first = PI; vpair[2].second = 0;
				vpair[3].first = PI; vpair[3].second = PI;

				for (int i = 0; i < 4; i++) {
					if (a.vDomain.has(vpair[i].first) && b.vDomain.has(vpair[i].second))
						validVpair[i] = true;
				}
			}

			// Check U domain
			buDomainA.set(buDomainA.beg() + PI, buDomainA.beg() + PI + buDomainA.width());
			a.uDomain.intersect(buDomainA, shareDomain, validShareDomain);
			for (int i = 0; i < 4; i++) {
				if (validShareDomain[i]) {
					bin.uA = shareDomain[i].beg();
					bin.uB = piDomain::regularize(bin.uA + (b.uDomain.beg() - buDomainA.beg()) + PI);

					for (int j = 0; j < 4; j++) {
						if (validVpair[j]) {
							bin.vA = vpair[j].first;
							bin.vB = vpair[j].second;
							bin.pointA = a.evaluate(bin.uA, bin.vA);
							bin.pointB = b.evaluate(bin.uB, bin.vB);
							bin.length = atob.apply(bin.pointA).dist(bin.pointB);
							bins.push_back(bin);
						}
					}
				}
				break;
			}
		}
		else {
			// [ b ]'s axis is in opposite direction with that of [ a ]
			// Type 3 binormal
			bin.type = 3;

			// Check V domain
			Real2 vpair[4];
			bool validVpair[4] = { false, false, false, false };
			{
				piDomain bvDomainA = piDomain::create(PI20 - b.vDomain.end(), PI20 - b.vDomain.end() + b.vDomain.width());

				a.vDomain.intersect(bvDomainA, shareDomain, validShareDomain);
				for (int i = 0; i < 4; i++) {
					if (validShareDomain[i]) {
						validVpair[0] = true;
						vpair[0].first = shareDomain[i].beg();
						vpair[0].second = PI20 - shareDomain[i].beg();
						break;
					}
				}

				bvDomainA.set(bvDomainA.beg() + PI, bvDomainA.width());
				a.vDomain.intersect(bvDomainA, shareDomain, validShareDomain);
				for (int i = 0; i < 4; i++) {
					if (validShareDomain[i]) {
						validVpair[1] = true;
						vpair[1].first = shareDomain[i].beg();
						vpair[1].second = piDomain::regularize(PI20 - shareDomain[i].beg() + PI);
						break;
					}
				}
			}

			// Check U domain
			Real buDomainBegin = b.uDomain.end();
			bpt = b.majorCircle().evaluate(buDomainBegin);
			bptA = btoa.apply(bpt);

			a.majorCircle().findMinDistParam(bptA, buDomainBegin);
			piDomain buDomainA = piDomain::create(buDomainBegin, buDomainBegin + b.uDomain.width());

			a.uDomain.intersect(buDomainA, shareDomain, validShareDomain);
			for (int i = 0; i < 4; i++) {
				if (validShareDomain[i]) {
					bin.uA = shareDomain[i].beg();
					Real uwidth = bin.uA - buDomainA.beg();
					if (uwidth < 0)
						uwidth += PI20;
					bin.uB = piDomain::regularize(b.uDomain.end() - uwidth);

					for (int j = 0; j < 2; j++) {
						if (validVpair[j]) {
							bin.vA = vpair[j].first;
							bin.vB = vpair[j].second;
							bin.pointA = a.evaluate(bin.uA, bin.vA);
							bin.pointB = b.evaluate(bin.uB, bin.vB);
							bin.length = atob.apply(bin.pointA).dist(bin.pointB);
							bins.push_back(bin);
						}
					}
				}
				break;
			}

			// Type 1 binormal
			bin.type = 1;

			// Check V domain
			{
				vpair[0].first = 0; vpair[0].second = 0;
				vpair[1].first = 0; vpair[1].second = PI;
				vpair[2].first = PI; vpair[2].second = 0;
				vpair[3].first = PI; vpair[3].second = PI;

				for (int i = 0; i < 4; i++) {
					if (a.vDomain.has(vpair[i].first) && b.vDomain.has(vpair[i].second))
						validVpair[i] = true;
				}
			}

			// Check U domain
			buDomainA.set(buDomainA.beg() + PI, buDomainA.beg() + PI + buDomainA.width());
			a.uDomain.intersect(buDomainA, shareDomain, validShareDomain);
			for (int i = 0; i < 4; i++) {
				if (validShareDomain[i]) {
					bin.uA = shareDomain[i].beg();
					Real uwidth = bin.uA - buDomainA.beg();
					if (uwidth < 0)
						uwidth += PI20;
					bin.uB = piDomain::regularize(b.uDomain.end() - uwidth + PI);

					for (int j = 0; j < 4; j++) {
						if (validVpair[j]) {
							bin.vA = vpair[j].first;
							bin.vB = vpair[j].second;
							bin.pointA = a.evaluate(bin.uA, bin.vA);
							bin.pointB = b.evaluate(bin.uB, bin.vB);
							bin.length = atob.apply(bin.pointA).dist(bin.pointB);
							bins.push_back(bin);
						}
					}
				}
				break;
			}
		}
	}
	static void exceptionAlignMajorCircle(const TorusPatch& a, const TorusPatch& b, const Transform& atob, const Transform& btoa, std::vector<TorusBinormal::Binormal>& bins) {
		bins.reserve(8);

		TorusBinormal::Binormal bin;
		Vec3 apt, aptB, bpt, bptA;
		piDomain shareDomain[4];
		bool validShareDomain[4];

		bin.type = 1;

		// Check U domain
		Real2 upair[2];
		bool validUpair[2] = { false, false };

		bool inverted = btoa.R[2][2] < 0;
		Real buDomainBegin;
		piDomain buDomainA;
		if (!inverted) {
			a.majorCircle().findMinDistParam(btoa.apply(b.majorCircle().evaluate(b.uDomain.beg())), buDomainBegin);
			buDomainA.set(buDomainBegin, buDomainBegin + b.uDomain.width());
			a.uDomain.intersect(buDomainA, shareDomain, validShareDomain);

			for (int i = 0; i < 4; i++) {
				if (validShareDomain[i]) {
					validUpair[0] = true;
					upair[0].first = shareDomain[i].beg();
					Real uwidth = shareDomain[i].beg() - buDomainBegin;
					if (uwidth < 0)
						uwidth += PI20;
					upair[0].second = piDomain::regularize(b.uDomain.beg() + uwidth);
					break;
				}
			}
		}
		else {
			a.majorCircle().findMinDistParam(btoa.apply(b.majorCircle().evaluate(b.uDomain.end())), buDomainBegin);
			buDomainA.set(buDomainBegin, buDomainBegin + b.uDomain.width());
			a.uDomain.intersect(buDomainA, shareDomain, validShareDomain);

			for (int i = 0; i < 4; i++) {
				if (validShareDomain[i]) {
					validUpair[1] = true;
					upair[1].first = shareDomain[i].beg();
					Real uwidth = shareDomain[i].beg() - buDomainBegin;
					if (uwidth < 0)
						uwidth += PI20;
					upair[1].second = piDomain::regularize(b.uDomain.end() - uwidth);
					break;
				}
			}
		}

		// Check V domain
		for (int i = 0; i < 2; i++) {
			if (validUpair[i]) {
				bin.uA = upair[i].first;
				bin.uB = upair[i].second;
				
				Real vA[2], vB[2];
				bool validVA[2], validVB[2];
				apt = a.majorCircle().evaluate(bin.uA);
				bpt = b.majorCircle().evaluate(bin.uB);
				aptB = atob.apply(apt);
				bptA = btoa.apply(bpt);

				int ares = a.findExtDistParamV(bptA, bin.uA, vA[0], vA[1]);
				int bres = b.findExtDistParamV(aptB, bin.uB, vB[0], vB[1]);

				validVA[0] = (ares == 3 || ares == 4);
				validVA[1] = (ares == 2 || ares == 4);
				validVB[0] = (bres == 3 || bres == 4);
				validVB[1] = (bres == 2 || bres == 4);

				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						if (validVA[j] && validVB[k]) {
							bin.vA = vA[j];
							bin.vB = vB[k];
							bin.pointA = a.evaluate(bin.uA, bin.vA);
							bin.pointB = b.evaluate(bin.uB, bin.vB);
							bin.length = atob.apply(bin.pointA).dist(bin.pointB);
							bins.push_back(bin);
						}
					}
				}
			}
		}
	}
	static void exceptionSameMinorCircleCenter(const TorusPatch& a, const TorusPatch& b, const Transform& atob, const Transform& btoa, Real uA, Real uB, std::vector<TorusBinormal::Binormal>& bins) {
		TorusBinormal::Binormal bin;
		Vec3 apt, bpt, aptB, bptA;
		
		apt = a.majorCircle().evaluate(uA);
		bpt = b.majorCircle().evaluate(uB);
		aptB = atob.apply(apt);
		bptA = btoa.apply(bpt);
		
		Transform tmp, minorBtoA, minorAtoB;	// Transform that takes minor circle B's local coordinates to minor circle A's local coordinates
		tmp = b.uTransform(uB);
		minorBtoA = tmp.inverse();
		minorBtoA.update(btoa);
		
		tmp = a.uTransform(uA);
		minorBtoA.update(tmp);
		minorAtoB = minorBtoA.inverse();

		// Minor circle's axis are parallel
		if (fabs(minorBtoA.R[2][2]) > SAME_ANGLE_EPS) {
			bool inverted = minorBtoA.R[2][2] < 0;
			piDomain shareDomain[4];
			bool validShareDomain[4];

			bin.type = 2;
			bin.uA = uA;
			bin.uB = uB;

			// Same 
			piDomain bDomainV;
			{
				Vec3 beg = minorBtoA.apply(b.minorCircle().evaluate(inverted ? b.vDomain.end() : b.vDomain.beg()));
				Real begParam;
				a.minorCircle().findMinDistParam(beg, begParam);
				bDomainV.set(begParam, begParam + b.vDomain.width());
			}
			a.vDomain.intersect(bDomainV, shareDomain, validShareDomain);
			for (int i = 0; i < 4; i++) {
				if (validShareDomain[i]) {
					bin.vA = shareDomain[i].beg();
					Real vwidth = shareDomain[i].beg() - bDomainV.beg();
					if (vwidth < 0)
						vwidth += PI20;
					bin.vB = piDomain::regularize((inverted ? b.vDomain.end() - vwidth : b.vDomain.beg() + vwidth));
					bin.pointA = a.evaluate(bin.uA, bin.vA);
					bin.pointB = b.evaluate(bin.uB, bin.vB);
					bin.length = atob.apply(bin.pointA).dist(bin.pointB);
					bins.push_back(bin);
					break;
				}
			}

			// Opposite
			bDomainV.set(bDomainV.beg() + PI, bDomainV.beg() + PI + bDomainV.width());
			a.vDomain.intersect(bDomainV, shareDomain, validShareDomain);
			for (int i = 0; i < 4; i++) {
				if (validShareDomain[i]) {
					bin.vA = shareDomain[i].beg();
					Real vwidth = shareDomain[i].beg() - bDomainV.beg();
					if (vwidth < 0)
						vwidth += PI20;
					bin.vB = piDomain::regularize((inverted ? b.vDomain.end() - vwidth + PI : b.vDomain.beg() + vwidth + PI));
					bin.pointA = a.evaluate(bin.uA, bin.vA);
					bin.pointB = b.evaluate(bin.uB, bin.vB);
					bin.length = atob.apply(bin.pointA).dist(bin.pointB);
					bins.push_back(bin);
					break;
				}
			}
		}
		// Not parallel : Circle - plane intersection
		else {
			bin.type = 0;
			bin.uA = uA;
			bin.uB = uB;

			int vaNum;
			Real vA[2], vB[2];
			Vec3 bNormal = { minorBtoA.R[0][2], minorBtoA.R[1][2], minorBtoA.R[2][2] };
			intersect(a.minorCircle(), bNormal, Vec3::zero(), vA, vaNum);

			if (a.vDomain.has(vA[0])) {
				apt = a.minorCircle().evaluate(vA[0]);
				aptB = minorAtoB.apply(apt);
				int bres = b.minorCircularArc().findExtDistParam(aptB, vB[0], vB[1]);
				
				bool validVB[2];
				validVB[0] = b.vDomain.has(vB[0]);
				validVB[1] = b.vDomain.has(vB[1]);

				for (int i = 0; i < 2; i++) {
					if (validVB[i]) {
						bin.vA = vA[0];
						bin.vB = vB[i];
						bin.pointA = a.evaluate(bin.uA, bin.vA);
						bin.pointB = b.evaluate(bin.uB, bin.vB);
						bin.length = atob.apply(bin.pointA).dist(bin.pointB);
						bins.push_back(bin);
					}
				}
			}

			if (a.vDomain.has(vA[1])) {
				apt = a.minorCircle().evaluate(vA[1]);
				aptB = minorAtoB.apply(apt);
				int bres = b.minorCircularArc().findExtDistParam(aptB, vB[0], vB[1]);

				bool validVB[2];
				validVB[0] = b.vDomain.has(vB[0]);
				validVB[1] = b.vDomain.has(vB[1]);

				for (int i = 0; i < 2; i++) {
					if (validVB[i]) {
						bin.vA = vA[1];
						bin.vB = vB[i];
						bin.pointA = a.evaluate(bin.uA, bin.vA);
						bin.pointB = b.evaluate(bin.uB, bin.vB);
						bin.length = atob.apply(bin.pointA).dist(bin.pointB);
						bins.push_back(bin);
					}
				}
			}
		}
	}
	static void exceptionAmajorBminorCircleAlign(const TorusPatch& a, const TorusPatch& b, const Transform& atob, const Transform& btoa, Real uA, Real uB, std::vector<TorusBinormal::Binormal>& bins) {
		//Vec3 apt, aptB;
		Vec3 bpt, bptA;
		TorusBinormal::Binormal bin;
		piDomain shareDomain[4];
		bool validShareDomain[4];

		bin.type = 4;
		bin.uB = uB;

		// Check V of A
		Real vA[2];
		bool validVA[2] = { false, false };
		if (a.vDomain.has(0)) {
			vA[0] = 0;
			validVA[0] = true;
		}
		if (a.vDomain.has(PI)) {
			vA[1] = PI;
			validVA[1] = true;
		}
		if (!validVA[0] && !validVA[1])
			return;

		// Check U of A and V of B
		Transform minorBmajorA, majorAminorB;
		minorBmajorA = b.uTransform(uB).inverse();
		minorBmajorA.update(btoa);
		majorAminorB = minorBmajorA.inverse();

		bool inverted = minorBmajorA.R[2][2] < 0;
		Real bvDomainBegin = (inverted ? b.vDomain.end() : b.vDomain.beg());
		bpt = b.minorCircle().evaluate(bvDomainBegin);
		bptA = minorBmajorA.apply(bpt);
		a.majorCircle().findMinDistParam(bptA, bvDomainBegin);
		piDomain bvDomainA = piDomain::create(bvDomainBegin, bvDomainBegin + b.vDomain.width());

		// Same
		a.uDomain.intersect(bvDomainA, shareDomain, validShareDomain);
		for (int i = 0; i < 4; i++) {
			if (validShareDomain[i]) {
				bin.uA = shareDomain[i].beg();
				Real uwidth = bin.uA - bvDomainBegin;
				if (uwidth < 0)
					uwidth += PI20;
				bin.vB = (inverted ? b.vDomain.end() - uwidth : b.vDomain.beg() + uwidth);
				for (int j = 0; j < 2; j++) {
					bin.vA = vA[j];
					bin.pointA = a.evaluate(bin.uA, bin.vA);
					bin.pointB = b.evaluate(bin.uB, bin.vB);
					bin.length = atob.apply(bin.pointA).dist(bin.pointB);
					bins.push_back(bin);
				}
				break;
			}
		}

		// Opposite
		bvDomainBegin = piDomain::regularize(bvDomainA.beg() + PI);
		bvDomainA.set(bvDomainBegin, bvDomainBegin + bvDomainA.width());
		a.uDomain.intersect(bvDomainA, shareDomain, validShareDomain);
		for (int i = 0; i < 4; i++) {
			if (validShareDomain[i]) {
				bin.uA = shareDomain[i].beg();
				Real uwidth = bin.uA - bvDomainBegin;
				if (uwidth < 0)
					uwidth += PI20;
				bin.vB = (inverted ? b.vDomain.end() - uwidth : b.vDomain.beg() + uwidth);
				for (int j = 0; j < 2; j++) {
					bin.vA = vA[j];
					bin.pointA = a.evaluate(bin.uA, bin.vA);
					bin.pointB = b.evaluate(bin.uB, bin.vB);
					bin.length = atob.apply(bin.pointA).dist(bin.pointB);
					bins.push_back(bin);
				}
				break;
			}
		}
	}
	static void exceptionAminorBmajorCircleAlign(const TorusPatch& a, const TorusPatch& b, const Transform& atob, const Transform& btoa, Real uA, Real uB, std::vector<TorusBinormal::Binormal>& bins) {
		//Vec3 bpt, bptA;
		Vec3 apt, aptB;
		TorusBinormal::Binormal bin;
		piDomain shareDomain[4];
		bool validShareDomain[4];

		bin.type = 5;
		bin.uA = uA;

		// Check V of B
		Real vB[2];
		bool validVB[2] = { false, false };
		if (b.vDomain.has(0)) {
			vB[0] = 0;
			validVB[0] = true;
		}
		if (b.vDomain.has(PI)) {
			vB[1] = PI;
			validVB[1] = true;
		}
		if (!validVB[0] && !validVB[1])
			return;

		// Check U of B and V of A
		Transform majorBminorA, minorAmajorB;
		minorAmajorB = a.uTransform(uA).inverse();
		minorAmajorB.update(atob);
		majorBminorA = minorAmajorB.inverse();

		bool inverted = minorAmajorB.R[2][2] < 0;
		Real avDomainBegin = (inverted ? a.vDomain.end() : a.vDomain.beg());
		apt = a.minorCircle().evaluate(avDomainBegin);
		aptB = minorAmajorB.apply(apt);
		b.majorCircle().findMinDistParam(aptB, avDomainBegin);
		piDomain avDomainB = piDomain::create(avDomainBegin, avDomainBegin + a.vDomain.width());

		// Same
		b.uDomain.intersect(avDomainB, shareDomain, validShareDomain);
		for (int i = 0; i < 4; i++) {
			if (validShareDomain[i]) {
				bin.uB = shareDomain[i].beg();
				Real uwidth = bin.uB - avDomainBegin;
				if (uwidth < 0)
					uwidth += PI20;
				bin.vA = (inverted ? a.vDomain.end() - uwidth : a.vDomain.beg() + uwidth);
				for (int j = 0; j < 2; j++) {
					bin.vB = vB[j];
					bin.pointA = a.evaluate(bin.uA, bin.vA);
					bin.pointB = b.evaluate(bin.uB, bin.vB);
					bin.length = atob.apply(bin.pointA).dist(bin.pointB);
					bins.push_back(bin);
				}
				break;
			}
		}

		// Opposite
		avDomainBegin = piDomain::regularize(avDomainB.beg() + PI);
		avDomainB.set(avDomainBegin, avDomainBegin + avDomainB.width());
		b.uDomain.intersect(avDomainB, shareDomain, validShareDomain);
		for (int i = 0; i < 4; i++) {
			if (validShareDomain[i]) {
				bin.uB = shareDomain[i].beg();
				Real uwidth = bin.uB - avDomainBegin;
				if (uwidth < 0)
					uwidth += PI20;
				bin.vA = (inverted ? a.vDomain.end() - uwidth : a.vDomain.beg() + uwidth);
				for (int j = 0; j < 2; j++) {
					bin.vB = vB[j];
					bin.pointA = a.evaluate(bin.uA, bin.vA);
					bin.pointB = b.evaluate(bin.uB, bin.vB);
					bin.length = atob.apply(bin.pointA).dist(bin.pointB);
					bins.push_back(bin);
				}
				break;
			}
		}
	}
	void TorusBinormal::solve(const TorusPatch& a, const TorusPatch& b, const Transform& ta, const Transform& tb, std::vector<Binormal>& bins) {
		Transform atob, btoa;
		atob = Transform::connect(ta, tb);
		btoa = Transform::connect(tb, ta);
		fSolve(a, b, atob, btoa, bins);
	}
	void TorusBinormal::fSolve(const TorusPatch& a, const TorusPatch& b, const Transform& atob, const Transform& btoa, std::vector<Binormal>& bins) {
		Vec3 apt, bpt, aptB, bptA;
		CircularArc arcA, arcB;
		arcA = a.majorCircularArc();
		arcB = b.majorCircularArc();
		Binormal bin;

		bins.clear();
		// Exception 1 : Same center ( on XY plane ), Same axis
		if (fabs(btoa.T[0]) < PROXIMITY_EPS && fabs(btoa.T[1]) < PROXIMITY_EPS) {
			// [ b ]'s center is on the axis of [ a ]
			if (fabs(btoa.R[2][2]) > 1 - PROXIMITY_EPS) {
				// [ b ]'s axis is parallel to that of [ a ]
				if (a.majorRadius == b.majorRadius && fabs(btoa.T[2]) < PROXIMITY_EPS) 
					// [ b ]'s major circle is same with that of [ a ]
					exceptionSameMajorCircle(a, b, atob, btoa, bins);
				else 
					// [ b ]'s major circle is aligned with that of [ a ]
					exceptionAlignMajorCircle(a, b, atob, btoa, bins);
				return;
			}
		}

		std::vector<CircleBinormal::Binormal> mcbins;	// Major circle binormals
		circleBinormal.solve(arcA, arcB, btoa, mcbins);
		
		bins.clear();
		bins.reserve(mcbins.size() * 4);

		for (auto& mcbin : mcbins) {
			apt = mcbin.pointA;
			bpt = mcbin.pointB;
			aptB = atob.apply(apt);
			bptA = btoa.apply(bpt);

			// Exception 2 : Minor circle's centers coincide
			if (apt.dist(bptA) < PROXIMITY_EPS) {
				exceptionSameMinorCircleCenter(a, b, atob, btoa, mcbin.paramA, mcbin.paramB, bins);
				continue;
			}

			// Exception 3 : mcbin's type is 2 or 3 ( cannot be type 1, because such cases are dealt with in Exception 1 )
			if (mcbin.type != 0) {
				if (mcbin.type == 2) 
					exceptionAmajorBminorCircleAlign(a, b, atob, btoa, mcbin.paramA, mcbin.paramB, bins);
				else if (mcbin.type == 3) 
					exceptionAminorBmajorCircleAlign(a, b, atob, btoa, mcbin.paramA, mcbin.paramB, bins);
				continue;
			}

			// Normal Case
			Real vA[2], vB[2];
			bool validA[2], validB[2];
			int ares = a.findExtDistParamV(bptA, mcbin.paramA, vA[0], vA[1]);
			int bres = b.findExtDistParamV(aptB, mcbin.paramB, vB[0], vB[1]);
			validA[0] = (ares == 3 || ares == 4);
			validA[1] = (ares == 2 || ares == 4);
			validB[0] = (bres == 3 || bres == 4);
			validB[1] = (bres == 2 || bres == 4);
			
			bin.type = 0;
			bin.uA = mcbin.paramA;
			bin.uB = mcbin.paramB;

			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					if (validA[i] && validB[j]) {
						bin.vA = vA[i];
						bin.vB = vB[j];
						bin.pointA = a.evaluate(bin.uA, bin.vA);
						bin.pointB = b.evaluate(bin.uB, bin.vB);
						bin.length = atob.apply(bin.pointA).dist(bin.pointB);
						bins.push_back(bin);
					}
				}
			}
		}
	}

	// Torus patch with gaussmap
	void TorusBinormal::solve(const TPatchGmap& a, const TPatchGmap& b, const Transform& ta, const Transform& tb, bool aOutward, bool aInward, bool bOutward, bool bInward, std::vector<Binormal>& bins) {
		Transform atob, btoa;
		atob = Transform::connect(ta, tb);
		btoa = Transform::connect(tb, ta);
		fSolve(a, b, atob, btoa, aOutward, aInward, bOutward, bInward, bins);
	}
	void TorusBinormal::fSolve(const TPatchGmap& a, const TPatchGmap& b, const Transform& atob, const Transform& btoa, bool aOutward, bool aInward, bool bOutward, bool bInward, std::vector<Binormal>& bins) {
		Vec3 apt, bpt, aptB, bptA;
		CircularArc arcA, arcB;
		arcA = a.patch.majorCircularArc();
		arcB = b.patch.majorCircularArc();
		Binormal bin;

		bins.clear();
		// Exception 1 : Same center ( on XY plane ), Same axis
		if (fabs(btoa.T[0]) < PROXIMITY_EPS && fabs(btoa.T[1]) < PROXIMITY_EPS) {
			// [ b ]'s center is on the axis of [ a ]
			if (fabs(btoa.R[2][2]) > 1 - PROXIMITY_EPS) {
				// [ b ]'s axis is parallel to that of [ a ]
				if (a.patch.majorRadius == b.patch.majorRadius && fabs(btoa.T[2]) < PROXIMITY_EPS)
					// [ b ]'s major circle is same with that of [ a ]
					exceptionSameMajorCircle(a.patch, b.patch, atob, btoa, bins);
				else
					// [ b ]'s major circle is aligned with that of [ a ]
					exceptionAlignMajorCircle(a.patch, b.patch, atob, btoa, bins);
				return;
			}
		}
		
		// Narrow down domain to find binormals
		std::vector<piDomain> buDomain;
		buDomain.reserve(3 * 4 * 4);

		piDomain gInterB[3];
		int gInterNumB;
		
		if (aOutward && bOutward) {
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					if (a.validGaussmap[i] && b.validGaussmap[j]) {
						intersect(a.iGaussmap[i], b.gaussmap[j], btoa, gInterB, gInterNumB);
						if (b.uInversion[j]) {
							for (int m = 0; m < gInterNumB; m++)
								gInterB[m].set(gInterB[m].beg() + PI, gInterB[m].beg() + PI+ gInterB[m].width());
						}
						for (int m = 0; m < gInterNumB; m++)
							buDomain.push_back(gInterB[m]);
					}
				}
			}
		}
		if (aOutward && bInward) {
			for (int i = 0; i < 2; i++) {
				for (int j = 2; j < 4; j++) {
					if (a.validGaussmap[i] && b.validGaussmap[j]) {
						intersect(a.iGaussmap[i], b.gaussmap[j], btoa, gInterB, gInterNumB);
						if (b.uInversion[j]) {
							for (int m = 0; m < gInterNumB; m++)
								gInterB[m].set(gInterB[m].beg() + PI, gInterB[m].beg() + PI + gInterB[m].width());
						}
						for (int m = 0; m < gInterNumB; m++)
							buDomain.push_back(gInterB[m]);
					}
				}
			}
		}
		if (aInward && bOutward) {
			for (int i = 2; i < 4; i++) {
				for (int j = 0; j < 2; j++) {
					if (a.validGaussmap[i] && b.validGaussmap[j]) {
						intersect(a.iGaussmap[i], b.gaussmap[j], btoa, gInterB, gInterNumB);
						if (b.uInversion[j]) {
							for (int m = 0; m < gInterNumB; m++)
								gInterB[m].set(gInterB[m].beg() + PI, gInterB[m].beg() + PI + gInterB[m].width());
						}
						for (int m = 0; m < gInterNumB; m++)
							buDomain.push_back(gInterB[m]);
					}
				}
			}
		}
		if (aInward && bInward) {
			for (int i = 2; i < 4; i++) {
				for (int j = 2; j < 4; j++) {
					if (a.validGaussmap[i] && b.validGaussmap[j]) {
						intersect(a.iGaussmap[i], b.gaussmap[j], btoa, gInterB, gInterNumB);
						if (b.uInversion[j]) {
							for (int m = 0; m < gInterNumB; m++)
								gInterB[m].set(gInterB[m].beg() + PI, gInterB[m].beg() + PI + gInterB[m].width());
						}
						for (int m = 0; m < gInterNumB; m++)
							buDomain.push_back(gInterB[m]);
					}
				}
			}
		}

		if (buDomain.size() == 0)
			return;

		std::vector<CircleBinormal::Binormal> mcbins;
		circleBinormal.solve(arcA, arcB, btoa, buDomain, mcbins);

		bins.clear();
		bins.reserve(mcbins.size() * 4);

		for (auto& mcbin : mcbins) {
			apt = mcbin.pointA;
			bpt = mcbin.pointB;
			aptB = atob.apply(apt);
			bptA = btoa.apply(bpt);

			// Exception 2 : Minor circle's centers coincide
			if (apt.dist(bptA) < PROXIMITY_EPS) {
				exceptionSameMinorCircleCenter(a.patch, b.patch, atob, btoa, mcbin.paramA, mcbin.paramB, bins);
				continue;
			}

			// Exception 3 : mcbin's type is 2 or 3 ( cannot be type 1, because such cases are dealt with in Exception 1 )
			if (mcbin.type != 0) {
				if (mcbin.type == 2)
					exceptionAmajorBminorCircleAlign(a.patch, b.patch, atob, btoa, mcbin.paramA, mcbin.paramB, bins);
				else if (mcbin.type == 3)
					exceptionAminorBmajorCircleAlign(a.patch, b.patch, atob, btoa, mcbin.paramA, mcbin.paramB, bins);
				continue;
			}

			// Normal Case
			Real vA[2], vB[2];
			bool validA[2], validB[2];
			int ares = a.patch.findExtDistParamV(bptA, mcbin.paramA, vA[0], vA[1]);
			int bres = b.patch.findExtDistParamV(aptB, mcbin.paramB, vB[0], vB[1]);
			validA[0] = (ares == 3 || ares == 4);
			validA[1] = (ares == 2 || ares == 4);
			validB[0] = (bres == 3 || bres == 4);
			validB[1] = (bres == 2 || bres == 4);

			bin.type = 0;
			bin.uA = mcbin.paramA;
			bin.uB = mcbin.paramB;

			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					if (validA[i] && validB[j]) {
						bin.vA = vA[i];
						bin.vB = vB[j];
						bin.pointA = a.patch.evaluate(bin.uA, bin.vA);
						bin.pointB = b.patch.evaluate(bin.uB, bin.vB);
						bin.length = atob.apply(bin.pointA).dist(bin.pointB);
						bins.push_back(bin);
					}
				}
			}
		}
	}
}

/*a.uDomain.intersect(b.uDomain, shareDomain, validShareDomain);
bin.uDomainA.set(0, PI20);
bin.uDomainB.set(0, PI20);

// Check for same V domain
bool inverted = (btoa.R[2][2] < 0);
piDomain bDomainV = inverted ? piDomain::create(PI20 - b.vDomain.end(), PI20 - b.vDomain.end() + b.vDomain.width()) : b.vDomain;
a.vDomain.intersect(bDomainV, shareDomain, validShareDomain);
for (int i = 0; i < 4; i++) {
	if (validShareDomain[i]) {
		bin.vDomainA = shareDomain[i];
		bin.vDomainB = inverted ?
			piDomain::create(PI20 - shareDomain[i].end(), PI20 - b.vDomain.end() + shareDomain[i].width()) :
			shareDomain[i];

		bin.uA = 0;
		bin.vA = piDomain::normalize((bin.vDomainA.begin() + bin.vDomainA.end()) * 0.5);
		bin.pointA = a.evaluate(bin.uA, bin.vA);

		Vec3 aptB = atob.apply(bin.pointA);
		b.projectOnU(aptB, bin.uB);
		b.projectOnV(aptB, bin.uB, bin.vB);
		bin.pointB = b.evaluate(bin.uB, bin.vB);

		bin.length = aptB.dist(bin.pointB);
		bins.push_back(bin);
	}
}

// Check for opposite V domain
bDomainV = inverted ? piDomain::create(PI20 - b.vDomain.end() + PI, PI20 - b.vDomain.end() + PI + b.vDomain.width()) :
	piDomain::create(b.vDomain.begin() + PI, b.vDomain.begin() + PI + b.vDomain.width());
a.vDomain.intersect(bDomainV, shareDomain, validShareDomain);
for (int i = 0; i < 4; i++) {
	if (validShareDomain[i]) {
		bin.vDomainA = shareDomain[i];
		bin.vDomainB = inverted ?
			piDomain::create(PI20 - shareDomain[i].end() + PI, PI20 - shareDomain[i].end() + PI + shareDomain[i].width()) :
			piDomain::create(shareDomain[i].begin() + PI, shareDomain[i].begin() + PI + shareDomain[i].width());

		bin.uA = 0;
		bin.vA = piDomain::normalize((bin.vDomainA.begin() + bin.vDomainA.end()) * 0.5);
		bin.pointA = a.evaluate(bin.uA, bin.vA);

		Vec3 aptB = atob.apply(bin.pointA);
		b.projectOnU(aptB, bin.uB);
		b.projectOnV(aptB, bin.uB, bin.vB);
		bin.vB = piDomain::normalize(bin.vB + PI);	// Opposite
		bin.pointB = b.evaluate(bin.uB, bin.vB);

		bin.length = aptB.dist(bin.pointB);
		bins.push_back(bin);
	}
}

// Opposite U
bin.type = 1;
Circle ac = a.majorCircle(), bc = b.majorCircle();
bin.uB = 0;
Vec3 bpt = bc.evaluate(bin.uB);
Vec3 bptA = btoa.apply(bpt);
ac.projectParam(bptA, bin.uA);
bin.uA = piDomain::normalize(bin.uA + PI);
Vec3 apt = ac.evaluate(bin.uA);
Vec3 aptB = atob.apply(apt);

Real av[2], bv[2];
bool avValid[2], bvValid[2];
a.projectAllOnV(bptA, bin.uA, av, avValid);
b.projectAllOnV(aptB, bin.uB, bv, bvValid);
for (int i = 0; i < 2; i++) {
	if (avValid[i]) {
		bin.vA = av[i];
		bin.pointA = a.evaluate(bin.uA, bin.vA);
		for (int j = 0; j < 2; j++) {
			if (bvValid[j]) {
				bin.vB = bv[j];
				bin.pointB = b.evaluate(bin.uB, bin.vB);
				bin.length = btoa.apply(bin.pointB).dist(bin.pointA);
				bins.push_back(bin);
			}
		}
	}
}
				}
				else {
				bins.reserve(8);	// 4 + 4

				Binormal bin;
				bin.type = 1;
				bin.uDomainA.set(0, PI20);
				bin.uDomainB.set(0, PI20);

				// Closest
				Circle ac = a.majorCircle(), bc = b.majorCircle();
				Real au, bu = 0;
				Vec3 bpt = bc.evaluate(bu);
				Vec3 bptA = btoa.apply(bpt);
				ac.projectParam(bptA, au);
				Vec3 apt = ac.evaluate(au);
				Vec3 aptB = atob.apply(apt);

				Real av[2], bv[2];
				bool avValid[2], bvValid[2];
				a.projectAllOnV(bptA, au, av, avValid);
				b.projectAllOnV(aptB, bu, bv, bvValid);

				bin.uA = au;
				bin.uB = bu;
				for (int i = 0; i < 2; i++) {
					if (avValid[i]) {
						bin.vA = av[i];
						bin.pointA = a.evaluate(bin.uA, bin.vA);
						for (int j = 0; j < 2; j++) {
							if (bvValid[j]) {
								bin.vB = bv[j];
								bin.pointB = b.evaluate(bin.uB, bin.vB);
								bin.length = btoa.apply(bin.pointB).dist(bin.pointA);
								bins.push_back(bin);
							}
						}
					}
				}

				// Farthest
				au = piDomain::normalize(au + PI);
				apt = ac.evaluate(au);
				aptB = atob.apply(apt);

				a.projectAllOnV(bptA, au, av, avValid);
				b.projectAllOnV(aptB, bu, bv, bvValid);

				bin.uA = au;
				bin.uB = bu;
				for (int i = 0; i < 2; i++) {
					if (avValid[i]) {
						bin.vA = av[i];
						bin.pointA = a.evaluate(bin.uA, bin.vA);
						for (int j = 0; j < 2; j++) {
							if (bvValid[j]) {
								bin.vB = bv[j];
								bin.pointB = b.evaluate(bin.uB, bin.vB);
								bin.length = btoa.apply(bin.pointB).dist(bin.pointA);
								bins.push_back(bin);
							}
						}
					}
				}
				}
				return;
			}
		}

		std::vector<CircleBinormal::Binormal> mcbins;	// Major circle binormals
		circleBinormal.solveFast(arcA, arcB, btoa, mcbins);

		Vec3 apt, bpt, aptB, bptA;
		bins.clear();
		bins.reserve(bins.size() * 4);

		for (auto& mcbin : mcbins) {
			apt = mcbin.pointA;
			bpt = mcbin.pointB;
			aptB = atob.apply(apt);
			bptA = btoa.apply(bpt);

			// Exception 2 : Minor circle's centers coincide
			if (apt.dist(bptA) < PROXIMITY_EPS) {
				Transform tmp, minorBtoA, minorAtoB;	// Transform that takes minor circle B's local coordinates to minor circle A's local coordinates
				b.uTransform(mcbin.paramB, tmp);
				minorBtoA = tmp.inverse();
				minorBtoA.update(btoa);
				a.uTransform(mcbin.paramA, tmp);
				minorBtoA.update(tmp);
				minorAtoB = minorBtoA.inverse();

				Circle circle;
				circle.radius = 1;
				// Minor circle's axis are parallel
				if (fabs(minorBtoA.R[2][2]) > 1 - 1e-10) {
					bool inverted = minorBtoA.R[2][2] < 0;
					piDomain shareDomain[4];
					bool validShareDomain[4];

					Binormal bin;
					bin.type = 2;
					bin.uA = mcbin.paramA;
					bin.uB = mcbin.paramB;

					// Same 
					piDomain bDomainV;
					{
						Vec3 beg = minorBtoA.apply(circle.evaluate(inverted ? b.vDomain.end() : b.vDomain.begin()));
						Real begParam;
						circle.projectParam(beg, begParam);
						bDomainV.set(begParam, begParam + b.vDomain.width());
					}
					a.vDomain.intersect(bDomainV, shareDomain, validShareDomain);
					for (int i = 0; i < 4; i++) {
						if (validShareDomain[i]) {
							bin.vDomainA = shareDomain[i];		// @TODO : Have to determine bin.vDomainB...
							bin.vA = piDomain::normalize((bin.vDomainA.begin() + bin.vDomainA.end()) * 0.5);
							bin.pointA = a.evaluate(bin.uA, bin.vA);
							Vec3 aptB = atob.apply(bin.pointA);
							b.projectOnV(aptB, bin.uB, bin.vB);
							bin.pointB = b.evaluate(bin.uB, bin.vB);
							bin.length = aptB.dist(bin.pointB);
							bins.push_back(bin);
						}
					}

					// Opposite
					{
						Vec3 beg = minorBtoA.apply(circle.evaluate(inverted ? b.vDomain.end() + PI : b.vDomain.begin() + PI));
						Real begParam;
						circle.projectParam(beg, begParam);
						bDomainV.set(begParam, begParam + b.vDomain.width());
					}
					a.vDomain.intersect(bDomainV, shareDomain, validShareDomain);
					for (int i = 0; i < 4; i++) {
						if (validShareDomain[i]) {
							bin.vDomainA = shareDomain[i];		// @TODO : Have to determine bin.vDomainB...
							bin.vA = piDomain::normalize((bin.vDomainA.begin() + bin.vDomainA.end()) * 0.5);
							bin.pointA = a.evaluate(bin.uA, bin.vA);
							Vec3 aptB = atob.apply(bin.pointA);
							b.projectOnV(aptB, bin.uB, bin.vB);
							bin.vB = piDomain::normalize(bin.vB + PI);	// Opposite
							bin.pointB = b.evaluate(bin.uB, bin.vB);
							bin.length = aptB.dist(bin.pointB);
							bins.push_back(bin);
						}
					}
				}
				// Not parallel : Circle - plane intersection
				else {
					Vec3 X, Y;
					X = { minorBtoA.R[0][0], minorBtoA.R[1][0], minorBtoA.R[2][0] };
					Y = { minorBtoA.R[0][1], minorBtoA.R[1][1], minorBtoA.R[2][1] };
					Real bv[2];
					if (X[2] == 0.0) {
						bv[0] = 0;
						bv[1] = PI;
					}
					else if (Y[2] == 0.0) {
						bv[0] = PI05;
						bv[1] = PI15;
					}
					else {	// Cannot be zero at the same time
						Real tanu = -X[2] / Y[2];
						bv[0] = piDomain::normalize(atan(tanu));
						bv[1] = piDomain::normalize(bv[0] + PI);
					}
					Real av[2];
					Vec3 tmp = minorBtoA.apply(circle.evaluate(bv[0]));
					circle.projectParam(tmp, av[0]);
					av[1] = piDomain::normalize(av[0] + PI);

					Binormal bin;
					bin.type = 0;
					bin.uA = mcbin.paramA;
					bin.uB = mcbin.paramB;

					for (int i = 0; i < 2; i++) {
						if (a.vDomain.has(av[i])) {
							bin.vA = av[i];
							bin.pointA = a.evaluate(bin.uA, bin.vA);
							for (int j = 0; j < 2; j++) {
								if (b.vDomain.has(bv[j])) {
									bin.vB = bv[j];
									bin.pointB = b.evaluate(bin.uB, bin.vB);
									bin.length = bin.pointA.dist(btoa.apply(bin.pointB));
									bins.push_back(bin);
								}
							}
						}
					}
				}
				continue;
			}

			// Exception 3 : mcbin's type is 2 or 3 ( cannot be type 1, because such cases are dealt with in Exception 1 )
			if (mcbin.type != 0) {
				Circle circle;
				circle.radius = 1;
				if (mcbin.type == 2) {
					if (!a.vDomain.has(0) && !a.vDomain.has(PI) && !a.vDomain.has(PI20))
						continue;

					Vec3 derivB = arcB.differentiate(mcbin.paramB) * -1.0;	// @BUGFIX : Axis follows inverted tangent vector
					derivB = btoa.applyR(derivB);
					bool inverted = derivB[2] < 0;

					// Same
					Vec3 beg = b.evaluate(mcbin.paramB, (inverted ? b.vDomain.end() : b.vDomain.begin()));
					beg = btoa.apply(beg);
					Real begParam;
					circle.projectParam(beg, begParam);

					piDomain bDomainV = piDomain::create(begParam, begParam + b.vDomain.width());

					Binormal bin;
					bin.type = 4;
					bin.uB = mcbin.paramB;
					bin.vDomainB = b.vDomain;
					bin.vB = piDomain::normalize((b.vDomain.begin() + b.vDomain.end()) * 0.5);
					bin.pointB = b.evaluate(bin.uB, bin.vB);

					if (a.vDomain.has(0) || a.vDomain.has(PI20)) {
						bin.uDomainA = bDomainV;
						bin.vA = 0;

						bin.uA = piDomain::normalize((bin.uDomainA.begin() + bin.uDomainA.end()) * 0.5);
						bin.pointA = a.evaluate(bin.uA, bin.vA);
						bin.length = btoa.apply(bin.pointB).dist(bin.pointA);
						bins.push_back(bin);
					}
					if (a.vDomain.has(PI)) {
						bin.uDomainA = bDomainV;
						bin.vA = PI;

						bin.uA = piDomain::normalize((bin.uDomainA.begin() + bin.uDomainA.end()) * 0.5);
						bin.pointA = a.evaluate(bin.uA, bin.vA);
						bin.length = btoa.apply(bin.pointB).dist(bin.pointA);
						bins.push_back(bin);
					}

					// Opposite
					beg = b.evaluate(mcbin.paramB, (inverted ? b.vDomain.end() + PI : b.vDomain.begin() + PI));
					beg = btoa.apply(beg);
					circle.projectParam(beg, begParam);
					bDomainV = piDomain::create(begParam, begParam + b.vDomain.width());

					if (a.vDomain.has(0) || a.vDomain.has(PI20)) {
						bin.uDomainA = bDomainV;
						bin.vA = 0;

						bin.uA = piDomain::normalize((bin.uDomainA.begin() + bin.uDomainA.end()) * 0.5);
						bin.pointA = a.evaluate(bin.uA, bin.vA);
						bin.length = btoa.apply(bin.pointB).dist(bin.pointA);
						bins.push_back(bin);
					}
					if (a.vDomain.has(PI)) {
						bin.uDomainA = bDomainV;
						bin.vA = PI;

						bin.uA = piDomain::normalize((bin.uDomainA.begin() + bin.uDomainA.end()) * 0.5);
						bin.pointA = a.evaluate(bin.uA, bin.vA);
						bin.length = btoa.apply(bin.pointB).dist(bin.pointA);
						bins.push_back(bin);
					}
				}
				else {
					if (!b.vDomain.has(0) && !b.vDomain.has(PI) && !b.vDomain.has(PI20))
						continue;

					Vec3 derivA = arcA.differentiate(mcbin.paramA) * -1.0;
					derivA = atob.applyR(derivA);
					bool inverted = derivA[2] < 0;

					// Same
					Vec3 beg = a.evaluate(mcbin.paramA, (inverted ? a.vDomain.end() : a.vDomain.begin()));
					beg = atob.apply(beg);
					Real begParam;
					circle.projectParam(beg, begParam);

					piDomain aDomainV = piDomain::create(begParam, begParam + a.vDomain.width());

					Binormal bin;
					bin.type = 4;
					bin.uA = mcbin.paramA;
					bin.vDomainA = a.vDomain;
					bin.vA = piDomain::normalize((a.vDomain.begin() + a.vDomain.end()) * 0.5);
					bin.pointA = a.evaluate(bin.uA, bin.vA);

					if (b.vDomain.has(0) || b.vDomain.has(PI20)) {
						bin.uDomainB = aDomainV;
						bin.vB = 0;

						bin.uB = piDomain::normalize((bin.uDomainB.begin() + bin.uDomainB.end()) * 0.5);
						bin.pointB = a.evaluate(bin.uB, bin.vB);
						bin.length = btoa.apply(bin.pointB).dist(bin.pointA);
						bins.push_back(bin);
					}
					if (b.vDomain.has(PI)) {
						bin.uDomainB = aDomainV;
						bin.vB = PI;

						bin.uB = piDomain::normalize((bin.uDomainB.begin() + bin.uDomainB.end()) * 0.5);
						bin.pointB = a.evaluate(bin.uB, bin.vB);
						bin.length = btoa.apply(bin.pointB).dist(bin.pointA);
						bins.push_back(bin);
					}

					// Opposite
					beg = a.evaluate(mcbin.paramA, (inverted ? a.vDomain.end() + PI : a.vDomain.begin() + PI));
					beg = atob.apply(beg);
					circle.projectParam(beg, begParam);
					aDomainV = piDomain::create(begParam, begParam + a.vDomain.width());

					if (b.vDomain.has(0) || b.vDomain.has(PI20)) {
						bin.uDomainB = aDomainV;
						bin.vB = 0;

						bin.uB = piDomain::normalize((bin.uDomainB.begin() + bin.uDomainB.end()) * 0.5);
						bin.pointB = a.evaluate(bin.uB, bin.vB);
						bin.length = btoa.apply(bin.pointB).dist(bin.pointA);
						bins.push_back(bin);
					}
					if (b.vDomain.has(PI)) {
						bin.uDomainB = aDomainV;
						bin.vB = PI;

						bin.uB = piDomain::normalize((bin.uDomainB.begin() + bin.uDomainB.end()) * 0.5);
						bin.pointB = a.evaluate(bin.uB, bin.vB);
						bin.length = btoa.apply(bin.pointB).dist(bin.pointA);
						bins.push_back(bin);
					}
				}
				continue;
			}

			// Normal Case
			Real vA[2], vB[2];
			bool validA[2], validB[2];
			b.projectAllOnV(aptB, mcbin.paramB, vB, validB);
			a.projectAllOnV(bptA, mcbin.paramA, vA, validA);

			Binormal bin;
			bin.type = 0;
			bin.uA = mcbin.paramA;
			bin.uB = mcbin.paramB;
			for (int i = 0; i < 2; i++) {
				if (validA[i]) {
					bin.vA = vA[i];
					for (int j = 0; j < 2; j++) {
						if (validB[j]) {
							bin.vB = vB[j];
							bin.pointA = a.evaluate(bin.uA, bin.vA);
							bin.pointB = b.evaluate(bin.uB, bin.vB);
							bin.length = bin.pointA.dist(btoa.apply(bin.pointB));
							bins.push_back(bin);
						}
					}
				}
			}
		}*/