#include "toFCurve_nodedef.h"

void TKCM::to_fcurve(
	const	Amino::Ptr<Amino::Array<Bifrost::Math::float2>>& cp,
	const	Amino::Ptr<Amino::Array<uint8_t>>&  spanInterpolation,
	const	Amino::Ptr<Amino::Array<uint8_t>>&  curveInterpolation,
	const	Amino::Ptr<Amino::Array<bool>>& lock,
	const	uint8_t& preExtrapolation,
	const	uint8_t& postExtrapolation,
			Amino::MutablePtr<Bifrost::Math::FCurve>& fcurve,
			bool& success
) {
	// 出力を用意
	fcurve = Amino::newMutablePtr<Bifrost::Math::FCurve>();

	// cpの入力数が適切か確認する
	Amino::Array<Bifrost::Math::FCurve::SegmentPoint> fcurvePoint;
	if (!cp || cp->empty()){
		//　不正
		for (int i = 0; i < 2; ++i) {
			Bifrost::Math::FCurve::ControlPoints fcurveCP{ 0, 0, 0, 0, 0, 0};
			fcurvePoint.push_back(Bifrost::Math::FCurve::SegmentPoint(fcurveCP));
		}
		success = false;
	}
	else if (cp->size() == 1) {
		// 不足
		for (int i = 0; i < 2; ++i) {
			Bifrost::Math::FCurve::ControlPoints fcurveCP{ 0, 0, float(i), cp->at(i).y, 0, 0};
			fcurvePoint.push_back(Bifrost::Math::FCurve::SegmentPoint(fcurveCP));
		}
		success = false;
	}
	else if (cp->size() == 2) {
		// ハンドルの値が不足
		for (int i = 0; i < 2; ++i) {
			Bifrost::Math::FCurve::ControlPoints fcurveCP{ 0, 0, cp->at(i).x, cp->at(i).y, 0, 0 };
			fcurvePoint.push_back(Bifrost::Math::FCurve::SegmentPoint(fcurveCP));
		}
		success = false;
	}
	else if (cp->size() < 6 || cp->size()%3 != 0) {
		// ハンドルの入力数が適数ではない
		for (int i = 0; i < cp->size(); ++i) {
			Bifrost::Math::FCurve::ControlPoints fcurveCP{ 0, 0, cp->at(i).x, cp->at(i).y, 0, 0 };
			fcurvePoint.push_back(Bifrost::Math::FCurve::SegmentPoint(fcurveCP));
		}
		success = false;
	}
	else {
		// ハンドルを含むcpの数が適数
		int cpCnt = cp->size() / 3;
		bool prepared = (cpCnt == spanInterpolation->size()) && (cpCnt == curveInterpolation->size()) && (cpCnt == lock->size());
		for (int i = 0; i < cp->size()/3; ++i) {
			Bifrost::Math::FCurve::ControlPoints fcurveCP{ cp->at(i).x, cp->at(i).y, cp->at(i+1).x, cp->at(i+1).y, cp->at(i+2).x, cp->at(i+2).y };
			if (prepared) {
				// cpごとのオプション値が適数入力されている場合は適用する
				Bifrost::Math::FCurve::SpanInterpolationMethod spanEr = static_cast<Bifrost::Math::FCurve::SpanInterpolationMethod>(spanInterpolation->at(i));
				Bifrost::Math::FCurve::CurveInterpolationMethod curveEr = static_cast<Bifrost::Math::FCurve::CurveInterpolationMethod>(curveInterpolation->at(i));
				fcurvePoint.push_back(Bifrost::Math::FCurve::SegmentPoint(fcurveCP, spanEr, curveEr, lock->at(i)));
			}
			else {
				// 適数ではない場合は初期値の[Curve, Bezier, lockOff]を適用する
				fcurvePoint.push_back(Bifrost::Math::FCurve::SegmentPoint(fcurveCP));
			}
		}
		success = prepared;
	}

	fcurve->setPoints(fcurvePoint);

	// 前方と後方のカーブのサイクル設定を適用する
	fcurve->setExtrapolationModes(
		static_cast<Bifrost::Math::FCurve::ExtrapolationMode>(preExtrapolation),
		static_cast<Bifrost::Math::FCurve::ExtrapolationMode>(postExtrapolation)
	);
}
