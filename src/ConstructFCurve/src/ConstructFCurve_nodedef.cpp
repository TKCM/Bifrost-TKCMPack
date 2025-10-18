#include "ConstructFCurve_nodedef.h"

/*
* FCurve.setPoints(FCurvePoints)でFカーブデータを作成
*	FCurvePoints = Amino::Array<SegmentPoint>
*		SegmentPoint：fカーブの頂点データ
　			ControlPoints：前方のハンドル座標xy、頂点位置座標xy、後方ハンドル座標xy
　			SpanInterpolationMethod：頂点のファンクションタイプ
 　				CurveInterpolationMethod：カーブタイプ（ファンクションタイプがCurveの場合に使用する）
			bool：ロック
*/

void TKCM::construct_fcurve_core(
	const	Amino::Ptr<Amino::Array<TKCM::FCurvePoint>>& cps,
	const	TKCM::FCurveExtrapolationMode& preExtrapolation,
	const	TKCM::FCurveExtrapolationMode& postExtrapolation,
	Amino::MutablePtr<Bifrost::Math::FCurve>& fcurve
) {
	// 出力を用意
	fcurve = Amino::newMutablePtr<Bifrost::Math::FCurve>();

	// cpの入力数が適切か確認する
	Amino::Array<Bifrost::Math::FCurve::SegmentPoint> fcurvePoint;
	if (!cps || cps->empty()) {
		//　不正
		for (int i = 0; i < 2; ++i) {
			Bifrost::Math::FCurve::ControlPoints fcurveCP{ 0, 1, 0, 1, 0, 1 };
			fcurvePoint.push_back(Bifrost::Math::FCurve::SegmentPoint(fcurveCP));
		}
	} else {
		for (int i = 0; i < cps->size(); ++i) {
			Bifrost::Math::FCurve::ControlPoints fcurveCP{
				cps->at(i).in_handle_cp.x,
				cps->at(i).in_handle_cp.y,
				cps->at(i).cp.x,
				cps->at(i).cp.y,
				cps->at(i).out_handle_cp.x,
				cps->at(i).out_handle_cp.y
			};
			fcurvePoint.push_back(
				Bifrost::Math::FCurve::SegmentPoint(
					fcurveCP,
					toAminoEnum(cps->at(i).spanInterpMethod),
					toAminoEnum(cps->at(i).curveMethod),
					cps->at(i).lock
				)
			);
		}
	}
	fcurve->setPoints(fcurvePoint);

	// 前方と後方のカーブのサイクル設定を適用する
	fcurve->setExtrapolationModes(
		static_cast<Bifrost::Math::FCurve::ExtrapolationMode>(preExtrapolation),
		static_cast<Bifrost::Math::FCurve::ExtrapolationMode>(postExtrapolation)
	);
}


