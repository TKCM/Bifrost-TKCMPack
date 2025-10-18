#if 1930 <= _MSC_VER // Visual Studio 2022 (MSVC 14.30 以降)
#define _ALLOW_COMPILER_AND_STL_VERSION_MISMATCH // Clangのバージョンエラーが出るのを無視する
#endif

#ifndef BIF_TOFCURVE_H
#define BIF_TOFCURVE_H
// Bifrost Amino
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>
#include <Amino/Core/Array.h>
#include <Amino/Core/String.h>

#include <Bifrost/Math/FCurve.h>

// TKCM
#include "../../core/TKCM_core.h"
#include "../../core/TKCM_fn.h"
//#include "functions.h"

namespace TKCM {
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	enum class AMINO_ANNOTATE("Amino::Enum") FCurveSpanInterpolationMethod {
		Curve,
		Constant,
		Linear
	};

	constexpr Bifrost::Math::FCurve::SpanInterpolationMethod toAminoEnum(FCurveSpanInterpolationMethod value) {
		switch (value) {
			case FCurveSpanInterpolationMethod::Curve:		return Bifrost::Math::FCurve::SpanInterpolationMethod::Curve;
			case FCurveSpanInterpolationMethod::Constant:	return Bifrost::Math::FCurve::SpanInterpolationMethod::Constant;
			case FCurveSpanInterpolationMethod::Linear:		return Bifrost::Math::FCurve::SpanInterpolationMethod::Linear;
		}
	}

	enum class AMINO_ANNOTATE("Amino::Enum") FCurveCurveInterpolationMethod {
		Bezier,
		Hermite,
		Sine,
		Parabolic,
		TangentLog
	};

	constexpr Bifrost::Math::FCurve::CurveInterpolationMethod toAminoEnum(FCurveCurveInterpolationMethod value) {
		switch (value) {
			case FCurveCurveInterpolationMethod::Bezier:		return Bifrost::Math::FCurve::CurveInterpolationMethod::Bezier;
			case FCurveCurveInterpolationMethod::Hermite:		return Bifrost::Math::FCurve::CurveInterpolationMethod::Hermite;
			case FCurveCurveInterpolationMethod::Sine:			return Bifrost::Math::FCurve::CurveInterpolationMethod::Sine;
			case FCurveCurveInterpolationMethod::Parabolic:		return Bifrost::Math::FCurve::CurveInterpolationMethod::Parabolic;
			case FCurveCurveInterpolationMethod::TangentLog:	return Bifrost::Math::FCurve::CurveInterpolationMethod::TangentLog;
		}
	}

	enum class AMINO_ANNOTATE("Amino::Enum") FCurveExtrapolationMode {
		Constant,
		Linear,
		Cycle,
		RelativeRepeat,
		Oscillate
	};

	struct AMINO_ANNOTATE("Amino::Struct") FCurvePoint {
		Bifrost::Math::float2 in_handle_cp, cp, out_handle_cp;
		FCurveSpanInterpolationMethod spanInterpMethod;
		FCurveCurveInterpolationMethod curveMethod;
		bool lock;
    };

	TKCM_DECL
	void construct_fcurve_core(
		const	Amino::Ptr<Amino::Array<FCurvePoint>>& cps,
		const	FCurveExtrapolationMode& preExtrapolation,
		const	FCurveExtrapolationMode& postExtrapolation,
		Amino::MutablePtr<Bifrost::Math::FCurve>& fcurve
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Math::construct_fcurve_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, false}] "
	);
}

#endif // BIF_TOFCURVE_H

