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
	// Bifrost nodedef
	TKCM_DECL
	void to_fcurve(
		const	Amino::Ptr<Amino::Array<Bifrost::Math::float2>>& cp,
		const	Amino::Ptr<Amino::Array<uint8_t>>& spanInterpolation,
		const	Amino::Ptr<Amino::Array<uint8_t>>& curveInterpolation,
		const	Amino::Ptr<Amino::Array<bool>>& lock,
		const	uint8_t& preExtrapolation,
		const	uint8_t& postExtrapolation,
				Amino::MutablePtr<Bifrost::Math::FCurve>& fcurve,
				bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Math::to_fcurve "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);
}

#endif // BIF_TOFCURVE_H