#pragma once
#include <vector>
#include <sstream>
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>

namespace TKCM {
	float Dot( const Bifrost::Math::float3& A, const Bifrost::Math::float3& B ) {
		return A.x * B.x + A.y * B.y + A.z * B.z;
	}

	float LengthSquared (const Bifrost::Math::float3& A, const Bifrost::Math::float3& B){
		return Dot ( A, B );
	}
}