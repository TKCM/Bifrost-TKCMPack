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

	float LengthSquared (const Bifrost::Math::float3& A){
		return Dot ( A, A );
	}
	
	float Length (const Bifrost::Math::float3& A){
		return std::sqrt (LengthSquared ( A ));
	}
	
	float DistanceSquared (const Bifrost::Math::float3& A, const Bifrost::Math::float3& B){
		Bifrost::Math::float3 v0;
		v0.x = A.x - B.x;
		v0.y = A.y - B.y;
		v0.z = A.z - B.z;
		
		return LengthSquared ( v0 );
	}
	
	float Distance (const Bifrost::Math::float3& A, const Bifrost::Math::float3& B){
		Bifrost::Math::float3 v0;
		v0.x = A.x - B.x;
		v0.y = A.y - B.y;
		v0.z = A.z - B.z;
		
		return Length ( v0 );
	}
	
	bool AlmostEqual(const float A, const float B, const float tolerance = 0.0001f){
		return std::abs(A - B) < tolerance;
	}
}