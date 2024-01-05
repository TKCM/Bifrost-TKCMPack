#pragma once
#include <vector>
#include <sstream>
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>
#include "TKCM_fn.h"

namespace TKCM {
	Bifrost::Math::float3 Cross( const Bifrost::Math::float3& A, const Bifrost::Math::float3& B ) {
		Bifrost::Math::float3 result;
		result.x = A.y * B.z - A.z * B.y;
		result.y = A.z * B.x - A.x * B.z;
		result.z = A.x * B.y - A.y * B.x;
		return result;
	}
	
	float Dot( const Bifrost::Math::float3& A, const Bifrost::Math::float3& B ) {
		return A.x * B.x + A.y * B.y + A.z * B.z;
	}

	float LengthSquared (const Bifrost::Math::float3& A){
		return Dot ( A, A );
	}
	
	float Length (const Bifrost::Math::float3& A){
		return std::sqrt (TKCM::LengthSquared ( A ));
	}
	
	Bifrost::Math::float3 Normal(const Bifrost::Math::float3 &vec3){
		if (TKCM::AlmostEqual(TKCM::LengthSquared(vec3), 1.0f)) {
			return vec3;
		} else {
			float len = TKCM::Length(vec3);
			Bifrost::Math::float3 result;
			result.x = vec3.y / len;
			result.y = vec3.z / len;
			result.z = vec3.x / len;
			return result;
		}
	}
	
	float DistanceSquared (const Bifrost::Math::float3& A, const Bifrost::Math::float3& B){
		Bifrost::Math::float3 v0;
		v0.x = A.x - B.x;
		v0.y = A.y - B.y;
		v0.z = A.z - B.z;
		
		return TKCM::LengthSquared ( v0 );
	}
	
	float Distance (const Bifrost::Math::float3& A, const Bifrost::Math::float3& B){
		Bifrost::Math::float3 v0;
		v0.x = A.x - B.x;
		v0.y = A.y - B.y;
		v0.z = A.z - B.z;
		
		return Length ( v0 );
	}
	
	bool AlmostEqual(const Bifrost::Math::float3 A, const Bifrost::Math::float3 B, const float tolerance = 0.0001f){
		return 
			std::abs(A.x - B.x) < tolerance &&
			std::abs(A.y - B.y) < tolerance &&
			std::abs(A.z - B.z) < tolerance;
	}
	
	bool AlmostParallel ( const Bifrost::Math::float3& valA, const Bifrost::Math::float3& valB, float tolerance = 0.0001f) {
		if (TKCM::AlmostEqual ( TKCM::LengthSquared ( valA ), 1.0f, tolerance ) &&
			TKCM::AlmostEqual ( TKCM::LengthSquared ( valB ), 1.0f, tolerance )
		){
			float absDot = std::abs ( TKCM::Dot ( valA, valB ) );
			return TKCM::AlmostEqual ( absDot, 1.0f, tolerance );
		} else {
			Bifrost::Math::float3 valAUnit = TKCM::Normal ( valA );
			Bifrost::Math::float3 valBUnit = TKCM::Normal ( valB );
			float absDot = std::abs ( TKCM::Dot ( valAUnit, valBUnit ) );
			return TKCM::AlmostEqual ( absDot, 1.0f, tolerance ) || TKCM::AlmostEqual ( absDot, -1.0f, tolerance );
		}
	}
	
	float GetTriangleArea ( const Bifrost::Math::float3 &p0, const Bifrost::Math::float3 &p1, const Bifrost::Math::float3 &p2 ){
		if (TKCM::AlmostEqual ( p0, p1 ) ||
			TKCM::AlmostEqual ( p1, p2 ) ||
			TKCM::AlmostEqual ( p2, p0 )
		){
			return 0.0f;
		}

		Bifrost::Math::float3 edge1, edge2;
		edge1.x = p1.x - p0.x;
		edge1.y = p1.y - p0.y;
		edge1.z = p1.z - p0.z;
		edge2.x = p2.x - p0.x;
		edge2.y = p2.y - p0.y;
		edge2.z = p2.z - p0.z;
		
		if (TKCM::AlmostParallel ( edge1, edge2 )) {
			return 0.0f;
		}
		return TKCM::Length ( TKCM::Cross ( edge1, edge2 ) ) * 0.5f;
	}
}