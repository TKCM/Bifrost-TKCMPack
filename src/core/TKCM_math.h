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
			result.x = vec3.x / len;
			result.y = vec3.y / len;
			result.z = vec3.z / len;
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
	
	bool AlmostEqual(const Bifrost::Math::float3 A, const Bifrost::Math::float3 B, const float tolerance = 10e-6f){
		return 
			std::abs(A.x - B.x) < tolerance &&
			std::abs(A.y - B.y) < tolerance &&
			std::abs(A.z - B.z) < tolerance;
	}
	
	bool AlmostParallel ( const Bifrost::Math::float3& valA, const Bifrost::Math::float3& valB, float tolerance = 10e-6f) {
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
	
	Bifrost::Math::float3 Lerp ( const Bifrost::Math::float3& A, const Bifrost::Math::float3& B, const float t ){
		Bifrost::Math::float3 result;
		result.x = TKCM::Lerp<float>(A.x, B.x, t);
		result.y = TKCM::Lerp<float>(A.y, B.y, t);
		result.z = TKCM::Lerp<float>(A.z, B.z, t);
		return result;
	}
	
	float UnitAngleTo ( const Bifrost::Math::float3& unitVec0, const Bifrost::Math::float3& unitVec1 ) {
		float dot = TKCM::Dot ( unitVec0, unitVec1 );
		dot = std::max ( dot, -1.0f );
		dot = std::min ( dot, 1.0f );
		if (TKCM::AlmostEqual ( dot, 1.0f )) {
			return 0.0f;
		} else {
			return std::acos ( dot );
		}
	}
	
	Bifrost::Math::float3 Sub(const Bifrost::Math::float3& A, const Bifrost::Math::float3& B){
		Bifrost::Math::float3 result;
		result.x = A.x - B.x;
		result.y = A.y - B.y;
		result.z = A.z - B.z;
		return result;
	}

	Bifrost::Math::float3 Add(const Bifrost::Math::float3& A, const Bifrost::Math::float3& B){
		Bifrost::Math::float3 result;
		result.x = A.x + B.x;
		result.y = A.y + B.y;
		result.z = A.z + B.z;
		return result;
	}
	
	Bifrost::Math::float3 Mul(const Bifrost::Math::float3& A, const float& B){
		Bifrost::Math::float3 result;
		result.x = A.x * B;
		result.y = A.y * B;
		result.z = A.z * B;
		return result;
	}

	Bifrost::Math::float3 GetClosestPositionOnLineSegment(
		const Bifrost::Math::float3& P,
		const Bifrost::Math::float3& segmentP0,
		const Bifrost::Math::float3& segmentP1,
		float& segment_t,
		const float tolerance
	){
		Bifrost::Math::float3 linVec = TKCM::Sub(segmentP1, segmentP0);
		Bifrost::Math::float3 pToP0 = TKCM::Sub(P, segmentP0);

		float dot1 = TKCM::Dot(pToP0, linVec);
		if (dot1 < 0.0f){
			segment_t = 0.0f;
		} else if (TKCM::AlmostEqual(dot1, 0.0f, tolerance)){
			segment_t = 0.0f;
		} else{
			float dot2 = TKCM::Dot(linVec, linVec);
			if (dot2 < dot1){
				segment_t = 1.0f;
			} else if (TKCM::AlmostEqual(dot2, dot1, tolerance)){
				segment_t = 1.0f;
			} else{
				segment_t = dot1 / dot2;
			}
		}

		return TKCM::Add(segmentP0, TKCM::Mul(linVec, segment_t));
	}

	Bifrost::Math::float4 GetTriBarycentric(
		const Bifrost::Math::float3& P,
		const Bifrost::Math::float3& triP0,
		const Bifrost::Math::float3& triP1,
		const Bifrost::Math::float3& triP2,
		const float tolerance = 0.0001f,
		const bool normalize = true
	){
		Bifrost::Math::float3 pt0ToStart = TKCM::Sub(P, triP0);
		Bifrost::Math::float3 pt0ToPt1 = TKCM::Sub(triP1, triP0);
		Bifrost::Math::float3 pt0ToPt2 = TKCM::Sub(triP2, triP0);

		Bifrost::Math::float3 triCross = TKCM::Cross(pt0ToPt1, pt0ToPt2);
		float triCrossSqLen = TKCM::Dot(triCross, triCross);

		Bifrost::Math::float3 pt0ToInters = TKCM::Sub(pt0ToStart, TKCM::Mul(triCross, (TKCM::Dot(pt0ToStart, triCross) / triCrossSqLen)));
		float val1;
		float d00, d01, d11, d20, d21;
		d00 = TKCM::Dot(pt0ToPt1, pt0ToPt1);
		d01 = TKCM::Dot(pt0ToPt1, pt0ToPt2);
		d11 = TKCM::Dot(pt0ToPt2, pt0ToPt2);
		d20 = TKCM::Dot(pt0ToInters, pt0ToPt1);
		d21 = TKCM::Dot(pt0ToInters, pt0ToPt2);
		val1 = d00 * d11 - d01 * d01;

		float invDenom = 1.0f / val1;
		float v = (d11 * d20 - d01 * d21) * invDenom;
		float w = (d00 * d21 - d01 * d20) * invDenom;
		float u = 1.0f - v - w;
		if (normalize){
			int neg[3];
			neg[0] = u < 0.0f ? 1 : 0;
			neg[1] = v < 0.0f ? 1 : 0;
			neg[2] = w < 0.0f ? 1 : 0;

			Bifrost::Math::float3 vec3Zero;
			if (neg[0] + neg[1] + neg[2]){
				float ratio;
				if (neg[0]){
					TKCM::GetClosestPositionOnLineSegment(pt0ToStart, pt0ToPt1, pt0ToPt2, ratio, tolerance);
					u = 0.0f;
					if (1.0f <= ratio){
						v = 0.0f;
					} else{
						v = 1.0f - ratio;
					}
				}
				if (neg[1]){
					TKCM::GetClosestPositionOnLineSegment(pt0ToStart, vec3Zero, pt0ToPt2, ratio, tolerance);
					if (1.0f <= ratio){
						u = 0.0f;
					} else{
						u = 1.0f - ratio;
					}
					v = 0.0f;
				}
				if (neg[2]){
					TKCM::GetClosestPositionOnLineSegment(pt0ToStart, vec3Zero, pt0ToPt1, ratio, tolerance);
					if (1.0f <= ratio){
						u = 0.0f;
					} else{
						u = 1.0f - ratio;
					}
					v = ratio;
				}
			}
		}

		float resultW;
		if (1.0f <= u + v){
			resultW = 0.0f;
		} else{
			resultW = 1.0f - (u + v);
		}
		Bifrost::Math::float4 result;
		result.x = u;
		result.y = v;
		result.z = resultW;
		result.w = 0.0f;
		return result;
	}
}