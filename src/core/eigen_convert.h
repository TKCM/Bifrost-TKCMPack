#pragma once
#include <vector>
#include <sstream>
#include <Amino/Bifrost/Math/Types.h>
#include <Amino/Core/AminoArray.h>
#include <Eigen/Dense>
#include "TKCM_define.h"

namespace TKCM {
	void Convert ( BIF_O Eigen::Vector3f &eigen,					BIF_I Bifrost::Math::float3 &bif );
	void Convert ( BIF_O Bifrost::Math::float3 &bif,				BIF_I Eigen::Vector3f &eigen );
	void Convert ( BIF_O std::vector<Eigen::Vector3f>& eigen,		BIF_I Amino::Array<Bifrost::Math::float3>& bif );
	void Convert ( BIF_O Amino::Array<Bifrost::Math::float3>& bif,	BIF_I std::vector<Eigen::Vector3f>& eigen );

	void Convert ( BIF_O Eigen::Matrix4f &eigen,					BIF_I Bifrost::Math::float4x4 &bif );
	void Convert ( BIF_O Bifrost::Math::float4x4 &bif,				BIF_I Eigen::Matrix4f &eigen );
	void Convert ( BIF_O std::vector<Eigen::Matrix4f> &eigen,		BIF_I Amino::Array<Bifrost::Math::float4x4> &bif );
	void Convert ( BIF_O Amino::Array<Bifrost::Math::float4x4> &bif,BIF_I std::vector<Eigen::Matrix4f> &eigen );
}