#ifndef BIF_EXPAND_H
#define BIF_EXPAND_H
// Bifrost Amino
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>

// TKCM
#include "../../core/TKCM_core.h"
#include "../../core/TKCM_math.h"
#include "../../core/TKCM_geometry.h"
#include "../../core/TKCM_geometry.cpp"

namespace TKCM {
	enum ExpandNormalType : uInt
	{
		Plane_Normal,
		Face_Normal,
		Connected_Edge_Direction
	};

	Bifrost::Math::uint4 ToVertexIndies(
		const uInt& faceID, 
		const Bifrost::Math::uint4& pointID, 
		const Amino::Ptr<Amino::Array<uInt>>& face_offset,
		const Amino::Ptr<Amino::Array<uInt>>& source_face_vertex ){
		Bifrost::Math::uint4 result;
		uInt p = pointID.x;
		for (int i = face_offset->at(faceID); i < face_offset->at(faceID+1); ++i){
			if (source_face_vertex->at(i) == p){
				result.x = i;
			}
		}
		p = pointID.y;
		for (int i = face_offset->at(faceID); i < face_offset->at(faceID + 1); ++i){
			if (source_face_vertex->at(i) == p){
				result.y = i;
			}
		}
		p = pointID.z;
		for (int i = face_offset->at(faceID); i < face_offset->at(faceID + 1); ++i){
			if (source_face_vertex->at(i) == p){
				result.z = i;
			}
		}
		p = pointID.w;
		for (int i = face_offset->at(faceID); i < face_offset->at(faceID + 1); ++i){
			if (source_face_vertex->at(i) == p){
				result.w = i;
			}
		}
		return result;
	}

	TKCM_DECL
	void poly_edge_expand_core(
		const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
		const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_face_normal,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_face,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_side,
		const	Amino::Ptr<Amino::Array<uInt>>& half_edge_to_expand,
		const	bool& apply_to_adjacent_half_edge,
		const	float& offset,
		const	uInt& division,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
				Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
				Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::uint4>>&	point_dst_to_source,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::float4>>& point_weight_to_source,
				Amino::MutablePtr<Amino::Array<uInt>>& face_dst_to_source,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::uint4>>&	face_vertex_dst_to_source,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::float4>>& face_vertex_weight_to_source,
				Amino::MutablePtr<Amino::Array<bool>>& new_face_tag_data,
				bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Modeling_Toolbox::Internal::poly_edge_expand_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, false}] "
	);
}

#endif // BIF_EXPAND_H
