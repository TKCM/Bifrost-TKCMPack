#ifndef BIF_SPLITEDGE_H
#define BIF_SPLITEDGE_H
// Bifrost Amino
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>
#include <Amino/Core/Array.h>
#include <Amino/Core/String.h>

// TKCM
#include "../../core/TKCM_core.h"
#include "../../core/TKCM_fn.h"
#include "../../core/TKCM_math.h"
#include "../../core/TKCM_geometry.h"
#include "../../core/TKCM_geometry.cpp"

namespace TKCM {
	TKCM_DECL
	void split_edge_core(
				Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position	AMINO_ANNOTATE("Amino::InOut outName=new_point_position"),
				Amino::Ptr<Amino::Array<uInt>>& source_face_vertex						AMINO_ANNOTATE("Amino::InOut outName=new_face_vertex"),
				Amino::Ptr<Amino::Array<uInt>>& source_face_offset						AMINO_ANNOTATE("Amino::InOut outName=new_face_offset"),
		const	Amino::Ptr<Amino::Array<uInt>>& face_vertex_adjacent_edge_face,
		const	Amino::Ptr<Amino::Array<uInt>>& face_vertex_adjacent_edge_side,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_edge_face,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_edge_side,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacency_index,
		const	Amino::Ptr<Amino::Array<uInt>>& half_edge_id,
		const	Amino::Ptr<Amino::Array<float>>& split_ratio,
		const	bool& apply_to_adjacent_half_edge,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::uint4>>& point_dst_to_source,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::float4>>& point_weight_dst_to_source,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::uint4>>& face_vertex_dst_to_source,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::float4>>& face_vertex_weight_dst_to_source,
				Amino::MutablePtr<Amino::Array<bool>>& new_point_tag_data,
				bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Modeling_Toolbox::Internal::split_edge_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);
}

#endif // BIF_SPLITEDGE_H
