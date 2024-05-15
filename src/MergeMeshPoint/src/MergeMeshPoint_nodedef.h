#ifndef BIF_MERGEMESHPOINT_H
#define BIF_MERGEMESHPOINT_H

// Bifrost Amino
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>
#include <Amino/Core/Array.h>
#include <Amino/Core/String.h>

// TKCM
#include "../../core/TKCM_core.h"

TKCM_DECL
void MergeMeshPoint (
	const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
	const	bool& use_tag_data,
	const	Amino::Ptr<Amino::Array<bool>>& point_tag_data_to_merge,
	const	bool& invert_tag,
	const	Amino::Ptr<Amino::Array<Amino::Ptr<Amino::Array<Amino::long_t>>>>& target_point_indices,
			Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
			Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
			Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
			Amino::MutablePtr<Amino::Array<uInt>>& point_dst_to_source,
			Amino::MutablePtr<Amino::Array<uInt>>& face_dst_to_source,
			Amino::MutablePtr<Amino::Array<uInt>>& face_vertex_dst_to_source,
			bool& success
) AMINO_ANNOTATE (
	"Amino::Node "
	"name=TKCM::Modeling_Toolbox::Internal::merge_mesh_point_core "
	"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] " 
);

#endif // BIF_MERGEMESHPOINT_H
