#ifndef BIF_DISSOLVE_H
#define BIF_DISSOLVE_H
// Bifrost Amino
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>
#include <Amino/Core/String.h>

// TKCM
#include "../../core/TKCM_core.h"
#include "../../core/TKCM_fn.h"
//#include "functions.h"

namespace TKCM {
	enum ComponentType : int
	{
		Point,
		Half_Edge,
		Face
	};

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// ã§í ä÷êî
	void dissolve_core(
		const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_face,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_side,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_face,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_side,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_index,
		const	TKCM::ComponentType& component_id_type,
		const	Amino::Ptr<Amino::Array<Amino::long_t>>& component_id_to_dissolve,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
				Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
				Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
				Amino::MutablePtr<Amino::Array<uInt>>& point_dst_to_source,
				Amino::MutablePtr<Amino::Array<uInt>>& face_dst_to_source,
				Amino::MutablePtr<Amino::Array<uInt>>& face_vertex_dst_to_source,
				bool& success
	);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Bifrost nodedef
	TKCM_DECL
	void dissolve_point_core(
		const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_face,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_side,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_face,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_side,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_index,
		const	Amino::Ptr<Amino::Array<Amino::long_t>>& point_to_dissolve,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
				Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
				Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
				Amino::MutablePtr<Amino::Array<uInt>>& point_dst_to_source,
				Amino::MutablePtr<Amino::Array<uInt>>& face_dst_to_source,
				Amino::MutablePtr<Amino::Array<uInt>>& face_vertex_dst_to_source,
				bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Modeling_Toolbox::Internal::dissolve_point_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);

	TKCM_DECL
	void dissolve_edge_core(
		const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_face,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_side,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_face,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_side,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_index,
		const	Amino::Ptr<Amino::Array<Amino::long_t>>& half_edge_to_dissolve,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
				Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
				Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
				Amino::MutablePtr<Amino::Array<uInt>>& point_dst_to_source,
				Amino::MutablePtr<Amino::Array<uInt>>& face_dst_to_source,
				Amino::MutablePtr<Amino::Array<uInt>>& face_vertex_dst_to_source,
				bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Modeling_Toolbox::Internal::dissolve_edge_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);

	TKCM_DECL
	void dissolve_face_core(
		const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_face,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_side,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_face,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_side,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_index,
		const	Amino::Ptr<Amino::Array<Amino::long_t>>& face_to_dissolve,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
				Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
				Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
				Amino::MutablePtr<Amino::Array<uInt>>& point_dst_to_source,
				Amino::MutablePtr<Amino::Array<uInt>>& face_dst_to_source,
				Amino::MutablePtr<Amino::Array<uInt>>& face_vertex_dst_to_source,
				bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Modeling_Toolbox::Internal::dissolve_face_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);
}

#endif // BIF_DISSOLVE_H