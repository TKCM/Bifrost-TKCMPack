#ifndef BIF_DETACHCOMPONENTS_H
#define BIF_DETACHCOMPONENTS_H
// Bifrost Amino
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>
#include <Amino/Core/String.h>

// TKCM
#include "../../core/TKCM_core.h"
#include "../../core/TKCM_math.h"
#include "../../core/TKCM_geometry.h"
#include "../../core/TKCM_geometry.cpp"

namespace TKCM {
	enum SplitComponentTagType : uInt
	{
		face_vertex,
		point,
		half_edge,
		face
	};
	

	TKCM_DECL
	void detach_mesh_components_core(
				Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position	AMINO_ANNOTATE("Amino::InOut outName=new_point_position"),
				Amino::Ptr<Amino::Array<uInt>>& source_face_vertex						AMINO_ANNOTATE("Amino::InOut outName=new_face_vertex"),
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_edge_face,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_edge_side,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacency_index,
		const	uInt& component_tag_data_type,
		const	Amino::Ptr<Amino::Array<bool>>& component_tag_data_to_split,
				Amino::MutablePtr<Amino::Array<uInt>>& point_dst_to_source,
				Amino::MutablePtr<Amino::Array<bool>>& new_point_tag_data,
				bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Modeling_Toolbox::Internal::detach_mesh_components_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);

}

#endif // BIF_DETACHCOMPONENTS_H
