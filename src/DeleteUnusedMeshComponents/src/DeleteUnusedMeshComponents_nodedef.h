#ifndef BIF_DELETEUNUSEDMESHCOMPONENTS_H
#define BIF_DELETEUNUSEDMESHCOMPONENTS_H
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

namespace TKCM {
	TKCM_DECL
	void delete_unused_mesh_components_core(
		const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
		const	bool& unused_point,
		const	bool& inline_point,
		const	float& line_threshold,
		const	bool& degenerate_polygon,
		const	bool& overlap_polygon,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
				Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
				Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
				Amino::MutablePtr<Amino::Array<uInt>>& point_dst_to_source,
				Amino::MutablePtr<Amino::Array<uInt>>& face_dst_to_source,
				Amino::MutablePtr<Amino::Array<uInt>>& face_vertex_dst_to_source,
				bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Modeling_Toolbox::Internal::delete_unused_mesh_components_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);
}

#endif // BIF_DELETEUNUSEDMESHCOMPONENTS_H