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

	TKCM_DECL
	void dissolve_point(
		const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_face,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_side,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_face,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_side,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_index,
		const	Amino::Ptr<Amino::Array<Amino::long_t>>& point_to_dissolve,
//		const	bool& use_tag_data,
//		const	Amino::Ptr<Amino::Array<bool>>& tag_data_to_dissolve,
		const	bool& unusedPoint,
		const	bool& inlinePoint,
		const	bool& degenerate,
		const	bool& overlapPolygon,
		const	bool& transferNormalProp,
		const	bool& transferUVsProp,
		const	bool& transferOtherProp,
		const	bool& unusedComponentTag,
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
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, false}] "
	);
}

#endif // BIF_DISSOLVE_H
