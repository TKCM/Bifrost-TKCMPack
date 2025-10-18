#if 1930 <= _MSC_VER // Visual Studio 2022 (MSVC 14.30 以降)
#define _ALLOW_COMPILER_AND_STL_VERSION_MISMATCH // Clangのバージョンエラーが出るのを無視する
#endif

#ifndef BIF_TAGS_H
#define BIF_TAGS_H

#include <vector>
#include <algorithm>

// Bifrost Amino
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>
#include <Amino/Core/Array.h>
#include <Amino/Core/String.h>

// TKCM
#include "../../core/TKCM_core.h"
#include "../../core/TKCM_geometry.h"
#include "../../core/TKCM_geometry.cpp"

namespace TKCM
{
	/////////////////////////////////////////////////////////////////////
	enum ComponentType : uInt
	{
		face_vertex_component = 0,
		point_component = 1,
		face_component = 2,
		half_edge_vertex_component = 3,
		non
	};

	enum LogicalOperation : uInt
	{
		AND = 0,
		OR = 1
	};

	/////////////////////////////////////////////////////////////////////

	TKCM::ComponentType ToComponentType(const Amino::String& component_type);

	bool NewTagPtr(
				Amino::MutablePtr<Amino::Array<bool>>& newData,
		const	uInt& point_count,
		const	Amino::Ptr<Amino::Array<uInt>>& face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& face_offset,
		const	Amino::String& component_type
	);

	/////////////////////////////////////////////////////////////////////
	TKCM_DECL
	void point_tag_promote(
		const	uInt& point_count,
		const	Amino::Ptr<Amino::Array<uInt>>& face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& face_offset,
		const	Amino::Ptr<Amino::Array<bool>>& in_tag_data,
		const	Amino::String& out_component_type,
		const	uInt& logical_operation,
				Amino::MutablePtr<Amino::Array<bool>>& out_tag_data,
				bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Geometry::Tags::Internal::point_tag_promote "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);

	TKCM_DECL
	void vertex_tag_promote(
		const	uInt& point_count,
		const	Amino::Ptr<Amino::Array<uInt>>& face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& face_offset,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_face,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_side,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_index,
		const	Amino::Ptr<Amino::Array<bool>>& in_tag_data,
		const	Amino::String& out_component_type,
		const	uInt& logical_operation,
		Amino::MutablePtr<Amino::Array<bool>>& out_tag_data,
		bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Geometry::Tags::Internal::vertex_tag_promote "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);

	TKCM_DECL
	void face_tag_promote(
		const	uInt& point_count,
		const	Amino::Ptr<Amino::Array<uInt>>& face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& face_offset,
		const	Amino::Ptr<Amino::Array<bool>>& in_tag_data,
		const	Amino::String& out_component_type,
				Amino::MutablePtr<Amino::Array<bool>>& out_tag_data,
				bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Geometry::Tags::Internal::face_tag_promote "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);

	TKCM_DECL
	void to_face_vertex(
				Amino::Ptr<Amino::Array<bool>>& tag_data AMINO_ANNOTATE("Amino::InOut inName=half_edge_vertex_tag_data outName=face_vertex_tag_data"),
		const	Amino::Ptr<Amino::Array<uInt>>& face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& face_offset,
				bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Geometry::Tags::Internal::to_face_vertex "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);

	TKCM_DECL
	void tag_mesh_border(
		const	uInt& point_count,
		const	Amino::Ptr<Amino::Array<uInt>>& face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& face_offset,
		const	Amino::Ptr<Amino::Array<uInt>>& adjacent_face,
		const	Amino::String& out_component_type,
				Amino::MutablePtr<Amino::Array<bool>>& tag_data,
				bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Geometry::Tags::Internal::tag_mesh_border_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);
}

#endif // BIF_TAGS_H
