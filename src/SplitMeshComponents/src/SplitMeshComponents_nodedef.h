#ifndef BIF_SPLITMESHCOMPONENTS_H
#define BIF_SPLITMESHCOMPONENTS_H
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
	void split_points_core(
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
		"name=Modeling_Toolbox::Internal::split_points_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, false}] "
	);

	/*TKCM_DECL
	void SplitEdges(
				Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position	AMINO_ANNOTATE("Amino::InOut outName=new_point_position"),
				Amino::Ptr<Amino::Array<uInt>>& source_face_vertex						AMINO_ANNOTATE("Amino::InOut outName=new_face_vertex"),
				Amino::Ptr<Amino::Array<uInt>>& source_face_offset						AMINO_ANNOTATE("Amino::InOut outName=new_face_offset"),
		const	Amino::Ptr<Amino::Array<uInt>>& half_edge_start_vertex_id,
		const	Amino::Ptr<Amino::Array<float>>& split_ratio,
		const	bool& apply_to_adjacent_half_edge,
				bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Modeling_Toolbox::Internal::split_edge_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, false}] "
	);*/
}


/*
const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
	const	Amino::Ptr<Amino::Array<bool>>& point_tag_data_to_split,
			Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
			Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
			Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
			Amino::MutablePtr<Amino::Array<uInt>>& point_dst_to_source,
			Amino::MutablePtr<Amino::Array<uInt>>& face_dst_to_source,
			Amino::MutablePtr<Amino::Array<uInt>>& face_vertex_dst_to_source,
* 
TKCM_DECL
void SplitEdges (
	BIF_I Amino::Ptr<Bifrost::Object>& mesh					AMINO_ANNOTATE ( "Amino::Port name=mesh" ),
	BIF_I Amino::Ptr<Amino::Array<UInt32>>& polygonIDs		AMINO_ANNOTATE ( "Amino::Port name=polygon_IDs" ),
	BIF_I Amino::Ptr<Amino::Array<UInt32>>& pointIDs		AMINO_ANNOTATE ( "Amino::Port name=point_IDs" ),
	BIF_I Amino::Ptr<Amino::Array<float>>& ratio			AMINO_ANNOTATE ( "Amino::Port name=split_ratio" ),
	BIF_I bool applyToAdjacentHalfEdge						AMINO_ANNOTATE ( "Amino::Port name=apply_to_adjacent_half_edge" ),
	BIF_I bool& transferNormalProp							AMINO_ANNOTATE ( "Amino::Port name=transfer_prop_normal" ),
	BIF_I bool& transferUVsProp								AMINO_ANNOTATE ( "Amino::Port name=transfer_prop_uvs" ),
	BIF_I bool& transferOtherProp							AMINO_ANNOTATE ( "Amino::Port name=transfer_prop_other" ),
	BIF_O Amino::Ptr<Bifrost::Object>& result				AMINO_ANNOTATE ( "Amino::Port name=result" )
) AMINO_ANNOTATE (
	"Amino::Node "
	"name=TKCM::Modeling_Toolbox::Internal::split_edge_core "
	"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
);

TKCM_DECL
void SplitPolygons (
	BIF_I Amino::Ptr<Bifrost::Object>& mesh					AMINO_ANNOTATE ( "Amino::Port name=mesh" ),
	BIF_I Amino::Ptr<Amino::Array<UInt32>>& polygonIDs		AMINO_ANNOTATE ( "Amino::Port name=polygon_IDs" ),
	BIF_I Amino::Ptr<Amino::Array<UInt32>>& pointIDs1		AMINO_ANNOTATE ( "Amino::Port name=point_IDs_1" ),
	BIF_I Amino::Ptr<Amino::Array<float>>& ratio1			AMINO_ANNOTATE ( "Amino::Port name=split_ratio_1" ),
	BIF_I Amino::Ptr<Amino::Array<UInt32>>& pointIDs2		AMINO_ANNOTATE ( "Amino::Port name=point_IDs_2" ),
	BIF_I Amino::Ptr<Amino::Array<float>>& ratio2			AMINO_ANNOTATE ( "Amino::Port name=split_ratio_2" ),
	BIF_I Amino::Ptr<Amino::Array<Amino::Ptr<Amino::Array<Bifrost::Math::float3>>>>& addPoints	AMINO_ANNOTATE ( "Amino::Port name=add_points" ),
	BIF_I bool applyToAdjacentPolygon						AMINO_ANNOTATE ( "Amino::Port name=apply_to_adjacent_polygon" ),
	BIF_I bool& transferNormalProp							AMINO_ANNOTATE ( "Amino::Port name=transfer_prop_normal" ),
	BIF_I bool& transferUVsProp								AMINO_ANNOTATE ( "Amino::Port name=transfer_prop_uvs" ),
	BIF_I bool& transferOtherProp							AMINO_ANNOTATE ( "Amino::Port name=transfer_prop_other" ),
	BIF_O Amino::Ptr<Bifrost::Object>& result					AMINO_ANNOTATE ( "Amino::Port name=result" )
) AMINO_ANNOTATE ( 
	"Amino::Node "
	"name=TKCM::Modeling_Toolbox::Internal::split_polygon_core "
	"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] " 
);
*/
#endif // BIF_SPLITMESHCOMPONENTS_H
