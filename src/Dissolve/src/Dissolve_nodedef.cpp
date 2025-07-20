#include "Dissolve_nodedef.h"

void TKCM::dissolve_point_core(
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
){
	TKCM::dissolve_core(
		source_point_position,
		source_face_vertex,
		source_face_offset,
		face_ver_adjacent_edge_face,
		face_ver_adjacent_edge_side,
		point_face_adjacecy_face,
		point_face_adjacecy_side,
		point_face_adjacecy_index,
		TKCM::ComponentType::Point,
		point_to_dissolve,
		true,
		point_position,
		face_vertex,
		face_offset,
		point_dst_to_source,
		face_dst_to_source,
		face_vertex_dst_to_source,
		success
	);
}

void TKCM::dissolve_edge_core(
	const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
	const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_face,
	const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_side,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_face,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_side,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_index,
	const	Amino::Ptr<Amino::Array<Amino::long_t>>& half_edge_to_dissolve,
	const	bool& delete_adjacency_point,
			Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
			Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
			Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
			Amino::MutablePtr<Amino::Array<uInt>>& point_dst_to_source,
			Amino::MutablePtr<Amino::Array<uInt>>& face_dst_to_source,
			Amino::MutablePtr<Amino::Array<uInt>>& face_vertex_dst_to_source,
			bool& success
){
	TKCM::dissolve_core(
		source_point_position,
		source_face_vertex,
		source_face_offset,
		face_ver_adjacent_edge_face,
		face_ver_adjacent_edge_side,
		point_face_adjacecy_face,
		point_face_adjacecy_side,
		point_face_adjacecy_index,
		TKCM::ComponentType::Half_Edge,
		half_edge_to_dissolve,
		delete_adjacency_point,
		point_position,
		face_vertex,
		face_offset,
		point_dst_to_source,
		face_dst_to_source,
		face_vertex_dst_to_source,
		success
	);
}

void TKCM::dissolve_face_core(
	const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
	const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_face,
	const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_side,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_face,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_side,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_index,
	const	Amino::Ptr<Amino::Array<Amino::long_t>>& face_to_dissolve,
	const	bool& delete_adjacency_point,
			Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
			Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
			Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
			Amino::MutablePtr<Amino::Array<uInt>>& point_dst_to_source,
			Amino::MutablePtr<Amino::Array<uInt>>& face_dst_to_source,
			Amino::MutablePtr<Amino::Array<uInt>>& face_vertex_dst_to_source,
			bool& success
){
	TKCM::dissolve_core(
		source_point_position,
		source_face_vertex,
		source_face_offset,
		face_ver_adjacent_edge_face,
		face_ver_adjacent_edge_side,
		point_face_adjacecy_face,
		point_face_adjacecy_side,
		point_face_adjacecy_index,
		TKCM::ComponentType::Face,
		face_to_dissolve,
		delete_adjacency_point,
		point_position,
		face_vertex,
		face_offset,
		point_dst_to_source,
		face_dst_to_source,
		face_vertex_dst_to_source,
		success
	);
}