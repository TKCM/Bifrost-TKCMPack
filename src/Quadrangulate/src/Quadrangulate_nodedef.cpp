#include "Quadrangulate_nodedef.h"

void TKCM::quadrangulate_core(
	const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& point_position,
	const	Amino::Ptr<Amino::Array<uInt>>& face_vertex,
	const	Amino::Ptr<Amino::Array<uInt>>& face_offset,
	const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_face,
	const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_side,
			Amino::MutablePtr<Amino::Array<bool>>& half_edge_tag_data_to_delete
) {
	half_edge_tag_data_to_delete = Amino::newMutablePtr<Amino::Array<bool>>();

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// 入力データの確認
	if (!face_offset || face_offset->empty () || 
		!point_position || point_position->empty () || 
		!face_vertex || face_vertex->empty () ||
		!face_ver_adjacent_edge_face || face_ver_adjacent_edge_face->empty() ||
		!face_ver_adjacent_edge_side || face_ver_adjacent_edge_side->empty()) 
	{
		return;
	}
	size_t vertexCount = face_vertex->size ();
	size_t polyCount = face_offset->size () - 1;
	if (face_ver_adjacent_edge_face->size() != vertexCount ||
		face_ver_adjacent_edge_side->size() != vertexCount)
	{
		return;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////

	Amino::Array<int> deleteLocID ( polyCount );
	#pragma omp parallel for
	for (size_t polyNum = 0; polyNum < polyCount; ++polyNum) {
		int polySize = face_offset->at ( polyNum + 1 ) - face_offset->at ( polyNum );
		if (polySize != 3) { continue; }

		int firstVertexNum = face_offset->at ( polyNum );
		const Bifrost::Math::float3& p0 = point_position->at ( face_vertex->at ( firstVertexNum ) );
		const Bifrost::Math::float3& p1 = point_position->at ( face_vertex->at ( firstVertexNum + 1 ) );
		const Bifrost::Math::float3& p2 = point_position->at ( face_vertex->at ( firstVertexNum + 2 ) );

		float len0 = TKCM::DistanceSquared(p0, p1);
		float len1 = TKCM::DistanceSquared(p1, p2);
		float len2 = TKCM::DistanceSquared(p2, p0);

		float len = len0;
		len = std::max ( len, len1 );
		len = std::max ( len, len2 );
		int localID = 0;
		if (TKCM::AlmostEqual ( len, len2 )) { localID = 2; }
		if (TKCM::AlmostEqual ( len, len1 )) { localID = 1; }

		deleteLocID[polyNum] = localID;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////

	half_edge_tag_data_to_delete->resize(vertexCount, false);
	#pragma omp parallel for
	for (size_t polyNum = 0; polyNum < polyCount; ++polyNum) {
		int edgeNum = face_offset->at ( polyNum ) + deleteLocID[polyNum];
		size_t fEdgeFace = face_ver_adjacent_edge_face->at(edgeNum);
		size_t fEdgeSide = face_ver_adjacent_edge_side->at(edgeNum);
		if (polyCount <= fEdgeFace) { continue; }

		if (deleteLocID[fEdgeFace] == fEdgeSide) {
			half_edge_tag_data_to_delete->at(edgeNum) = true;
		}
	}
}
