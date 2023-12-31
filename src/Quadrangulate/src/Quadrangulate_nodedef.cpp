#include "Quadrangulate_nodedef.h"

void Quadrangulate (
	BIF_I Amino::Ptr < Bifrost::Object>& mesh,
	BIF_O Amino::Ptr<Amino::Array<bool>>& deleteEdge
) {
	/////////////////////////////////////////////////////////////////////////////////////////////////////

	auto polySizePtr = BifrostExp::getDataGeoPropValues<BifGeoIndex> ( *mesh, "face_offset" );
	auto polyPoiListPtr = BifrostExp::getDataGeoPropValues<BifGeoIndex> ( *mesh, "face_vertex" );
	auto poiPosPtr = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3> ( *mesh, "point_position" );
	auto faceAdjPtr = BifrostExp::getDataGeoPropValues<BifFaceEdge> ( *mesh, "face_vertex_adjacent_edge" );
	if (!polySizePtr || polySizePtr->empty () || !poiPosPtr || poiPosPtr->empty () || !faceAdjPtr || faceAdjPtr->empty ()) {
		Amino::Array<bool> result ( 0 );
		deleteEdge = Amino::newClassPtr<Amino::Array<bool>> ( std::move(result) );
		return;
	}
	UInt32 vertexCount = UInt32 (faceAdjPtr->size ());
	UInt32 polyCount = UInt32 (polySizePtr->size () - 1);

	/////////////////////////////////////////////////////////////////////////////////////////////////////

	Amino::Array<int> deleteLocID ( polyCount );
	BIF::parallel_for ( tbb::blocked_range<size_t> ( 0, polyCount ), [&]( const tbb::blocked_range<size_t>& range ) {
		for (size_t polyNum = range.begin (); polyNum != range.end (); ++polyNum) {
			int polySize = polySizePtr->at ( polyNum + 1 ) - polySizePtr->at ( polyNum );
			if (polySize != 3) { continue; }

			int firstVertexNum = polySizePtr->at ( polyNum );
			const Bifrost::Math::float3& p0 = poiPosPtr->at ( polyPoiListPtr->at ( firstVertexNum ) );
			const Bifrost::Math::float3& p1 = poiPosPtr->at ( polyPoiListPtr->at ( firstVertexNum + 1 ) );
			const Bifrost::Math::float3& p2 = poiPosPtr->at ( polyPoiListPtr->at ( firstVertexNum + 2 ) );

			float len0 = Bifrost::Math::mag2 ( p0 - p1 );
			float len1 = Bifrost::Math::mag2 ( p1 - p2 );
			float len2 = Bifrost::Math::mag2 ( p2 - p0 );

			float len = len0;
			len = std::max ( len, len1 );
			len = std::max ( len, len2 );
			int localID;
			if (TKCM::Math::AlmostEqual ( len, len2 )) { localID = 2; }
			if (TKCM::Math::AlmostEqual ( len, len1 )) { localID = 1; }
			if (TKCM::Math::AlmostEqual ( len, len0 )) { localID = 0; }

			deleteLocID[polyNum] = localID;
		}
	} );

	/////////////////////////////////////////////////////////////////////////////////////////////////////

	Amino::Array<bool> result ( vertexCount, false );
	BIF::parallel_for ( tbb::blocked_range<size_t> ( 0, polyCount ), [&]( const tbb::blocked_range<size_t>& range ) {
		for (size_t polyNum = range.begin (); polyNum != range.end (); ++polyNum) {
			int edgeNum = polySizePtr->at ( polyNum ) + deleteLocID[polyNum];
			const BifFaceEdge& fEdge = faceAdjPtr->at ( edgeNum );
			if (polyCount <= fEdge.face) { continue; }

			if (deleteLocID[fEdge.face] == fEdge.side) {
				result[edgeNum] = true;
			}
		}
	} );

	// Deleteノードに渡すエッジ番号のリストを出力する
	deleteEdge = Amino::newClassPtr<Amino::Array<bool>> ( std::move(result) );
}
