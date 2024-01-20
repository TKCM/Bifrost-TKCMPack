#include "SplitMeshComponents_nodedef.h"

void TKCM::split_points_core(
			Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& io_point_position,
			Amino::Ptr<Amino::Array<uInt>>& io_face_vertex,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_edge_face,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_edge_side,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacency_index,
	const	uInt& component_tag_data_type,
	const	Amino::Ptr<Amino::Array<bool>>& component_tag_data_to_split,
			Amino::MutablePtr<Amino::Array<uInt>>& point_dst_to_source,
			Amino::MutablePtr<Amino::Array<bool>>& new_point_tag_data,
			bool& success
){
	// 出力データの準備
	point_dst_to_source = Amino::newMutablePtr<Amino::Array<uInt>>();
	new_point_tag_data = Amino::newMutablePtr<Amino::Array<bool>>();
	success = false;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 入力データのチェック

	if (!io_point_position || io_point_position->empty() ||
		!io_face_vertex || io_face_vertex->empty() ||
		!source_face_offset || source_face_offset->empty() || 
		!component_tag_data_to_split || component_tag_data_to_split->empty() ||
		!point_face_adjacency_index || point_face_adjacency_index->empty() ||
		!point_face_adjacent_edge_face || point_face_adjacent_edge_face->empty() ||
		!point_face_adjacent_edge_side || point_face_adjacent_edge_side->empty()){
		return;
	}
	size_t poiCount = io_point_position->size();
	size_t polyCount = source_face_offset->size()-1;
	size_t verCount = io_face_vertex->size();
	if (point_face_adjacent_edge_face->size() != verCount){ return; }
	if (point_face_adjacent_edge_side->size() != verCount){ return; }

	TKCM::SplitComponentTagType comp_type = static_cast<TKCM::SplitComponentTagType>(component_tag_data_type); ;
	switch (comp_type){
		case TKCM::SplitComponentTagType::face_vertex:
			if (component_tag_data_to_split->size() != verCount){ return; }
			break;
		case TKCM::SplitComponentTagType::point:
			if (component_tag_data_to_split->size() != poiCount){ return; }
			break;
		case TKCM::SplitComponentTagType::half_edge:
			if (component_tag_data_to_split->size() != verCount){ return; }
			break;
		case TKCM::SplitComponentTagType::face:
			if (component_tag_data_to_split->size() != polyCount){ return; }
			break;
		default:
			return;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>> point_position_MPtr = io_point_position.toMutable();
	Amino::MutablePtr<Amino::Array<uInt>> face_vertex_MPtr = io_face_vertex.toMutable();

	point_dst_to_source->reserve(poiCount * 4);
	for (size_t i = 0; i < poiCount; ++i){
		point_dst_to_source->push_back(i); // prop
	}

	switch (comp_type){
		case TKCM::SplitComponentTagType::face_vertex: // タグ付けしたバーテックスを別頂点に分離する
		{
			for (size_t verID = 0; verID < verCount; ++verID){
				if (component_tag_data_to_split->at(verID) == false){ continue; }

				uInt poiID = face_vertex_MPtr->at(verID);

				uInt sharedVerCount = point_face_adjacency_index->at(poiID + 1) - point_face_adjacency_index->at(poiID);
				if (sharedVerCount <= 1){ 
					continue; // 頂点を構成するバーテックスの数が１の場合は分割処理を行わずにスキップ
				} else{
					point_position_MPtr->push_back(point_position_MPtr->at(poiID));
					face_vertex_MPtr->at(verID) = point_position_MPtr->size() - 1;

					point_dst_to_source->push_back(poiID); // prop
				}				
			}
			break;
		}

		case TKCM::SplitComponentTagType::point: // タグ付けした頂点を構成する全てのバーテックスを別頂点に分離する
		{
			for (size_t i = 0; i < poiCount; ++i){
				if (component_tag_data_to_split->at(i) == false){ continue; }

				for (uInt j = point_face_adjacency_index->at(i) + 1; j < point_face_adjacency_index->at(i + 1); ++j){
					uInt faceID = point_face_adjacent_edge_face->at(j);
					uInt face_local_verID = point_face_adjacent_edge_side->at(j);
					uInt verID = source_face_offset->at(faceID) + face_local_verID;
					uInt poiID = face_vertex_MPtr->at(verID);

					point_position_MPtr->push_back(point_position_MPtr->at(poiID));
					face_vertex_MPtr->at(verID) = point_position_MPtr->size() - 1;

					point_dst_to_source->push_back(poiID); // prop
				}
			}
			break;
		}

		case TKCM::SplitComponentTagType::half_edge: // タグ付けしたハーフエッジの始点・終点のバーテックスを別頂点に分離する
		{
			for (size_t i = 0; i < verCount; ++i){
				uInt startVerID = i;
				uInt startPoiID = face_vertex_MPtr->at(startVerID);
				if (component_tag_data_to_split->at(i) == false){ continue; }

				////////////////////////////////////////////////////////////////////////////////
				// エッジのスタート側のバーテックスを別頂点に分割する
				uInt sharedVerCount = point_face_adjacency_index->at(startPoiID + 1) - point_face_adjacency_index->at(startPoiID);
				if (sharedVerCount <= 1){
					continue; // 頂点を構成するバーテックスの数が１の場合は分割処理を行わずにスキップ
				} else{
					point_position_MPtr->push_back(point_position_MPtr->at(startPoiID));
					face_vertex_MPtr->at(startVerID) = point_position_MPtr->size() - 1;

					point_dst_to_source->push_back(startPoiID); // prop
				}

				////////////////////////////////////////////////////////////////////////////////
				// エッジのエンド側のバーテックスを別頂点に分割する
				uInt polyID;
				for (uInt j = point_face_adjacency_index->at(startPoiID) + 1; j < point_face_adjacency_index->at(startPoiID + 1); ++j){
					uInt faceID = point_face_adjacent_edge_face->at(j);
					uInt face_local_verID = point_face_adjacent_edge_side->at(j);
					uInt verID = source_face_offset->at(faceID) + face_local_verID;
					if (startVerID == verID){
						polyID = faceID;
						break;
					}
				}

				uInt endVerID = TKCM::GetNextFaceVertexIndex(polyID, startVerID, source_face_offset);
				uInt endPoiID = face_vertex_MPtr->at(endVerID);

				sharedVerCount = point_face_adjacency_index->at(endPoiID + 1) - point_face_adjacency_index->at(endPoiID);
				if (sharedVerCount <= 1){
					continue; // 頂点を構成するバーテックスの数が１の場合は分割処理を行わずにスキップ
				} else{
					point_position_MPtr->push_back(point_position_MPtr->at(endPoiID));
					face_vertex_MPtr->at(endVerID) = point_position_MPtr->size() - 1;

					point_dst_to_source->push_back(endPoiID); // prop
				}
			}
			break;
		}
		case TKCM::SplitComponentTagType::face: // 指定したフェースを構成する全てのバーテックスを別頂点に分離する
		{
			for (size_t polyID = 0; polyID < polyCount; ++polyID){
				if (component_tag_data_to_split->at(polyID) == false){ continue; }

				for (size_t verID = source_face_offset->at(polyID); verID < source_face_offset->at(polyID + 1); ++verID){
					uInt poiID = face_vertex_MPtr->at(verID);
					
					uInt sharedVerCount = point_face_adjacency_index->at(poiID + 1) - point_face_adjacency_index->at(poiID);
					if (sharedVerCount <= 1){
						continue; // 頂点を構成するバーテックスの数が１の場合は分割処理を行わずにスキップ
					} else{
						point_position_MPtr->push_back(point_position_MPtr->at(poiID));
						face_vertex_MPtr->at(verID) = point_position_MPtr->size() - 1;

						point_dst_to_source->push_back(poiID); // prop
					}
				}
			}
			break;
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 出力
	io_point_position = std::move(point_position_MPtr);
	io_face_vertex = std::move(face_vertex_MPtr);

	new_point_tag_data->resize(io_point_position->size(), false);
	for (int i = poiCount; i < io_point_position->size(); ++i){
		new_point_tag_data->at(i) = true;
	}

	point_dst_to_source->shrink_to_fit();
	success = true;
}