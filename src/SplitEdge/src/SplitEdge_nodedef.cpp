#include "SplitEdge_nodedef.h"

void TKCM::split_edge_core(
			Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& io_point_position,
			Amino::Ptr<Amino::Array<uInt>>& io_face_vertex,
			Amino::Ptr<Amino::Array<uInt>>& io_face_offset,
	const	Amino::Ptr<Amino::Array<uInt>>& face_vertex_adjacent_edge_face,
	const	Amino::Ptr<Amino::Array<uInt>>& face_vertex_adjacent_edge_side,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_edge_face,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_edge_side,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacency_index,
	const	Amino::Ptr<Amino::Array<uInt>>& half_edge_id,
	const	Amino::Ptr<Amino::Array<float>>& split_ratio,
	const	bool& apply_to_adjacent_half_edge,
			Amino::MutablePtr<Amino::Array<Bifrost::Math::uint4>>& point_dst_to_source,
			Amino::MutablePtr<Amino::Array<Bifrost::Math::float4>>& point_weight_dst_to_source,
			Amino::MutablePtr<Amino::Array<Bifrost::Math::uint4>>& face_vertex_dst_to_source,
			Amino::MutablePtr<Amino::Array<Bifrost::Math::float4>>& face_vertex_weight_dst_to_source,
			Amino::MutablePtr<Amino::Array<bool>>& new_point_tag_data,
			bool& success
){
	// 出力データの準備
	point_dst_to_source = Amino::newMutablePtr<Amino::Array<Bifrost::Math::uint4>>();
	point_weight_dst_to_source = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float4>>();
	face_vertex_dst_to_source = Amino::newMutablePtr<Amino::Array<Bifrost::Math::uint4>>();
	face_vertex_weight_dst_to_source = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float4>>();
	new_point_tag_data = Amino::newMutablePtr<Amino::Array<bool>>();
	success = false;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 入力データのチェック

	if (!io_point_position || io_point_position->empty() ||
		!io_face_vertex || io_face_vertex->empty() ||
		!io_face_offset || io_face_offset->empty() ||
		!face_vertex_adjacent_edge_face || face_vertex_adjacent_edge_face->empty() ||
		!face_vertex_adjacent_edge_side || face_vertex_adjacent_edge_side->empty() ||
		!point_face_adjacent_edge_face || point_face_adjacent_edge_face->empty() ||
		!point_face_adjacent_edge_side || point_face_adjacent_edge_side->empty() ||
		!point_face_adjacency_index || point_face_adjacency_index->empty() ||
		!half_edge_id || half_edge_id->empty() ||
		!split_ratio || split_ratio->empty()
		){
		return;
	}
	size_t poiCount = io_point_position->size();
	size_t polyCount = io_face_offset->size() - 1;
	size_t verCount = io_face_vertex->size();
	if (half_edge_id->size() != split_ratio->size()){ return; }

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	
	// 隣接エッジの分割指定を追加したリストを作成し直す 
	// 同時に分割指定が全く同じモノやratioが0.0や1.0のものはリストから省く
	Amino::Array<uInt> newEdgeIDs;
	Amino::Array<float> newRatios;
	newEdgeIDs.reserve(verCount * 2);
	newRatios.reserve(verCount * 2);
	for (uInt i = 0; i < half_edge_id->size(); i++){
		float ratio;
		ratio = std::min(1.0f, split_ratio->at(i));
		ratio = std::max(0.0f, split_ratio->at(i));
		// ratioが0.0もしくは1.0場合はスキップ
		if (TKCM::AlmostEqual(ratio, 0.0f) || TKCM::AlmostEqual(ratio, 1.0f)){
			continue;
		}
		// エッジIDが存在しない場合はスキップ
		if (verCount <= half_edge_id->at(i)){ continue; }

		// newリストに同じ設定値が既に存在する場合はスキップ
		bool add = true;
		for (uInt j = 0; j < newEdgeIDs.size(); j++){
			if (half_edge_id->at(i) == newEdgeIDs[j] && TKCM::AlmostEqual(ratio, newRatios[j])){
				add = false;
				break;
			}
		}
		if (add == false){ continue; }

		// newリストに追加
		// エッジを複数回分割することになる場合、レシオが降順になるように配列に登録する
		float restRatio = ratio;
		if (newEdgeIDs.size() == 0){
			newEdgeIDs.push_back(half_edge_id->at(i));
			newRatios.push_back(restRatio);
		} else{
			for (uInt j = 0; j < newEdgeIDs.size(); ++j){
				if (half_edge_id->at(i) == newEdgeIDs[j] && newRatios[j] < restRatio){
					float r = newRatios[j];
					newRatios[j] = restRatio;
					restRatio = r;
				}

				if (j == newEdgeIDs.size() - 1){
					newEdgeIDs.push_back(half_edge_id->at(i));
					newRatios.push_back(restRatio);
					break;
				}
			}
		}

		// 隣接するハーフエッジも分割対象にする場合
		if (apply_to_adjacent_half_edge == true){
			uInt adjFaceID = face_vertex_adjacent_edge_face->at(half_edge_id->at(i));
			if (adjFaceID < polyCount){
				// エッジが境界エッジでない場合は隣接するハーフエッジもnewリストに追加
				uInt adjEdgeID = io_face_offset->at(adjFaceID) + face_vertex_adjacent_edge_side->at(half_edge_id->at(i));// verID = edgeID
				// エッジを複数回分割することになる場合、レシオが降順になるように配列に登録する
				restRatio = 1.0f - ratio;
				for (uInt j = 0; j < newEdgeIDs.size(); ++j){
					if (adjEdgeID == newEdgeIDs[j] && newRatios[j] < restRatio){
						float r = newRatios[j];
						newRatios[j] = restRatio;
						restRatio = r;
					}

					if (j == newEdgeIDs.size() - 1){
						newEdgeIDs.push_back(adjEdgeID);
						newRatios.push_back(restRatio);
						break;
					}
				}
			}
		}
	}
	uInt addPoiCnt = uInt(newEdgeIDs.size());
	if (addPoiCnt == 0){ return; }

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//

	Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>> point_position_MPtr = io_point_position.toMutable();
	point_position_MPtr->reserve(poiCount + addPoiCnt);

	new_point_tag_data->resize(poiCount, false);
	new_point_tag_data->reserve(poiCount + addPoiCnt);

	// prop
	Bifrost::Math::uint4 propId;
	Bifrost::Math::float4 propWeight;
	point_dst_to_source->reserve(poiCount + addPoiCnt);
	point_weight_dst_to_source->reserve(poiCount + addPoiCnt);
	for (int i = 0; i < poiCount; ++i){
		propId.x = i;
		propWeight.x = 1.0f;
		point_dst_to_source->push_back(propId);
		point_weight_dst_to_source->push_back(propWeight);
;	}

	// 追加する頂点のリストと、頂点IDリストを作成する
	uInt cnt = 0;
	bool firstLoop = true;
	Amino::Array<bool> done(addPoiCnt, false);
	Amino::Array<Amino::Array<uInt>> newPointIDs(polyCount);
	Amino::Array<Amino::Array<Bifrost::Math::uint4>> newFaceVerIDs(polyCount); // for prop
	Amino::Array<Amino::Array<float>> propVerWeight(polyCount); // for prop
	uInt addCount = 0;
	while (true){
		for (uInt i = 0; i < addPoiCnt; i++){
			if (firstLoop){ // 初回の処理
				// ハーフエッジで構成しているポリゴンIDを取得し、頂点IDリストを抜き出して新規リストの下地を作る
				// 既にリストが抜き出されている場合はスキップ（同じエッジに対して複数回分割が行われる）
				uInt polyID = TKCM::GetEdgeLeftPolygonID(io_face_vertex, io_face_offset, point_face_adjacent_edge_face, point_face_adjacent_edge_side, point_face_adjacency_index, newEdgeIDs[i]);
				uInt localID;
				if (newPointIDs[polyID].size() != 0){
					continue;
				}
				for (uInt j = io_face_offset->at(polyID); j < io_face_offset->at(polyID + 1); ++j){
					// ポイントを追加するポリゴンの頂点番号リストを複製しておく
					newPointIDs[polyID].push_back(io_face_vertex->at(j));
					// エッジの始点の頂点をポリゴン内のローカルIDとして取得しておく
					if (newEdgeIDs[i] == j){ localID = j - io_face_offset->at(polyID); }

					// prop
					propId.x = newEdgeIDs[i];
					newFaceVerIDs[polyID].push_back(propId);
					propVerWeight[polyID].push_back(1.0f);
				}

				// 分割位置
				Bifrost::Math::uint4 poiID;
				uInt nextVer = TKCM::GetNextFaceVertexIndex(polyID, newEdgeIDs[i], io_face_offset);
				poiID.x = io_face_vertex->at(newEdgeIDs[i]);
				poiID.y = io_face_vertex->at(nextVer);
				const Bifrost::Math::float3& fPos = point_position_MPtr->at(poiID.x);
				const Bifrost::Math::float3& ePos = point_position_MPtr->at(poiID.y);
				Bifrost::Math::float3 newPos = TKCM::LerpVec3(fPos, ePos, newRatios[i]);

				// 既に登録している新規ポイント位置と一致するものが無いか確認
				int find = -1;
				for (int j = 0; j < addCount; j++){
					if (TKCM::AlmostEqual(point_position_MPtr->at(poiCount+j), newPos)){
						find = j;
						break;
					}
				}
				uInt newPoiID;
				if (find == -1){ // リスト内に一致するポジションが無かった場合はリストに追加
					newPoiID = poiCount + addCount;
					point_position_MPtr->push_back(newPos);
					new_point_tag_data->push_back(true); // tag
					addCount++;
					
					// prop
					point_dst_to_source->push_back(poiID);
					propWeight.x = newRatios[i];
					propWeight.y = 1.0f - newRatios[i];
					point_weight_dst_to_source->push_back(propWeight);
				} else{ // リスト内に存在した場合はIDを取得
					newPoiID = poiCount + find;
				}

				// 新規ポイントIDを頂点IDリストに挿入する
				// 挿入位置はエッジのスタートポイントIDで判断
				for (uInt j = 0; j < newPointIDs[polyID].size(); j++){
					if (newPointIDs[polyID][j] == poiID.x){
						newPointIDs[polyID].insert(newPointIDs[polyID].begin() + j + 1, newPoiID);

						// prop
						propId.x = newEdgeIDs[i];
						propId.y = TKCM::GetNextFaceVertexIndex(polyID, newEdgeIDs[i], io_face_offset);
						newFaceVerIDs[polyID].insert(newFaceVerIDs[polyID].begin() + j + 1, propId);
						propVerWeight[polyID].insert(propVerWeight[polyID].begin() + j + 1, newRatios[i]);
					}
				}

				done[i] = true;
				cnt++;
			} else{ // 同じエッジを複数回分割している場合の処理
				if (done[i] == true){ continue; }

				uInt polyID = TKCM::GetEdgeLeftPolygonID(io_face_vertex, io_face_offset, point_face_adjacent_edge_face, point_face_adjacent_edge_side, point_face_adjacency_index, newEdgeIDs[i]);
				// エッジの始点の頂点をポリゴン内のローカルIDとして取得しておく
				uInt localID;
				for (uInt j = io_face_offset->at(polyID); j < io_face_offset->at(polyID + 1); ++j){
					if (newEdgeIDs[i] == j){ localID = j - io_face_offset->at(polyID); }
				}
				// 分割位置
				Bifrost::Math::uint4 poiID;
				uInt nextVer = TKCM::GetNextFaceVertexIndex(polyID, newEdgeIDs[i], io_face_offset);		
				poiID.x = io_face_vertex->at(newEdgeIDs[i]);
				poiID.y = io_face_vertex->at(nextVer);
				const Bifrost::Math::float3& fPos = point_position_MPtr->at(poiID.x);
				const Bifrost::Math::float3& ePos = point_position_MPtr->at(poiID.y);
				Bifrost::Math::float3 newPos = TKCM::LerpVec3(fPos, ePos, newRatios[i]);

				// 既に登録している新規ポイント位置と一致するものが無いか確認
				int find = -1;
				for (int j = 0; j < addCount; j++){
					if (TKCM::AlmostEqual(point_position_MPtr->at(poiCount+j), newPos)){
						find = j;
						break;
					}
				}
				uInt newPoiID;
				if (find == -1){ // リスト内に一致するポジションが無かった場合はリストに追加
					newPoiID = poiCount + addCount;
					point_position_MPtr->push_back(newPos);
					new_point_tag_data->push_back(true); // tag
					addCount++;

					// prop
					point_dst_to_source->push_back(poiID);
					propWeight.x = newRatios[i];
					propWeight.y = 1.0f - newRatios[i];
					point_weight_dst_to_source->push_back(propWeight);
				} else{ // リスト内に存在した場合はIDを取得
					newPoiID = poiCount + uInt(find);
				}

				// 新規ポイントIDを頂点IDリストに挿入する
				// 挿入位置は対象ポリゴンのエッジのスタートエンドのポイントへのベクトルで判断する
				for (uInt j = 0; j < newPointIDs[polyID].size(); j++){
					uInt sID = newPointIDs[polyID][j];
					uInt eID = j == newPointIDs[polyID].size() - 1 ? newPointIDs[polyID][0] : newPointIDs[polyID][j + 1];
					Bifrost::Math::float3 sPos, ePos;
					sPos = point_position_MPtr->at(sID);
					ePos = point_position_MPtr->at(eID);

					Bifrost::Math::float3 dir1 = TKCM::Normal(TKCM::Sub(ePos, sPos));		// スタートからエンドまでのベクトル
					Bifrost::Math::float3 dir2 = TKCM::Normal(TKCM::Sub(ePos, newPos));	// 追加ポイントからエンドまでのベクトル
					if (TKCM::AlmostSameDirection(dir1, dir2)){ // 同じ向きか確認
						dir1 = TKCM::Mul(dir1, -0.1f);					// エンドからスタートまでのベクトル
						dir2 = TKCM::Normal(TKCM::Sub(sPos,newPos));			// 追加ポイントからスタートまでのベクトル
						if (TKCM::AlmostSameDirection(dir1, dir2) && newPointIDs[polyID][j + 1] != newPoiID && newPointIDs[polyID][j] != newPoiID){ // 同じ向きか確認
							newPointIDs[polyID].insert(newPointIDs[polyID].begin() + j + 1, newPoiID);

							// prop
							propId.x = newEdgeIDs[i];
							propId.y = TKCM::GetNextFaceVertexIndex(polyID, newEdgeIDs[i], io_face_offset);
							newFaceVerIDs[polyID].insert(newFaceVerIDs[polyID].begin() + j + 1, propId);
							propVerWeight[polyID].insert(propVerWeight[polyID].begin() + j + 1, newRatios[i]);
							break;
						}
					}
				}
				done[i] = true;
				cnt++;
			}

			propId.y = -1; // 流用変数をリセット
			propWeight.y = 0.0f; // 流用変数をリセット
		}

		if (cnt < addPoiCnt){
			firstLoop = false;
		} else{
			break;
		}		
	}
	point_position_MPtr->shrink_to_fit();
	new_point_tag_data->shrink_to_fit(); // tag
	point_dst_to_source->shrink_to_fit(); // prop
	point_weight_dst_to_source->shrink_to_fit(); // prop

	Amino::Array<uInt> new_face_offset(io_face_offset->size());
	new_face_offset.at(0) = 0;
	Amino::Array<uInt> new_face_vertex;
	new_face_vertex.reserve(verCount + addPoiCnt * 2);
	face_vertex_dst_to_source->reserve(verCount + addPoiCnt * 2); // prop
	face_vertex_weight_dst_to_source->reserve(verCount + addPoiCnt * 2); // prop
	uInt origVerCounter = 0;
	for (uInt i = 0; i < polyCount; i++){
		uInt polySize = io_face_offset->at(i + 1) - io_face_offset->at(i);
		if (newPointIDs[i].size() == 0){ // オリジナルの値を継続
			for (uInt j = io_face_offset->at(i); j < io_face_offset->at(i + 1); j++){
				new_face_vertex.push_back(io_face_vertex->at(j));
				
				// prop
				propId.x = io_face_vertex->at(j);
				face_vertex_dst_to_source->push_back(propId);
				propWeight.x = 1.0f;
				point_weight_dst_to_source->push_back(propWeight);
				origVerCounter++;
			}
		} else{ // 新しいリストの値を適用
			polySize = newPointIDs[i].size();
			for (uInt j = 0; j < polySize; j++){
				new_face_vertex.push_back(newPointIDs[i][j]);

				// prop
				if (newPointIDs[i][j] < poiCount){
					propId.x = io_face_vertex->at(origVerCounter);
					propWeight.x = 1.0f;
				} else{
					propId = newFaceVerIDs[i][j];
					propWeight.x = propVerWeight[i][j];
					propWeight.y = 1.0f - propWeight.x;
				}
				face_vertex_dst_to_source->push_back(propId);
				point_weight_dst_to_source->push_back(propWeight);

				propId.y = 0; // 流用変数をリセット
				propWeight.y = 0.0f; // 流用変数をリセット
				origVerCounter++;
			}
		}
		new_face_offset.at(i + 1) = new_face_offset.at(i) + polySize;
	}
	new_face_vertex.shrink_to_fit();
	face_vertex_dst_to_source->shrink_to_fit();
	point_weight_dst_to_source->shrink_to_fit();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 

	Amino::MutablePtr<Amino::Array<uInt>> face_offset_MPtr = io_face_offset.toMutable();
	Amino::MutablePtr<Amino::Array<uInt>> face_vertex_MPtr = io_face_vertex.toMutable();
	face_offset_MPtr->swap(new_face_offset);
	face_vertex_MPtr->swap(new_face_vertex);

	io_face_offset = std::move(face_offset_MPtr);
	io_face_vertex = std::move(face_vertex_MPtr);
	io_point_position = std::move(point_position_MPtr);

	success = true;
}