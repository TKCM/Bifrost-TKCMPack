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
	// �o�̓f�[�^�̏���
	point_dst_to_source = Amino::newMutablePtr<Amino::Array<Bifrost::Math::uint4>>();
	point_weight_dst_to_source = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float4>>();
	face_vertex_dst_to_source = Amino::newMutablePtr<Amino::Array<Bifrost::Math::uint4>>();
	face_vertex_weight_dst_to_source = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float4>>();
	new_point_tag_data = Amino::newMutablePtr<Amino::Array<bool>>();
	success = false;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// ���̓f�[�^�̃`�F�b�N

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
	
	// �אڃG�b�W�̕����w���ǉ��������X�g���쐬������ 
	// �����ɕ����w�肪�S���������m��ratio��0.0��1.0�̂��̂̓��X�g����Ȃ�
	Amino::Array<uInt> newEdgeIDs;
	Amino::Array<float> newRatios;
	newEdgeIDs.reserve(verCount * 2);
	newRatios.reserve(verCount * 2);
	for (uInt i = 0; i < half_edge_id->size(); i++){
		float ratio;
		ratio = std::min(1.0f, split_ratio->at(i));
		ratio = std::max(0.0f, split_ratio->at(i));
		// ratio��0.0��������1.0�ꍇ�̓X�L�b�v
		if (TKCM::AlmostEqual(ratio, 0.0f) || TKCM::AlmostEqual(ratio, 1.0f)){
			continue;
		}
		// �G�b�WID�����݂��Ȃ��ꍇ�̓X�L�b�v
		if (verCount <= half_edge_id->at(i)){ continue; }

		// new���X�g�ɓ����ݒ�l�����ɑ��݂���ꍇ�̓X�L�b�v
		bool add = true;
		for (uInt j = 0; j < newEdgeIDs.size(); j++){
			if (half_edge_id->at(i) == newEdgeIDs[j] && TKCM::AlmostEqual(ratio, newRatios[j])){
				add = false;
				break;
			}
		}
		if (add == false){ continue; }

		// new���X�g�ɒǉ�
		// �G�b�W�𕡐��񕪊����邱�ƂɂȂ�ꍇ�A���V�I���~���ɂȂ�悤�ɔz��ɓo�^����
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

		// �אڂ���n�[�t�G�b�W�������Ώۂɂ���ꍇ
		if (apply_to_adjacent_half_edge == true){
			uInt adjFaceID = face_vertex_adjacent_edge_face->at(half_edge_id->at(i));
			if (adjFaceID < polyCount){
				// �G�b�W�����E�G�b�W�łȂ��ꍇ�͗אڂ���n�[�t�G�b�W��new���X�g�ɒǉ�
				uInt adjEdgeID = io_face_offset->at(adjFaceID) + face_vertex_adjacent_edge_side->at(half_edge_id->at(i));// verID = edgeID
				// �G�b�W�𕡐��񕪊����邱�ƂɂȂ�ꍇ�A���V�I���~���ɂȂ�悤�ɔz��ɓo�^����
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

	// �ǉ����钸�_�̃��X�g�ƁA���_ID���X�g���쐬����
	uInt cnt = 0;
	bool firstLoop = true;
	Amino::Array<bool> done(addPoiCnt, false);
	Amino::Array<Amino::Array<uInt>> newPointIDs(polyCount);
	Amino::Array<Amino::Array<Bifrost::Math::uint4>> newFaceVerIDs(polyCount); // for prop
	Amino::Array<Amino::Array<float>> propVerWeight(polyCount); // for prop
	uInt addCount = 0;
	while (true){
		for (uInt i = 0; i < addPoiCnt; i++){
			if (firstLoop){ // ����̏���
				// �n�[�t�G�b�W�ō\�����Ă���|���S��ID���擾���A���_ID���X�g�𔲂��o���ĐV�K���X�g�̉��n�����
				// ���Ƀ��X�g�������o����Ă���ꍇ�̓X�L�b�v�i�����G�b�W�ɑ΂��ĕ����񕪊����s����j
				uInt polyID = TKCM::GetEdgeLeftPolygonID(io_face_vertex, io_face_offset, point_face_adjacent_edge_face, point_face_adjacent_edge_side, point_face_adjacency_index, newEdgeIDs[i]);
				uInt localID;
				if (newPointIDs[polyID].size() != 0){
					continue;
				}
				for (uInt j = io_face_offset->at(polyID); j < io_face_offset->at(polyID + 1); ++j){
					// �|�C���g��ǉ�����|���S���̒��_�ԍ����X�g�𕡐����Ă���
					newPointIDs[polyID].push_back(io_face_vertex->at(j));
					// �G�b�W�̎n�_�̒��_���|���S�����̃��[�J��ID�Ƃ��Ď擾���Ă���
					if (newEdgeIDs[i] == j){ localID = j - io_face_offset->at(polyID); }

					// prop
					propId.x = newEdgeIDs[i];
					newFaceVerIDs[polyID].push_back(propId);
					propVerWeight[polyID].push_back(1.0f);
				}

				// �����ʒu
				Bifrost::Math::uint4 poiID;
				uInt nextVer = TKCM::GetNextFaceVertexIndex(polyID, newEdgeIDs[i], io_face_offset);
				poiID.x = io_face_vertex->at(newEdgeIDs[i]);
				poiID.y = io_face_vertex->at(nextVer);
				const Bifrost::Math::float3& fPos = point_position_MPtr->at(poiID.x);
				const Bifrost::Math::float3& ePos = point_position_MPtr->at(poiID.y);
				Bifrost::Math::float3 newPos = TKCM::LerpVec3(fPos, ePos, newRatios[i]);

				// ���ɓo�^���Ă���V�K�|�C���g�ʒu�ƈ�v������̂��������m�F
				int find = -1;
				for (int j = 0; j < addCount; j++){
					if (TKCM::AlmostEqual(point_position_MPtr->at(poiCount+j), newPos)){
						find = j;
						break;
					}
				}
				uInt newPoiID;
				if (find == -1){ // ���X�g���Ɉ�v����|�W�V���������������ꍇ�̓��X�g�ɒǉ�
					newPoiID = poiCount + addCount;
					point_position_MPtr->push_back(newPos);
					new_point_tag_data->push_back(true); // tag
					addCount++;
					
					// prop
					point_dst_to_source->push_back(poiID);
					propWeight.x = newRatios[i];
					propWeight.y = 1.0f - newRatios[i];
					point_weight_dst_to_source->push_back(propWeight);
				} else{ // ���X�g���ɑ��݂����ꍇ��ID���擾
					newPoiID = poiCount + find;
				}

				// �V�K�|�C���gID�𒸓_ID���X�g�ɑ}������
				// �}���ʒu�̓G�b�W�̃X�^�[�g�|�C���gID�Ŕ��f
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
			} else{ // �����G�b�W�𕡐��񕪊����Ă���ꍇ�̏���
				if (done[i] == true){ continue; }

				uInt polyID = TKCM::GetEdgeLeftPolygonID(io_face_vertex, io_face_offset, point_face_adjacent_edge_face, point_face_adjacent_edge_side, point_face_adjacency_index, newEdgeIDs[i]);
				// �G�b�W�̎n�_�̒��_���|���S�����̃��[�J��ID�Ƃ��Ď擾���Ă���
				uInt localID;
				for (uInt j = io_face_offset->at(polyID); j < io_face_offset->at(polyID + 1); ++j){
					if (newEdgeIDs[i] == j){ localID = j - io_face_offset->at(polyID); }
				}
				// �����ʒu
				Bifrost::Math::uint4 poiID;
				uInt nextVer = TKCM::GetNextFaceVertexIndex(polyID, newEdgeIDs[i], io_face_offset);		
				poiID.x = io_face_vertex->at(newEdgeIDs[i]);
				poiID.y = io_face_vertex->at(nextVer);
				const Bifrost::Math::float3& fPos = point_position_MPtr->at(poiID.x);
				const Bifrost::Math::float3& ePos = point_position_MPtr->at(poiID.y);
				Bifrost::Math::float3 newPos = TKCM::LerpVec3(fPos, ePos, newRatios[i]);

				// ���ɓo�^���Ă���V�K�|�C���g�ʒu�ƈ�v������̂��������m�F
				int find = -1;
				for (int j = 0; j < addCount; j++){
					if (TKCM::AlmostEqual(point_position_MPtr->at(poiCount+j), newPos)){
						find = j;
						break;
					}
				}
				uInt newPoiID;
				if (find == -1){ // ���X�g���Ɉ�v����|�W�V���������������ꍇ�̓��X�g�ɒǉ�
					newPoiID = poiCount + addCount;
					point_position_MPtr->push_back(newPos);
					new_point_tag_data->push_back(true); // tag
					addCount++;

					// prop
					point_dst_to_source->push_back(poiID);
					propWeight.x = newRatios[i];
					propWeight.y = 1.0f - newRatios[i];
					point_weight_dst_to_source->push_back(propWeight);
				} else{ // ���X�g���ɑ��݂����ꍇ��ID���擾
					newPoiID = poiCount + uInt(find);
				}

				// �V�K�|�C���gID�𒸓_ID���X�g�ɑ}������
				// �}���ʒu�͑Ώۃ|���S���̃G�b�W�̃X�^�[�g�G���h�̃|�C���g�ւ̃x�N�g���Ŕ��f����
				for (uInt j = 0; j < newPointIDs[polyID].size(); j++){
					uInt sID = newPointIDs[polyID][j];
					uInt eID = j == newPointIDs[polyID].size() - 1 ? newPointIDs[polyID][0] : newPointIDs[polyID][j + 1];
					Bifrost::Math::float3 sPos, ePos;
					sPos = point_position_MPtr->at(sID);
					ePos = point_position_MPtr->at(eID);

					Bifrost::Math::float3 dir1 = TKCM::Normal(TKCM::Sub(ePos, sPos));		// �X�^�[�g����G���h�܂ł̃x�N�g��
					Bifrost::Math::float3 dir2 = TKCM::Normal(TKCM::Sub(ePos, newPos));	// �ǉ��|�C���g����G���h�܂ł̃x�N�g��
					if (TKCM::AlmostSameDirection(dir1, dir2)){ // �����������m�F
						dir1 = TKCM::Mul(dir1, -0.1f);					// �G���h����X�^�[�g�܂ł̃x�N�g��
						dir2 = TKCM::Normal(TKCM::Sub(sPos,newPos));			// �ǉ��|�C���g����X�^�[�g�܂ł̃x�N�g��
						if (TKCM::AlmostSameDirection(dir1, dir2) && newPointIDs[polyID][j + 1] != newPoiID && newPointIDs[polyID][j] != newPoiID){ // �����������m�F
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

			propId.y = -1; // ���p�ϐ������Z�b�g
			propWeight.y = 0.0f; // ���p�ϐ������Z�b�g
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
		if (newPointIDs[i].size() == 0){ // �I���W�i���̒l���p��
			for (uInt j = io_face_offset->at(i); j < io_face_offset->at(i + 1); j++){
				new_face_vertex.push_back(io_face_vertex->at(j));
				
				// prop
				propId.x = io_face_vertex->at(j);
				face_vertex_dst_to_source->push_back(propId);
				propWeight.x = 1.0f;
				point_weight_dst_to_source->push_back(propWeight);
				origVerCounter++;
			}
		} else{ // �V�������X�g�̒l��K�p
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

				propId.y = 0; // ���p�ϐ������Z�b�g
				propWeight.y = 0.0f; // ���p�ϐ������Z�b�g
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