#include "Dissolve_nodedef.h"

void TKCM::dissolve_point(
	const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
	const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_face,
	const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_side,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_face,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_side,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_index,
	const	Amino::Ptr<Amino::Array<Amino::long_t>>& point_to_dissolve,
//	const	bool& use_tag_data,
//	const	Amino::Ptr<Amino::Array<bool>>& tag_data_to_dissolve,
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
){
	// �o�̓f�[�^�̏���
	point_position = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float3>>();
	face_vertex = Amino::newMutablePtr<Amino::Array<uInt>>();
	face_offset = Amino::newMutablePtr<Amino::Array<uInt>>();
	point_dst_to_source = Amino::newMutablePtr<Amino::Array<uInt>>();
	face_dst_to_source = Amino::newMutablePtr<Amino::Array<uInt>>();
	face_vertex_dst_to_source = Amino::newMutablePtr<Amino::Array<uInt>>();
	success = false;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// ���̓f�[�^�̃`�F�b�N

	if (!source_point_position || source_point_position->empty() ||
		!source_face_vertex || source_face_vertex->empty() ||
		!source_face_offset || source_face_offset->empty() ||
		!face_ver_adjacent_edge_face || face_ver_adjacent_edge_face->empty() ||
		!face_ver_adjacent_edge_side || face_ver_adjacent_edge_side->empty() ||
		!point_face_adjacecy_face || point_face_adjacecy_face->empty() ||
		!point_face_adjacecy_side || point_face_adjacecy_side->empty() ||
		!point_face_adjacecy_index || point_face_adjacecy_index->empty() ||
		!point_to_dissolve || point_to_dissolve->empty()
		){
		return;
	}
	size_t poiCount = source_point_position->size();
	size_t polyCount = source_face_offset->size() - 1;
	size_t verCount = source_face_vertex->size();

//	if (use_tag_data){
//		if (!tag_data_to_dissolve || tag_data_to_dissolve->empty()){ return; }
//		if (tag_data_to_dissolve->size() != verCount){ return; }
//	}


	////////////////
//test
	TKCM::ComponentType compType = TKCM::ComponentType::Half_Edge;
	/*
	step-1 : mesh����g�|���W�[�f�[�^���擾����
	step-2 : �R���|�[�l���g�̍폜�����Ŏg�p���郊�X�g����������
	step-3 : �폜�\�胊�X�g���쐬����
	step-4 : �폜�\�胊�X�g������step-2�ŏ����������X�g���X�V����
		Delete = �R���|�[�l���g�v�f���폜�������X�g�ɂ���i���b�V���Ɍ����J���j
		Dissolve = ���̗��̂�ێ�����悤�ɃR���|�[�l���g�v�f���폜�������X�g�ɂ���i�폜�����R���|�[�l���g�����L����|���S���̓}�[�W�j
	step-5 : ���b�V���̃N���[���i�b�v����
		overlapPolygon = �S���������_���X�g�Ő�������|���S�����폜����悤�Ƀ��X�g���X�V����
		degenerate = �ʐς��O�̃|���S�����폜����悤�Ƀ��X�g���X�V����
		inlinePoint = ���̂̍\���Ŏg�p���Ă��Ȃ�������̒��_���폜����悤�Ƀ��X�g���X�V����
		unusedPoint = �|���S���\���Ŏg�p����Ă��Ȃ����_���폜����悤�Ƀ��X�g���X�V����
	step-6 : ���X�g�����Ƀ��b�V���g�|���W�[�f�[�^��V�K�쐬����
	step-7 : ���b�V���f�[�^���쐬����
	step-8 : ���b�V���̎�v�ȃv���p�e�B�\�Q���p������iFaceEdge�⃆�[�U�[�f�[�^�^�C�v�̃v���p�e�B�\�͏Ȃ��j
	step-9 : ���b�V���f�[�^���o�͂���
	*/

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// �R���|�[�l���g�̍폜�����Ŏg�p���郊�X�g����������
	
	// ���_�̍폜�ɂ�钸�_�ԍ��̃Y�����L�^���Ă������X�g (0=�ŏI�o�͂Ɋ܂߂钸�_�Œ��_�ԍ��͂��̂܂܁A1~�ԍ��̂���鐔, -2=�폜�\��)
	Amino::Array<int> poiCondition ( poiCount, 0 );
	// �|���S���t�F�[�X�̏������e���L�^���Ă������X�g�@( false=�폜�����̑ΏۊO(�ŏI�o�͂Ɋ܂߂�j�Atrue=�폜�\��)
	Amino::Array<bool> polyCondition ( polyCount, false );
	// edgeCondition[polyNum][PolyPoiSize] -1=�폜�����̑ΏۊO�@-2=�폜�\��̃G�b�W�{�אڃ|���S���Ȃ� 0~�폜�\��̃G�b�W�{�אڃ|���S��ID
	Amino::Array<Amino::Array<int>> edgeCondition ( polyCount );
	// �|���S�����\�����钸�_�ԍ��̃��X�g(�ҏW���s���₷�����邽�߂�Bifrost�̃��[���Ƃ͈قȂ�Q�����z��̃��X�g)
	Amino::Array<Amino::Array<int>> polyPoiList ( polyCount );
	#pragma omp parallel for
	for (int i = 0; i < polyCount; ++i) {
		int polySize = source_face_offset->at(i+1) - source_face_offset->at(i);
		polyPoiList[i].resize ( polySize );
		for (int j = 0; j < polySize; ++j) {
			int vertexPoiID = source_face_vertex->at ( source_face_offset->at ( i ) + j );
			polyPoiList[i][j] = vertexPoiID;
		}
	}

	// �v���p�e�B�\�]�����L���̏ꍇ�Ɏg�p����ϐ�������
	Amino::Array<bool> propPolyCondition; // ���b�V���̃N���[���A�b�v�������s���O��polyCondition��ێ����Ă������߂̕ϐ�
	Amino::Array<Amino::Array<int>> origPolyPoiList; // ���b�V���̃N���[���A�b�v�������s���O��polyPoiList��ێ����Ă������߂̕ϐ�
	if (transferNormalProp || transferUVsProp || transferOtherProp){
		origPolyPoiList = polyPoiList;
	}
	
	// �폜�\��G�b�W�ŗאڂ���|���S���������true �� �n�ڏ������K�v�ȃ|���S��
	Amino::Array<bool> processPoly ( polyCount, false );

	///////////////////////////////////////////////////////////////////////////////////////////////////// step-3
	
	// �R���|�[�l���g�̍폜���X�g���쐬����
	Amino::Array<bool> deleteVertexID;
	bool main = true;
/*	bool main = false;
	if (tag_data_to_dissolve || tag_data_to_dissolve->empty () == false) {
		switch (compType) {
			case TKCM::ComponentType::Point:		if (tag_data_to_dissolve->size () == poiCount) { main = true; break; }
			case TKCM::ComponentType::Half_Edge:	if (tag_data_to_dissolve->size () == verCount) { main = true; break; }
			case TKCM::ComponentType::Face:			if (tag_data_to_dissolve->size () == polyCount) { main = true; break; }
		}
		if (main) {
			deleteID = *tag_data_to_dissolve;
		}
	}
*/
	///////////////////////////////////////////////////////////////////////////////////////////////////// step-4
	
	// �R���|�[�l���g�^�C�v���Ƃɍ폜�������s��
	if ( main ){
		bool haveAdjProp = true;		
		if (1==0) {
			return;
		} else { // mode == TKCM::Mode::Dissolve
			//////////////////////////////////////////////////// dessolve mode ////////////////////////////////////////////////////
			bool valid = false;
			switch (compType) {
/*				case TKCM::ComponentType::Point: {
					if (haveAdjProp == false || deleteID.size () != poiCount) { break; }

					// restDeleteID��deleteID�̏��𕡐�������AdeleteID�̔z����o�[�e�b�N�X���Ń��Z�b�g����
					Amino::Array<bool> restDeleteID( deleteID.size() );
					restDeleteID = deleteID;
					
					deleteID.assign ( verCount, false );
					
					// ���Z�b�g����deleteID�ɑ΂��A�폜�Ώۂ̒��_���\������o�[�e�b�N�X��true�ɂ���
					#pragma omp parallel for
					for (size_t poiNum = 0; poiNum < poiCount; ++poiNum) {
						if (restDeleteID[poiNum] == false) { continue; }

						// ���_���\������o�[�e�b�N�X��true�ɂ���			
						for (uInt i = point_face_adjacecy_index->at ( poiNum ); i < point_face_adjacecy_index->at ( poiNum + 1 ); ++i) {
							uInt polyNum = point_face_adjacecy_face->at(i);
							uInt localID = point_face_adjacecy_side->at(i);
							uInt vertexID = source_face_offset->at ( polyNum ) + localID;
							deleteID[vertexID] = true;
						}
						valid = true;
					}

					// �|���S�����\������S�Ẵo�[�e�b�N�X���폜�Ώۂ̏ꍇ��polyCondition��true�ɂ���
					#pragma omp parallel for
					for (size_t polyNum = 0; polyNum < polyCount; ++polyNum) {
						bool delPoly = true;
						for (uInt verNum = source_face_offset->at ( polyNum ); verNum < source_face_offset->at ( polyNum + 1 ); ++verNum) {
							// �|���S�����\������S�Ẵo�[�e�b�N�X���폜�\��̏ꍇ�̓|���S�����̂��폜�\��ɂ���
							if (deleteID[verNum] == false) { delPoly = false; break; }
						}
						// ���Ń|���S�����폜�\��ƂȂ����ꍇ�A�|���S�����ɘA������{�[�_�[�G�b�W���������`�F�b�N����
						if (delPoly) {
							bool border = false;
							for (uInt verNum1 = source_face_offset->at ( polyNum ); verNum1 < source_face_offset->at ( polyNum + 1 ) + 1; ++verNum1) {
								// ����+1�Ɓ��̏����́A�Ō�̃G�b�W�̌�����ɂ�����x�����ŏ��̃G�b�W���������邽�߂̑Ή�
								int verNum = verNum1;
								if (verNum1 == source_face_offset->at ( polyNum + 1 )) { verNum = source_face_offset->at ( polyNum ); }

								if (polyCount < face_ver_adjacent_edge_face->at(verNum)) {
									// �{�[�_�[�G�b�W�������ꍇ
									if (border) {
										// 1�O�̃G�b�W���{�[�_�[�������ꍇ�̓|���S���̍폜�\����������ďI��
										delPoly = false;
										break;
									} else {
										// ����̏����̂��߂ɁA���̃G�b�W���{�[�_�[�ł������ƋL���Ă���
										border = true;
									}
								} else {
									// ����̏����̂��߂ɁA���̃G�b�W���{�[�_�[�ł͖��������ƋL���Ă���
									border = false;
								}
							}
						}
						// ���̌����̌��ʂ��󂯂āA�|���S�����폜�\��ɂ���
						if (delPoly) { polyCondition[polyNum] = true; }
					}
					break;
				}
*/
				case TKCM::ComponentType::Half_Edge:
					deleteVertexID.resize(verCount, false);
					for (int i = 0; i < point_to_dissolve->size(); ++i){
						deleteVertexID[point_to_dissolve->at(i)] = true;
					}
					break;

/*				case TKCM::ComponentType::Face: {
					if (haveAdjProp == false || deleteID.size () != polyCount) { break; }

					// restDeleteID��deleteID�̏��𕡐�������AdeleteID�̔z����o�[�e�b�N�X���Ń��Z�b�g����
					Amino::Array<bool> restDeleteID ( deleteID.size () );
					restDeleteID = deleteID;
					
					deleteID.assign ( verCount, false );

					// ���Z�b�g����deleteID�ɑ΂��A�폜�Ώۂ̃|���S�����\������o�[�e�b�N�X��true�ɂ���
					#pragma omp parallel for
					for (size_t polyNum = 0; polyNum < polyCount; ++polyNum) {
						if (restDeleteID[polyNum] == false) { continue; }

						int polySize = source_face_offset->at(polyNum + 1) - source_face_offset->at(polyNum);
						int polyFirstVertex = source_face_offset->at ( polyNum );
						for (int i = 0; i < polySize; ++i) {
							deleteID[polyFirstVertex + i] = true;
						}
						valid = true;
					}
					break;
				}
				*/
			}
			// ���b�V���̃N���[���A�b�v�������L���̏ꍇ�ɂ��̎��_��polyCondition�𕡐����Ă���
			if (transferNormalProp || transferUVsProp || transferOtherProp){
				propPolyCondition.resize(polyCount);
				propPolyCondition = propPolyCondition;
			}

			if(valid = true){
				// �e�|���S���̊e�G�b�W���Ƃɏ�Ԃ��L�^���Ă���
				#pragma omp parallel for
				for (size_t polyNum = 0; polyNum < polyCount; ++polyNum) {
					int polySize = int(polyPoiList[polyNum].size ());
					edgeCondition[polyNum].resize ( polySize, -1 );
					for (int i = 0; i < polySize; ++i) {
						// �|���S�����\�����钸�_�ԍ�
						int poiID = polyPoiList[polyNum][i];
						// �|���S�����̎��_�̒��_�ԍ����擾
						int nextPoiID = TKCM::ArrayValueNextElement<int>( polyPoiList[polyNum], poiID );
						// �n�[�t�G�b�W�ԍ�(=�o�[�e�b�N�X�ԍ�)���擾
						int halfEdgeNum = source_face_offset->at ( polyNum ) + i;

						// �G�b�W�ŗאڂ���|���S���ԍ��ƃn�[�t�G�b�W�ԍ����擾����
						int adjPolyNum = -1;
						int adjHalfEdgeNum;
						if (face_ver_adjacent_edge_face->at(halfEdgeNum) < polyCount) {
							adjPolyNum = face_ver_adjacent_edge_face->at(halfEdgeNum);
							adjHalfEdgeNum = source_face_offset->at ( adjPolyNum ) + face_ver_adjacent_edge_side->at(halfEdgeNum);
						}

						if (deleteVertexID[halfEdgeNum] == true) {
							// �n�[�t�G�b�W���폜�Ώۂ̏ꍇ
							if (adjPolyNum < 0) {
								// �n�[�t�G�b�W�ŗאڂ���|���S���������ꍇ
								edgeCondition[polyNum][i] = -2;
							} else {
								// �n�[�t�G�b�W�ŗאڂ���|���S��������ꍇ
								edgeCondition[polyNum][i] = adjPolyNum;
								processPoly[polyNum] = true;
							}
						} else {
							// �n�[�t�G�b�W���폜�Ώۂł͂Ȃ��ꍇ
							if (adjPolyNum < 0) {
								// �n�[�t�G�b�W�ŗאڂ���|���S���������ꍇ
								edgeCondition[polyNum][i] = -1;
							} else {
								// �n�[�t�G�b�W�ŗאڂ���|���S��������ꍇ
								if (deleteVertexID[adjHalfEdgeNum] == true) {
									// ���Α��̃n�[�t�G�b�W���폜�Ώۂ̏ꍇ
									edgeCondition[polyNum][i] = adjPolyNum;
									processPoly[polyNum] = true;
								} else {
									// ���Α��̃n�[�t�G�b�W���폜�ΏۊO�̏ꍇ
									edgeCondition[polyNum][i] = -1;
								}
							}
						}
					}
				}

				// �|���S���|�C���g���X�g�����Ƀ}�[�W����
				for (uInt polyNum = 0; polyNum < polyCount; ++polyNum) {
					if (polyCondition[polyNum] == true) { continue; } // ���݂̃|���S�����폜�Ώۂ̏ꍇ�̓X�L�b�v
					if (processPoly[polyNum] == false) { continue; } // ���݂̃|���S�����ɍ폜����G�b�W�����݂��Ȃ���΃X�L�b�v

					// ���_�����Ƀ}�[�W���Ă������߂̐V�����|�C���g���X�g
					Amino::Array<int> newPolyPoiList;
					// �}�[�W�������s�����|���S���ԍ����L�^���郊�X�g
					Amino::Array<int> processedPolyNum;

					int startPoiNum = polyPoiList[polyNum][0];
					int thisPoiNum = polyPoiList[polyNum][0];
					int thisPolyNum = polyNum;
					int thisLocalID = 0;
				
					while (true) {
						// �G�b�W�ɗאڂ���|���S�������擾����
						int adjacentPolyNum = edgeCondition[thisPolyNum][thisLocalID];
						if (adjacentPolyNum == -1 || adjacentPolyNum == -2 || polyCondition[adjacentPolyNum] == true) {
							// �G�b�W���폜�Ώۂł͂Ȃ��A�폜����G�b�W���{�[�_�[�G�b�W�̂��ߗאڃ|���S�������݂��Ȃ��A�אڃ|���S���͍폜�\��̂����ꂩ�ꍇ

							if (adjacentPolyNum == -1 || adjacentPolyNum == -2) {
								// �G�b�W���폜�Ώۂł͂Ȃ��A�폜����G�b�W�̗אڃ|���S�������݂��Ȃ��ꍇ
								// ���݂̒��_�ԍ������X�g�Ɋi�[����
								newPolyPoiList.push_back ( thisPoiNum );
							}
							// ���_�ԍ�/���[�J���ԍ������_�֐i�߂�
							thisPoiNum = TKCM::ArrayValueNextElement<int>( polyPoiList[thisPolyNum], thisPoiNum );
							thisLocalID = edgeCondition[thisPolyNum].size () - 1 == thisLocalID ? 0 : thisLocalID + 1;
						}
						else {
							// �אڃ|���S�����}�[�W�Ώۂ̏ꍇ

							// ���݂̃|���S���ԍ����L�^���Ă���
							processedPolyNum.push_back ( thisPolyNum );

							// �|���S���ԍ���אڃ|���S���ԍ��ɍX�V����
							thisPolyNum = adjacentPolyNum;
							// �אڃ|���S�����Ō��݂̒��_�ԍ��̃��[�J��ID���擾����
							for (thisLocalID = 0; thisLocalID < polyPoiList[thisPolyNum].size (); ++thisLocalID) {
								if (polyPoiList[thisPolyNum][thisLocalID] == thisPoiNum) { break; }
							}
						}

						// ���̏����Ώۂ��X�^�[�g�ɖ߂�����I��
						if (thisPolyNum == polyNum && startPoiNum == thisPoiNum) { break; }
						
					}

					if (newPolyPoiList.empty () || newPolyPoiList.size () == 0) { newPolyPoiList.resize ( 1 ); }
					polyPoiList[polyNum].resize ( newPolyPoiList.size () );
					polyPoiList[polyNum] = newPolyPoiList;
					
					// ����}�[�W�������s�����|���S�����폜�\��ɂ��� (�}�[�W�̃x�[�X�ɂ����|���S���̓X�L�b�v����)
					for (int i = 0; i < processedPolyNum.size (); ++i) {
						if (processedPolyNum[i] == polyNum) { continue; }
						polyCondition[processedPolyNum[i]] = true;
					}
				}
				
				// �|���S���|�C���g���X�g�𐮗�
				for (uInt polyNum = 0; polyNum < polyCount; ++polyNum) {
					if (polyCondition[polyNum] == true) { continue; } // �|���S�����폜�����̑Ώۂɐݒ肳��Ă���ꍇ�̓X�L�b�v
					if (polyPoiList[polyNum].size () < 3) { continue; }

					Amino::Array<int> restPoiNum; // �|���S���|�C���g�ԍ��̈ꎞ�ۑ���
					restPoiNum.reserve ( polyPoiList[polyNum].size () );
					for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
						int poiID = polyPoiList[polyNum][i];
						// �폜�\��̒��_�̓J�E���g���Ȃ�
						if (poiCondition[poiID] < 0) { continue; }
						// �������_�ԍ��������g�p�\��̏ꍇ�̓X�L�b�v
						if (2 <= restPoiNum.size () && poiID == restPoiNum[restPoiNum.size () - 2]) {
							restPoiNum.pop_back ();
							continue;
						}
						// ���_�ԍ����ꎞ�ۑ����X�g�ɒǉ�
						restPoiNum.push_back ( poiID );
					}
					// ���_�ԍ����X�g�̍ŏ��ƍŌ�œ������_�ԍ��������g�p����Ă��Ȃ����m�F����
					// ���_�ԍ����X�g��[0,6,4,5,7,6]�̏ꍇ��[6,4,5,7]�ƂȂ�
					bool delArrayVal0 = false;
					for (int i=0; i< restPoiNum.size(); ++i) {
						if (restPoiNum[1 + i] == restPoiNum[restPoiNum.size () - 1 - i]) {
							if (i == 0) { delArrayVal0 = true; }
							restPoiNum.pop_back ();
						} else {
							break;
						}
					}
					if (delArrayVal0) {
						for (int i = 0; i < restPoiNum.size (); ++i) {
							restPoiNum[i] = restPoiNum[i + 1];
						}
						restPoiNum.pop_back ();
					}
					restPoiNum.shrink_to_fit ();
					if (restPoiNum.empty () || restPoiNum.size () == 0) { restPoiNum.resize ( 1 ); }
					polyPoiList[polyNum].resize ( restPoiNum.size () );
					polyPoiList[polyNum] = restPoiNum;
				}
				
				// ���_�͑S�č폜�ΏۊO�ɂ���
				poiCondition.assign ( poiCondition.size (), 0 );
			}
		}
	}


	// ���C���̍폜�����͂����܂�
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// ��������̓��b�V���̃N���[���i�b�v����

	///////////////////////////////////////////////////////////////////////////////////////////////////// step-5
	// �S���������_���X�g�̃|���S�����폜�i�\���̏ꍇ�͍폜���Ȃ��j
/*	if (overlapPolygon){
		// �|���S�����\�����钸�_�̔ԍ������v�������X�g���쐬����
		Amino::Array<int> polyPoiNumSumList ( polyCount );
		BIF::parallel_for ( tbb::blocked_range<size_t> ( 0, polyCount ), [&]( const tbb::blocked_range<size_t>& range ) {
			for (int polyNum = int(range.begin ()); polyNum != range.end (); ++polyNum) {
				if (polyCondition[polyNum] == true) {
					// �ʂ̃|���S���փ}�[�W�\��or�|���S���폜�\��ɐݒ肳��Ă���ꍇ�͕s���l���Z�b�g���Ă���X�L�b�v
					polyPoiNumSumList[polyNum] = -1;
					continue;
				}

				polyPoiNumSumList[polyNum] = 0;
				for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
					polyPoiNumSumList[polyNum] += polyPoiList[polyNum][i];
				}
			}
		} );

		// ���g��菬�����|���S���ԍ���Ώۂɂ��āA���S�ɏd�Ȃ��Ă���|���S����T��
		BIF::parallel_for ( tbb::blocked_range<size_t> ( 0, polyCount ), [&]( const tbb::blocked_range<size_t>& range ) {
			for (int polyNum = int(range.begin ()); polyNum != range.end (); ++polyNum) {
				int polyNSum = polyPoiNumSumList[polyNum]; // �|���S�����\�����钸�_�ԍ��̍��v�l���擾
				if (polyNSum <= 2) { continue; } // �l���s���l�̏ꍇ�̓X�L�b�v

				for (int i = 0; i < polyNum; ++i) {
					// 0�Ԃ��珇��polyNum�������|���S����T�� �������g�܂ŒH�蒅������I��
					if (polyNSum == polyPoiNumSumList[i]) {
						if (polyPoiList[polyNum].size () != polyPoiList[i].size ()) { continue; } // polyNum�������ł��|���S���̉搔���قȂ�ꍇ�̓X�L�b�v

						// �|���S�����\�����钸�_�ԍ��Ƃ��̕��я����S�ē������`�F�b�N����
						bool same = true;
						for (int j = 0; j < polyPoiList[polyNum].size (); ++j) {
							int thisP = polyPoiList[polyNum][j];
							int nextP = TKCM::Math::GetNextVal ( polyPoiList[polyNum], thisP );

							if (nextP != TKCM::Math::GetNextVal ( polyPoiList[i], thisP )) {
								same = false;
								break;
							}
						}
						if (same == false) { continue; } // ���_�ԍ����X�g���قȂ�ꍇ�̓X�L�b�v

						// ���_�ԍ����X�g�����v�����ꍇ
						// ��̏����ł��̃|���S���͍폜�Ώۂł���Ɣ��肳���悤�ɁAmergeTargetFaceId��true���Z�b�g����
						polyCondition[polyNum] = true;
					}
				}
			}
		} );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// �ʐς��O�̃|���S�����폜
	if (degenerate) {
		BIF::parallel_for ( tbb::blocked_range<size_t> ( 0, polyCount ), [&]( const tbb::blocked_range<size_t>& range ) {
			for (int polyNum = int(range.begin ()); polyNum != range.end (); ++polyNum) {
				if (polyCondition[polyNum] == true) { continue; } // �ʂ̃|���S���փ}�[�W�\��or�|���S���폜�\��ɐݒ肳��Ă���ꍇ�̓X�L�b�v

				int polySize = int ( polyPoiList[polyNum].size () ); // �|���S���̉搔���擾
				if (polySize <= 2) { polyCondition[polyNum] = true; continue; }

				const Bifrost::Math::float3& p0 = source_point_position->at ( polyPoiList[polyNum][0] ); // �ŏ��̒��_�̈ʒu���擾
				// �T�u�g���C�A���O���̖ʐς����ɎZ�o���Ă���
				bool zero = true;
				for (int j = 0; j < polySize - 2; ++j) {
					const Bifrost::Math::float3& p1 = source_point_position->at ( polyPoiList[polyNum][j + 1] );
					const Bifrost::Math::float3& p2 = source_point_position->at ( polyPoiList[polyNum][j + 2] );

					float area = TKCM::Math::GetTriangleArea ( p0, p1, p2 ); // �ʐ�
					if (0.0001f < area) {
						// �ʐς��O�ł͂Ȃ��T�u�g���C�A���O�����P�ł�����΃��[�v�I��
						zero = false;
						break;
					}
				}
				if (zero) {
					// zero = true�̂܂܂̏ꍇ�̓[���ʐσ|���S���Ɣ��肵�A�폜�����̑ΏۂƂ��邽�߂�true���Z�b�g����
					polyCondition[polyNum] = true;
				}
			}
			} );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// �s�v�Ȓ��_���폜
	if (inlinePoint || unusedPoint) {
		// ���_�ō\�����Ă���|���S��ID��Z�߂����X�g
		Amino::Array< Amino::Array < int >> poiPolyNum ( poiCount );
		for (uInt polyNum = 0; polyNum < polyCount; ++polyNum) {
			if (polyCondition[polyNum] == true) { continue; } // �|���S�����폜�����̑Ώۂɐݒ肳��Ă���ꍇ�̓X�L�b�v

			// �|���S���Ŏg�p���钸�_�ԍ��𐮗�����
			Amino::Array < int > restPoiIDlist ( 0 );
			for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
				int poiID = polyPoiList[polyNum][i];
				if (poiCondition[poiID] < 0) { continue; } // �폜�\��̒��_�̓��X�g�ɒǉ����Ȃ�
				restPoiIDlist.push_back ( poiID );
			}
			// �|���S���̉搔���s���l�̏ꍇ�̓X�L�b�v
			if (restPoiIDlist.size() < 3) { continue; } 
			// ���_���ƂɎg�p���Ă���|���S���ԍ����L�^���Ă���
			for (int i = 0; i < restPoiIDlist.size(); ++i) {
				int poiID = restPoiIDlist[i]; // �|���S�����\�����钸�_�ԍ�
				poiPolyNum[poiID].push_back ( polyNum ); // ���_�ԍ��̔z��Ƀ|���S���ԍ����L�^����
			}
		}

		if (unusedPoint) {
			BIF::parallel_for ( tbb::blocked_range<size_t> ( 0, poiCount ), [&]( const tbb::blocked_range<size_t>& range ) {
				for (size_t poiID = range.begin (); poiID != range.end (); ++poiID) {
					if (0 == poiPolyNum[poiID].size ()) { poiCondition[poiID] = -2; }
				}
			} );
		}

		if (inlinePoint) {
			BIF::parallel_for ( tbb::blocked_range<size_t> ( 0, poiCount ), [&]( const tbb::blocked_range<size_t>& range ) {
				for (int poiID = int(range.begin ()); poiID != range.end (); ++poiID) {
					if (poiCondition[poiID] < 0) { continue; }
					if (3 <= poiPolyNum[poiID].size ()) { continue; } // �R�ȏ�̃|���S���Ŏg�p����Ă���ꍇ�̓X�L�b�v

					switch (poiPolyNum[poiID].size ()) {
						case 2:
						{
							// ���_�ō\������Ă���Q�̃|���S�����אڂ��Ă��邩�m�F����
							int polyNum0 = poiPolyNum[poiID][0];
							int nextPoiID0 = TKCM::Math::GetNextVal ( polyPoiList, polyNum0, poiID, poiCondition );
							int prevPoiID0 = TKCM::Math::GetPreviousVal ( polyPoiList, polyNum0, poiID, poiCondition );
							if (nextPoiID0 == prevPoiID0 || poiID == nextPoiID0 || poiID == prevPoiID0) { break; }
							int polyNum1 = poiPolyNum[poiID][1];
							int nextPoiID1 = TKCM::Math::GetNextVal ( polyPoiList, polyNum1, poiID, poiCondition );
							int prevPoiID1 = TKCM::Math::GetPreviousVal ( polyPoiList, polyNum1, poiID, poiCondition );
							if (nextPoiID1 == prevPoiID1 || poiID == nextPoiID1 || poiID == prevPoiID1) { break; }

							if (nextPoiID0 != prevPoiID1 && prevPoiID0 != nextPoiID1) { break; }// �אڂ��Ă��Ȃ��ꍇ�̓X�L�b�v
							// �אڂ��Ă���ꍇ��case=1�ɐi��
						}
						case 1:
						{
							// ������̒��_�i=�`����\�����邽�߂Ɏg�p���Ă��Ȃ����_�j���m�F����
							int polyNum = poiPolyNum[poiID][0];
							int nextPoiID = TKCM::Math::GetNextVal ( polyPoiList[polyNum], poiID );
							int prevPoiID = TKCM::Math::GetPreviousVal ( polyPoiList[polyNum], poiID );
							if (nextPoiID == prevPoiID || poiID == nextPoiID || poiID == prevPoiID) { break; }

							const Bifrost::Math::float3& p0 = source_point_position->at ( nextPoiID );
							const Bifrost::Math::float3& p1 = source_point_position->at ( poiID );
							const Bifrost::Math::float3& p2 = source_point_position->at ( prevPoiID );
							Bifrost::Math::float3 p0p1UnitVec3 = Bifrost::Math::normalize ( p0 - p1 );
							Bifrost::Math::float3 p1p2UnitVec3 = Bifrost::Math::normalize ( p1 - p2 );
							float dot = Bifrost::Math::dot ( p0p1UnitVec3, p1p2UnitVec3 );

							if (TKCM::Math::AlmostEqual ( dot, 1.0f )) {
								poiCondition[poiID] = -2;
							}
							break;
						}
					}
				}
			} );
		}
	}
*/
	//////////////////////////////////////////////////////////////////////////////////////////////////////////// step-6
	// ���_�ԍ��̕ύX���X�g���쐬
	int offset = 0;
	for (uInt i = 0; i < poiCount; ++i) {
		if (poiCondition[i] < 0) {
			offset++;
		}
		else {
			poiCondition[i] = offset;
		}
	}

	Amino::Array<uInt> newPolySize; // �|���S���̉搔�̃��X�g (�z��v�f��0�X�^�[�g�̉��Z��)
	Amino::Array<uInt> newPolyPoiID; // �|���S�����\�����钸�_�ԍ����X�g
	Amino::Array<Bifrost::Math::float3> newPoiPos; // ���_�̈ʒu�f�[�^
	face_offset->reserve ( polyCount );
	face_vertex->reserve ( verCount );
	point_position->reserve ( poiCount );

	// �|���S���̉搔�̃��X�g/���_�ԍ����X�g
	face_offset->push_back ( 0 );
	for (uInt polyNum = 0; polyNum < polyCount; ++polyNum) {
		// �|���S�����폜�\��ɐݒ肳��Ă���ꍇ�̓X�L�b�v
		if (polyCondition[polyNum] == true) { 
			polyPoiList[polyNum].resize ( 0 );
			continue; 
		} 

		int polySize = 0;// �|���S���̉搔
		for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
			int poiID = polyPoiList[polyNum][i];

			// �폜�\��̒��_�̓J�E���g���Ȃ�
			if (poiCondition[poiID] < 0) { continue; }

			// �|���S���̉搔���C���N�������g
			polySize++;
		}

		// �|���S�����\�����钸�_�����s�����Ă���ꍇ�̓X�L�b�v
		if (polySize <= 2) { 
			polyPoiList[polyNum].resize ( 0 ); 
			continue; 
		}

		// �|���S���̉搔�̃��X�g�ɒǉ��i���Z�j
		face_offset->push_back (face_offset->back () + polySize );
		// �|�C���g�ԍ��𐮗����ă|���S���|�C���g���X�g�ɒǉ�
		for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
			int poiNum = polyPoiList[polyNum][i];
			if (poiCondition[poiNum] < 0) { continue; }
			face_vertex->push_back ( poiNum - poiCondition[poiNum] );
		}
	}

	// ���_�̈ʒu�f�[�^
	for (uInt i = 0; i < poiCount; ++i) {
		if (poiCondition[i] < 0) { continue; }
		point_position->push_back ( source_point_position->at ( i ) );
	}	

	point_position->shrink_to_fit();
	face_vertex->shrink_to_fit();
	face_offset->shrink_to_fit();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////// step-7

	//////////////////////////////////////////////////////////////////////////////////////////////////////////// step-8
	// �v���p�e�B�\�������p��

/*	if (transferNormalProp || transferUVsProp || transferOtherProp){
		// ���_���邢�̓|���S���̍폜�ɂ��o�[�e�b�N�X�̍폜�\���Ԃ𒲂ׂĂ���
		if (mode == TKCM::Delete::Mode::Delete) {
			// �폜�\��̃o�[�e�b�N�X�����X�g�ɂ��� (false=�폜���Ȃ��Atrue=�폜�\��)
			Amino::Array<bool> verCondition ( verCount, false );

			for (uInt polyNum = 0; polyNum < polyCount; ++polyNum) {
				if (polyCondition.at ( polyNum ) == true || polyPoiList[polyNum].size () < 3) {
					// �|���S�����o�͑Ώۂł͂Ȃ������ꍇ
					for (uInt verNum = source_face_offset->at ( polyNum ); verNum < source_face_offset->at ( polyNum + 1 ); ++verNum) {
						// �|���S�����\������S�Ẵo�[�e�b�N�X���폜�\��ɂ���
						verCondition[verNum] = true;
					}
				} else {
					// �|���S����resukt�ɏo�͂���Ă���ꍇ
					for (uInt verNum = source_face_offset->at ( polyNum ); verNum < source_face_offset->at ( polyNum + 1 ); ++verNum) {
						int polySize = 0;// �|���S���̉搔�̃J�E���^�[
						for (size_t i = 0; i < polyPoiList[polyNum].size (); ++i) {
							int poiID = polyPoiList[polyNum][i];
							// ���_���폜�\��̏ꍇ�̓o�[�e�b�N�X���폜�\��ɂ��ăJ�E���g�A�b�v���X�L�b�v
							if (poiCondition[poiID] < 0) {
								verCondition[source_face_offset->at ( polyNum ) + i] = true;
								continue;
							}
							polySize++;// �|���S���̉搔���C���N�������g
						}
						// �|���S���̉搔���s���l�̏ꍇ
						if (polySize < 3) {
							// �|���S�����\������S�Ẵo�[�e�b�N�X���폜�\��ɂ���
							for (uInt verNum = source_face_offset->at ( polyNum ); verNum < source_face_offset->at ( polyNum + 1 ); ++verNum) {
								verCondition[verNum] = true;
							}
						}
					}
				}
			}
			if (transferNormalProp) { // �@��
				TKCM::Delete::Mesh::CreateNormalProp ( *resultBifMesh, *source, poiCondition, verCondition, polyCondition );
			}
			if (transferUVsProp) { // UV
				TKCM::Delete::Mesh::CreateUvProp ( *resultBifMesh, *source, verCondition );
			}
			if (transferOtherProp) {
				TKCM::Delete::Mesh::CreatePointComponentProp ( *resultBifMesh, *source, poiCondition ); // point_component�^�C�v
				TKCM::Delete::Mesh::CreateFaceComponentProp ( *resultBifMesh, *source, polyCondition ); // face_component�^�C�v
				TKCM::Delete::Mesh::CreateVertexComponentProp ( *resultBifMesh, *source, verCondition ); // face_vertex_component�^�C�v

				// ���g�p�i=�S�Ă̒l��false�̃R���|�[�l���g�^�O�j�̃v���p�e�B�\���폜����
				if (unusedComponentTag) {
					TKCM::Delete::Mesh::EraseUnusedComponentTag ( *resultBifMesh );
				}
			}
		}else{ // mode == TKCM::Delete::Mode::Dissolve
			// �I���W�i������]������o�[�e�b�N�XID�����Ɋi�[���郊�X�g
			Amino::Array<uInt> vertexOrder;
			vertexOrder.reserve ( verCount );
			// Dissolve�̏������ēx���s����
			for (uInt polyNum = 0; polyNum < polyCount; ++polyNum) {
				if (propPolyCondition[polyNum] == true) { continue; } // ���݂̃|���S�����폜�Ώۂ̏ꍇ�̓X�L�b�v
				if (polyPoiList[polyNum].size () < 3) { continue; } // �|���S�����o�͑Ώۂł͂Ȃ������ꍇ�̓X�L�b�v
				
				// ���݂̃|���S�����ɍ폜����G�b�W�����݂��Ȃ���΃I���W�i���̃o�[�e�b�N�X��S�ē]�����X�g�ɓo�^����
				if (processPoly[polyNum] == false) { 
					for (uInt verNum = source_face_offset->at ( polyNum ); verNum < source_face_offset->at ( polyNum + 1 ); ++verNum) {
						vertexOrder.push_back ( verNum );
					}
					continue; 
				} 
				
				// �}�[�W�������s�����|���S���ԍ����L�^���郊�X�g
				Amino::Array<int> processedPolyNum;

				int startPoiNum = origPolyPoiList[polyNum][0];
				int thisPoiNum = origPolyPoiList[polyNum][0];
				int thisPolyNum = polyNum;
				int thisLocalID = 0;
				int nThisPoi = polyPoiList[polyNum][0];
				int nLocalID = 0;

				while (true) {
					// polyPoiList�̒��_�ԍ��ƌ��ݍs���Ă���Dissolve�����̒��_�ԍ�����v�����ꍇ
					if (nThisPoi == thisPoiNum ) {
						if (0 <= poiCondition[thisPoiNum]) {
							// ���_���o�͂Ɋ܂܂��ꍇ�͂��̃o�[�e�b�N�X���v���p�e�B�\�]���ΏۂƂ��ċL�^����
							int verNum = source_face_offset->at ( thisPolyNum ) + thisLocalID;
							vertexOrder.push_back ( verNum );
						}
						// polyPoiList�̒��_�ԍ������ɐi�߂�
						nLocalID++;
						nThisPoi = polyPoiList[polyNum][nLocalID];
					}

					// �G�b�W�ɗאڂ���|���S�������擾����
					int adjacentPolyNum = edgeCondition[thisPolyNum][thisLocalID];
					// �G�b�W�ɗאڂ���|���S�����폜�\�肩�m�F����
					if (adjacentPolyNum == -1 || adjacentPolyNum == -2 || propPolyCondition[adjacentPolyNum] == true) {
						// �G�b�W���폜�Ώۂł͂Ȃ� || �폜����G�b�W���{�[�_�[�G�b�W�̂��ߗאڃ|���S�������݂��Ȃ� || �אڃ|���S���͍폜�\��̂����ꂩ�ꍇ
						// = ���̃G�b�W�͏o�͂Ɋ܂߂�\��������

						// ���_�ԍ�/���[�J���ԍ������_�֐i�߂�
						thisPoiNum = TKCM::Math::GetNextVal ( origPolyPoiList[thisPolyNum], thisPoiNum );
						thisLocalID = origPolyPoiList[thisPolyNum].size () - 1 == thisLocalID ? 0 : thisLocalID + 1;
					} else {
						// �אڃ|���S�����}�[�W�Ώۂ̏ꍇ = ���̃G�b�W�͏o�͂Ɋ܂܂Ȃ�

						// ���݂̃|���S���ԍ����L�^���Ă���
						processedPolyNum.push_back ( thisPolyNum );
						// �|���S���ԍ���אڃ|���S���ԍ��ɍX�V����
						thisPolyNum = adjacentPolyNum;
						// �אڃ|���S�����Ō��݂̒��_�ԍ��̃��[�J��ID���擾����
						for (thisLocalID = 0; thisLocalID < origPolyPoiList[thisPolyNum].size (); ++thisLocalID) {
							if (origPolyPoiList[thisPolyNum][thisLocalID] == thisPoiNum) { break; }
						}
					}
					
					// ���̏����Ώۂ̒��_���X�^�[�g�ɖ߂�����I��
					if (thisPolyNum == polyNum && startPoiNum == thisPoiNum) { break; }
				}
				// ����}�[�W�������s�����|���S�����ȍ~�̃��[�v�����Ŗ��g�p�ɂ���
				for (int i = 0; i < processedPolyNum.size (); ++i) {
					if (processedPolyNum[i] == polyNum) { continue; }
					propPolyCondition[processedPolyNum[i]] = true;
				}
			}
			if (transferNormalProp) { // �@��
				TKCM::Dissolve::Mesh::CreateNormalProp ( *resultBifMesh, *source, poiCondition, vertexOrder, polyCondition );
			}
			if (transferUVsProp) { // UV
				TKCM::Dissolve::Mesh::CreateUvProp ( *resultBifMesh, *source, vertexOrder );
			}
			if (transferOtherProp) {
				TKCM::Delete::Mesh::CreatePointComponentProp ( *resultBifMesh, *source, poiCondition ); // point_component�^�C�v
				TKCM::Delete::Mesh::CreateFaceComponentProp ( *resultBifMesh, *source, polyCondition ); // face_component�^�C�v
				TKCM::Dissolve::Mesh::CreateVertexComponentProp ( *resultBifMesh, *source, vertexOrder ); // face_vertex_component�^�C�v
				
				// ���g�p�i=�S�Ă̒l��false�̃R���|�[�l���g�^�O�j�̃v���p�e�B�\���폜����
				if (unusedComponentTag) {
					TKCM::Delete::Mesh::EraseUnusedComponentTag ( *resultBifMesh );
				}
			}
		}
	}
*/
	success = true;
}
