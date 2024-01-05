#include "DeleteUnusedMeshComponents_nodedef.h"
#include "DeleteUnusedMeshComponents_function.h"

void TKCM::delete_unused_mesh_components_core(
	const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
	const	bool& unusedPoint,
	const	bool& inlinePoint,
	const	bool& degenerate,
	const	bool& overlapPolygon,
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
		!source_face_offset || source_face_offset->empty()
		){
		return;
	}
	size_t poiCount = source_point_position->size();
	size_t polyCount = source_face_offset->size() - 1;
	size_t verCount = source_face_vertex->size();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// �R���|�[�l���g�̍폜�����Ŏg�p���郊�X�g����������

	// ���_�̍폜�ɂ�钸�_�ԍ��̃Y�����L�^���Ă������X�g (0=�ŏI�o�͂Ɋ܂߂钸�_�Œ��_�ԍ��͂��̂܂܁A1~�ԍ��̂���鐔, -2=�폜�\��)
	Amino::Array<int> poiCondition ( poiCount, 0 );
	// �|���S���t�F�[�X�̏������e���L�^���Ă������X�g�@( false=�폜�����̑ΏۊO(�ŏI�o�͂Ɋ܂߂�j�Atrue=�폜�\��)
	Amino::Array<bool> polyCondition ( polyCount, false );
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

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// �S���������_���X�g�̃|���S�����폜�i�\���̏ꍇ�͍폜���Ȃ��j
	if (overlapPolygon){
		// �|���S�����\�����钸�_�̔ԍ������v�������X�g���쐬����
		Amino::Array<int> polyPoiNumSumList ( polyCount );
		#pragma omp parallel for
		for (int polyNum = 0; polyNum < polyCount; ++polyNum) {
			polyPoiNumSumList[polyNum] = 0;
			for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
				polyPoiNumSumList[polyNum] += polyPoiList[polyNum][i];
			}
		}

		// ���g��菬�����|���S���ԍ���Ώۂɂ��āA���S�ɏd�Ȃ��Ă���|���S����T��
		#pragma omp parallel for
		for (int polyNum = 0; polyNum < polyCount; ++polyNum) {
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
						int nextP = TKCM::ArrayValueNextElement<int>( polyPoiList[polyNum], thisP );

						if (nextP != TKCM::ArrayValueNextElement<int>( polyPoiList[i], thisP )) {
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
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// �ʐς��O�̃|���S�����폜
	if (degenerate) {
		#pragma omp parallel for
		for (int polyNum = 0; polyNum < polyCount; ++polyNum) {
			if (polyCondition[polyNum] == true) { continue; } // �ʂ̃|���S���փ}�[�W�\��or�|���S���폜�\��ɐݒ肳��Ă���ꍇ�̓X�L�b�v

			int polySize = int ( polyPoiList[polyNum].size () ); // �|���S���̉搔���擾
			if (polySize <= 2) { polyCondition[polyNum] = true; continue; }

			const Bifrost::Math::float3& p0 = source_point_position->at ( polyPoiList[polyNum][0] ); // �ŏ��̒��_�̈ʒu���擾
			// �T�u�g���C�A���O���̖ʐς����ɎZ�o���Ă���
			bool zero = true;
			for (int j = 0; j < polySize - 2; ++j) {
				const Bifrost::Math::float3& p1 = source_point_position->at ( polyPoiList[polyNum][j + 1] );
				const Bifrost::Math::float3& p2 = source_point_position->at ( polyPoiList[polyNum][j + 2] );

				float area = TKCM::GetTriangleArea ( p0, p1, p2 ); // �ʐ�
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
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
			#pragma omp parallel for
			for (size_t poiID = 0; poiID < poiCount; ++poiID) {
				if (0 == poiPolyNum[poiID].size ()) { poiCondition[poiID] = -2; }
			}
		}

		if (inlinePoint) {
			#pragma omp parallel for
			for (int poiID = 0; poiID < poiCount; ++poiID) {
				if (poiCondition[poiID] < 0) { continue; }
				if (3 <= poiPolyNum[poiID].size ()) { continue; } // �R�ȏ�̃|���S���Ŏg�p����Ă���ꍇ�̓X�L�b�v

				switch (poiPolyNum[poiID].size ()) {
					case 2:
					{
						// ���_�ō\������Ă���Q�̃|���S�����אڂ��Ă��邩�m�F����
						int polyNum0 = poiPolyNum[poiID][0];
						int nextPoiID0 = TKCM::Array2DValueNextElement( polyPoiList, polyNum0, poiID, poiCondition );
						int prevPoiID0 = TKCM::Array2DValuePreviousElement( polyPoiList, polyNum0, poiID, poiCondition );
						if (nextPoiID0 == prevPoiID0 || poiID == nextPoiID0 || poiID == prevPoiID0) { break; }
						int polyNum1 = poiPolyNum[poiID][1];
						int nextPoiID1 = TKCM::Array2DValueNextElement( polyPoiList, polyNum1, poiID, poiCondition );
						int prevPoiID1 = TKCM::Array2DValuePreviousElement( polyPoiList, polyNum1, poiID, poiCondition );
						if (nextPoiID1 == prevPoiID1 || poiID == nextPoiID1 || poiID == prevPoiID1) { break; }

						if (nextPoiID0 != prevPoiID1 && prevPoiID0 != nextPoiID1) { break; }// �אڂ��Ă��Ȃ��ꍇ�̓X�L�b�v
						// �אڂ��Ă���ꍇ��case=1�ɐi��
					}
					case 1:
					{
						// ������̒��_�i=�`����\�����邽�߂Ɏg�p���Ă��Ȃ����_�j���m�F����
						int polyNum = poiPolyNum[poiID][0];
						int nextPoiID = TKCM::ArrayValueNextElement<int>( polyPoiList[polyNum], poiID );
						int prevPoiID = TKCM::ArrayValuePreviousElement<int>( polyPoiList[polyNum], poiID );
						if (nextPoiID == prevPoiID || poiID == nextPoiID || poiID == prevPoiID) { break; }

						const Bifrost::Math::float3& p0 = source_point_position->at ( nextPoiID );
						const Bifrost::Math::float3& p1 = source_point_position->at ( poiID );
						const Bifrost::Math::float3& p2 = source_point_position->at ( prevPoiID );
						Bifrost::Math::float3 p01, p12;
						p01.x = p0.x - p1.x;
						p01.y = p0.y - p1.y;
						p01.z = p0.z - p1.z;
						p12.x = p1.x - p2.x;
						p12.y = p1.y - p2.y;
						p12.z = p1.z - p2.z;
						Bifrost::Math::float3 p0p1UnitVec3 = TKCM::Normal ( p01 );
						Bifrost::Math::float3 p1p2UnitVec3 = TKCM::Normal ( p12 );
						float dot = TKCM::Dot ( p0p1UnitVec3, p1p2UnitVec3 );

						if (TKCM::AlmostEqual ( dot, 1.0f )) {
							poiCondition[poiID] = -2;
						}
						break;
					}
				}
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// ���b�V���̃X�g���N�g�f�[�^���o�͂���
	// �v���p�e�B�\�]�������Ŏg�p����ID���X�g�������ɍ쐬����i�R�����g"Prop"�̉ӏ��j
	point_position->reserve(poiCount);
	face_vertex->reserve(verCount);
	face_offset->reserve(polyCount);
	point_dst_to_source->reserve(poiCount);
	face_vertex_dst_to_source->reserve(verCount);
	face_dst_to_source->reserve(polyCount);

	int offset = 0;
	for (uInt i = 0; i < poiCount; ++i) {
		if (poiCondition[i] < 0) {
			offset++;
		} else {
			poiCondition[i] = offset;
		}
	}

	// �|���S���̉搔�̃��X�g/���_�ԍ����X�g
	face_offset->push_back ( 0 );
	for (uInt polyNum = 0; polyNum < polyCount; ++polyNum) {
		// �|���S�����폜�\��ɐݒ肳��Ă���ꍇ�̓X�L�b�v
		if (polyCondition[polyNum] == true){ continue; }

		int polySize = 0;// �|���S���̉搔
		for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
			int poiID = polyPoiList[polyNum][i];
			if (poiCondition[poiID] < 0) { continue; } // �폜�\��̒��_�̓J�E���g���Ȃ�
			polySize++; // �|���S���̉搔���C���N�������g
		}
		// �|���S�����\�����钸�_�����s�����Ă���ꍇ�̓X�L�b�v
		if (polySize <= 2){ continue; }

		// �|���S���̉搔�̃��X�g�ɒǉ��i���Z�j
		face_offset->push_back (face_offset->back () + polySize );
		// �|�C���g�ԍ��𐮗����ă|���S���|�C���g���X�g�ɒǉ�
		for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
			int poiNum = polyPoiList[polyNum][i];
			if (poiCondition[poiNum] < 0) { continue; }
			face_vertex->push_back ( poiNum - poiCondition[poiNum] );
			face_vertex_dst_to_source->push_back(source_face_offset->at(polyNum) + i); // prop
		}
		face_dst_to_source->push_back(polyNum); // prop
	}

	// ���_�̈ʒu�f�[�^
	for (uInt i = 0; i < poiCount; ++i) {
		if (poiCondition[i] < 0) { continue; }
		point_position->push_back ( source_point_position->at ( i ) );
		point_dst_to_source->push_back(i); // prop
	}	

	point_position->shrink_to_fit();
	face_vertex->shrink_to_fit();
	face_offset->shrink_to_fit();
	point_dst_to_source->shrink_to_fit();
	face_vertex_dst_to_source->shrink_to_fit();
	face_dst_to_source->shrink_to_fit();	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	success = true;
}
