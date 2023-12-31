#include "MergeMeshPoint_nodedef.h"

void MergeMeshPoint (
	const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
	const	bool& use_tag_data,
	const	Amino::Ptr<Amino::Array<bool>>& point_tag_data_to_merge,
	const	bool& invert_tag,
	const	Amino::Ptr<Amino::Array<Amino::Ptr<Amino::Array<Amino::long_t>>>>& target_point_indices,
			Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
			Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
			Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
			Amino::MutablePtr<Amino::Array<uInt>>& point_dst_to_source,
			Amino::MutablePtr<Amino::Array<uInt>>& face_dst_to_source,
			Amino::MutablePtr<Amino::Array<uInt>>& face_vertex_dst_to_source,
			bool& success
) {
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

	if (!source_point_position || source_point_position->empty () || !source_face_vertex || source_face_vertex->empty () || !source_face_offset || source_face_offset->empty ()) {
		return; 
	}
	size_t poiCount = source_point_position->size();
	size_t polyCount = source_face_offset->size () - 1;
	size_t verCount = source_face_vertex->size ();

	if (use_tag_data){
		if (!point_tag_data_to_merge || point_tag_data_to_merge->empty()){ return; }
		if (point_tag_data_to_merge->size() != poiCount){ return; }
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// �}�[�W�Ώۂ̒��_�ԍ����X�gretargetMap[sauce�̒��_��]���쐬����
	// �l��0�` �F �}�[�W��̒��_�ԍ�(�\�ߌ������Ă�������苗�����̒��_���X�g������ł����������_�ԍ�)���Z�b�g����
	// �l��-1�@�F �}�[�W�����̑ΏۊO��}�[�W��ƂȂ钸�_
	Amino::Array<int> retargetMap( poiCount );
	#pragma omp parallel for
	for ( int poiNum = 0; poiNum < poiCount; ++poiNum){
		// �}�X�N���L���ȏꍇ
		if (use_tag_data) {
			// ���_�������ΏۂɂȂ��Ă��Ȃ����-1���Z�b�g���ďI��
			if (invert_tag == false) {
				if (point_tag_data_to_merge->at ( poiNum ) == false) { retargetMap[poiNum] = -1; continue; }
			} else {
				if (point_tag_data_to_merge->at ( poiNum ) == true) { retargetMap[poiNum] = -1; continue; }
			}
		}
		
		// ��苗�����Ɏ��g�ȊO�̒��_��������Ȃ������ꍇ��-1���Z�b�g���ďI��
		if (target_point_indices->at(poiNum)->size() <= 1){
			retargetMap[poiNum] = -1;
			continue;
		}

		// ��苗�����Ɋ܂܂�钸�_�̃��X�g����ł��Ⴂ�ԍ����擾����
		size_t targetID = poiCount + 1; // �_�~�[�l���Z�b�g
		for (int i = 0; i < target_point_indices->at(poiNum)->size(); ++i){
			size_t poiID = target_point_indices->at(poiNum)->at(i);
			if (use_tag_data){
				// ���_�������ΏۂłȂ���ΏI��
				if (invert_tag == false){
					if (point_tag_data_to_merge->at(poiID) == false){ continue; }
				} else{
					if (point_tag_data_to_merge->at(poiID) == true){ continue; }
				}
			}
			if (poiID < targetID){ targetID = poiID; }
		}
		// ���_���X�g���̍ł��Ⴂ�ԍ������g�̔ԍ�����̏ꍇ�͏I��
		if (poiNum <= targetID ){
			retargetMap[poiNum] = -1;
			continue;
		}

		// �}�[�W��̒��_�ԍ���o�^����
		retargetMap[poiNum] = int(targetID);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// �}�[�W�ɂ���Ĕ������钸�_�ԍ��̃Y�������O�ɎZ�o����[idOffsetMap]�ɃZ�b�g���Ă���
	Amino::Array<int> idOffsetMap( poiCount );
	int offset = 0;
	for ( int i = 0; i < poiCount; ++i ){
		idOffsetMap[i] = offset;
		offset += retargetMap[i] == -1 ? 0 : 1;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// [retargetMap]��[idOffsetMap]���g���Ē��_�}�[�W��̐V����face_vertex��polygon_size�f�[�^���쐬����
	// merge face vertex
	Amino::Array<uInt> mergedPolyPoiIDs(verCount);
	// �|���S�����\�����钸�_���X�g�𑖍�����
	#pragma omp parallel for
	for ( int vertexID = 0; vertexID < verCount; ++vertexID){
		int poiID = source_face_vertex->at( vertexID );
		// �o�^����Ă��钸�_�ԍ����}�[�W�ɂ���ĕʂ̔ԍ��ɒu�������Ɣ��肳��Ă��邩�m�F
		if ( retargetMap[poiID] == -1 ){
			// �ύX���Ȃ��ꍇ�͔ԍ��̃I�t�Z�b�g�����K������
			mergedPolyPoiIDs[vertexID] = poiID - idOffsetMap[poiID];
		} else{
			// �}�[�W��̒��_���ʂ̃}�[�W�悪�ݒ肳��Ă����ꍇ�͖��[�܂ŒH��ŏI�I�ȃ}�[�W��̒��_�ԍ���������
			while ( true ){
				if ( retargetMap[retargetMap[poiID]] == -1 ){ break; }
				poiID = retargetMap[poiID];
			}
			mergedPolyPoiIDs[vertexID] = retargetMap[poiID] - idOffsetMap[retargetMap[poiID]];
		}
	}
		
	// ���_�̃}�[�W�ɂ���ă|���S���̉搔���ω�����ꍇ�����邽�߁Aface_vertex��polygon_size�̃��X�g�����Ȃ���
	Amino::Array<Amino::Array<uInt>> restPolyPoiIDList( polyCount );
	Amino::Array<Amino::Array<bool>> vertexPropCond ( polyCount ); // step-8�p�̃��X�g
	// ���_�}�[�W������̒��_���X�g���|���S�����ɑ�������
	#pragma omp parallel for
	for ( int polyNum = 0; polyNum < polyCount; ++polyNum){
		// �|���S�����\�����钸�_�𑖍����A�V�������_���X�g���쐬����
		for ( uInt thisVerID = source_face_offset->at ( polyNum ); thisVerID < source_face_offset->at ( polyNum + 1); ++thisVerID){
			uInt poiID = mergedPolyPoiIDs.at( thisVerID );
			// �������_�ԍ����|���S���̍\�����X�g�ɓo�^����Ă��邩�m�F����
			bool find = false;
			for (uInt checkVerID = source_face_offset->at ( polyNum ); checkVerID < source_face_offset->at ( polyNum + 1 ); ++checkVerID) {
				if (thisVerID == checkVerID) { continue; }
				if (mergedPolyPoiIDs[checkVerID] == poiID) {
					// �������_�ԍ����������ꍇ�A���_�ԍ��̕ҏW���s����O�̌��̒��_�ԍ����r����
					if (source_face_vertex->at ( checkVerID ) < source_face_vertex->at ( thisVerID )) {
						// ���̔ԍ������g���Ⴂ���̂��������ꍇ�́A���������݂����ƋL�^����
						find = true;
						break;
					}				
				}
			}
			if (find == false) {
				// �������_�ԍ���������Ȃ������ꍇ�͏o�͂Ɋ܂߂钸�_�Ƃ��ċL�^����
				restPolyPoiIDList[polyNum].push_back ( poiID );
			}
			vertexPropCond[polyNum].push_back ( find ); // true=�폜�ƂȂ�o�[�e�b�N�X
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// [restPolyPoiIDList]��[vertexPropCond]���烁�b�V���g�|���W�[�f�[�^��V�K�쐬����
	// �v���p�e�B�\�]�������Ŏg�p����ID���X�g�������ɍ쐬����i�R�����g"Prop"�̉ӏ��j
	
	// new face vertex / polygon size
	face_offset->reserve( polyCount );
	face_vertex->reserve( mergedPolyPoiIDs.size() );
	face_dst_to_source->reserve(polyCount);
	face_vertex_dst_to_source->reserve(mergedPolyPoiIDs.size());

	face_offset->push_back( 0 );
	uInt vID = 0; // Prop
	for (size_t i = 0; i < polyCount; ++i ){
		size_t polySize = restPolyPoiIDList[i].size();
		if ( polySize <= 2 ){ 
			size_t sourcePolySize = source_face_offset->at(i + 1) - source_face_offset->at(i); // Prop
			vID += sourcePolySize; // Prop
			continue; // �|���S�����\�����钸�_�����s�����Ă���ꍇ�̓X�L�b�v
		}

		face_offset->push_back( uInt(face_offset->back() + polySize) );
		face_dst_to_source->push_back(i); // Prop
		for (size_t j = 0; j < polySize; ++j ){
			face_vertex->push_back( restPolyPoiIDList[i][j] );
			face_vertex_dst_to_source->push_back(vID); // Prop
			vID++;
		}
	}
	face_offset->shrink_to_fit ();
	face_vertex->shrink_to_fit ();
	if (face_offset->size() <= 1 || face_vertex->size() < 3){ return; }

	// new point position
	point_position->reserve( poiCount );
	point_dst_to_source->reserve(poiCount);
	for ( int i = 0; i < poiCount; ++i ){
		if ( retargetMap[i] == -1 ){
			point_position->push_back( source_point_position->at( i ) );
			point_dst_to_source->push_back(i); // Prop
		}
	}
	point_position->shrink_to_fit();
	point_dst_to_source->shrink_to_fit();
	if (point_position->size() < 3){ return; }
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	success = true;
}
