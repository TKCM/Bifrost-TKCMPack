#include "Snap_nodedef.h"

void TKCM::snap_grid(
			Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& point_position,
	const	bool& use_tag_data,
	const	Amino::Ptr<Amino::Array<bool>>& point_tag_data_to_snap,
	const	bool& invert_tag,
	const	Bifrost::Math::uint3& grid_scale,
			Amino::MutablePtr<Amino::Array<bool>>& processed_tag_data,
			bool& success
) {
	processed_tag_data = Amino::newMutablePtr<Amino::Array<bool>>();
	success = false;
	
	if (!point_position || point_position->empty()){ return; }
	if (grid_scale.x == 0 && grid_scale.y == 0){ return; }
	if (grid_scale.y == 0 && grid_scale.z == 0){ return; }
	if (grid_scale.z == 0 && grid_scale.x == 0){ return; }

	Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>> point_position_MPtr = point_position.toMutable();
	assert(!point_position);

	uInt poiCount = point_position_MPtr->size();
	processed_tag_data->resize(poiCount, false);

	bool maskEnable = false;
	if (use_tag_data &&
		point_tag_data_to_snap && 
		point_tag_data_to_snap->empty () == false &&
		point_tag_data_to_snap->size() == poiCount) {
		maskEnable = true;
	}

	float eX, eY, eZ, sX, sY, sZ;
	eX = grid_scale.x == 0 ? 0.0f : 1.0f / float(grid_scale.x);
	eY = grid_scale.y == 0 ? 0.0f : 1.0f / float(grid_scale.y);
	eZ = grid_scale.z == 0 ? 0.0f : 1.0f / float(grid_scale.z);
	sX = float(grid_scale.x);
	sY = float(grid_scale.y);
	sZ = float(grid_scale.z);

	#pragma omp parallel for
	for (uInt poiID = 0; poiID < poiCount; ++poiID){
		if (maskEnable){
			if (invert_tag == false){
				if (point_tag_data_to_snap->at(poiID) == false){ continue; }
			} else{
				if (point_tag_data_to_snap->at(poiID) == true){ continue; }
			}
		}
		point_position_MPtr->at(poiID).x = std::round((point_position_MPtr->at(poiID).x) * eX) * sX;
		point_position_MPtr->at(poiID).y = std::round((point_position_MPtr->at(poiID).y) * eY) * sY;
		point_position_MPtr->at(poiID).z = std::round((point_position_MPtr->at(poiID).z) * eZ) * sZ;
		processed_tag_data->at(poiID) = true;
	}

	point_position = std::move(point_position_MPtr);
	assert(point_position);
	assert(point_position.unique());

	success = true;
}

void TKCM::snap_geometry(
			Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& point_position,
	const	bool& use_tag_data,
	const	Amino::Ptr<Amino::Array<bool>>& point_tag_data_to_snap,
	const	bool& invert_tag,
	const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& target_point_position,
	const	Amino::Ptr<Amino::Array<Amino::Ptr<Amino::Array<Amino::long_t>>>>& target_point_indices,
	const	bool& use_target_tag_data,
	const	Amino::Ptr<Amino::Array<bool>>& point_tag_data_of_target,
	const	bool& invert_target_tag,
	const	uInt& snap_using,
			Amino::MutablePtr<Amino::Array<bool>>& processed_point_tag_data,
			bool& success
){
	processed_point_tag_data = Amino::newMutablePtr<Amino::Array<bool>>();
	success = false;

	if (!point_position || point_position->empty() ||
		!target_point_position || target_point_position->empty() ||
		!target_point_indices || target_point_indices->empty() ||
		point_position->size() != target_point_indices->size()
		){ return; }

	Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>> point_position_MPtr = point_position.toMutable();
	assert(!point_position);

	uInt poiCount = point_position_MPtr->size();
	uInt targetPoiCount = target_point_position->size();
	processed_point_tag_data->resize(poiCount, false);

	bool maskEnable = false;
	if (use_tag_data &&
		point_tag_data_to_snap &&
		point_tag_data_to_snap->empty() == false &&
		point_tag_data_to_snap->size() == poiCount){
		maskEnable = true;
	}
	bool targetMaskEnable = false;
	if (use_target_tag_data &&
		point_tag_data_of_target &&
		point_tag_data_of_target->empty() == false &&
		point_tag_data_of_target->size() == targetPoiCount){
		targetMaskEnable = true;
	}

	TKCM::SnapTarget snapTargetType = TKCM::SnapTarget(snap_using);

	#pragma omp parallel for
	for (size_t poiID = 0; poiID < poiCount; ++poiID){
		// この頂点がマスキングされている場合はスナップ処理をスキップ
		if (maskEnable){
			if (invert_tag == false){
				if (point_tag_data_to_snap->at(poiID) == false){ continue; }
			} else{
				if (point_tag_data_to_snap->at(poiID) == true){ continue; }
			}
		}

		// 一定距離内に頂点が見つかっていない場合は終了
		uInt findedPointCount = target_point_indices->at(poiID)->size();
		if (findedPointCount == 0){ continue; }

		// スナップ先のポジション番号を検出する
		int targetID = targetPoiCount; // まず範囲外の値（実在する頂点数以上の数）をセット
		switch (snapTargetType){
			case TKCM::SnapTarget::Least_Target_Point_ID: // 頂点リストから最も若い番号を取得する
			{
				for (uInt i = 0; i < findedPointCount; ++i){
					size_t targetPoiID = target_point_indices->at(poiID)->at(i);
					// ターゲットの頂点がマスキングされている場合は、スナップ先候補リストへの追加をスキップ
					if (targetMaskEnable){
						if (invert_target_tag == false ){
							if (point_tag_data_of_target->at(targetPoiID) == false){ continue; }
						} else {
							if (point_tag_data_of_target->at(targetPoiID) == true){ continue; }
						}
					}
					if (targetPoiID < targetID){
						targetID = targetPoiID;
					}
				}
				break;
			}
			case TKCM::SnapTarget::Closest_Target_Point: // 頂点リストから最も距離の近い頂点番号を取得する
			{
				float sqLen = std::numeric_limits<float>::infinity(); // まず無限値をセット
				for (uInt i = 0; i < findedPointCount; ++i){
					size_t targetPoiID = target_point_indices->at(poiID)->at(i);
					// ターゲットの頂点がマスキングされている場合は、スナップ先候補リストへの追加をスキップ
					if (targetMaskEnable){
						if (invert_target_tag == false){
							if (point_tag_data_of_target->at(targetPoiID) == false){ continue; }
						} else{
							if (point_tag_data_of_target->at(targetPoiID) == true){ continue; }
						}
					}
					float sqDistance = TKCM::LengthSquared(point_position_MPtr->at(poiID), target_point_position->at(poiID));
					if (sqDistance < sqLen){
						targetID = targetPoiID;
						sqLen = sqDistance;
					}
				}
				break;
			}
		}
		// スナップ先候補が見つからなかった場合はスキップ
		if (targetPoiCount <= targetID ){ continue; }

		point_position_MPtr->at(poiID) = target_point_position->at(targetID);
		processed_point_tag_data->at(poiID) = true;
	}
	
	point_position = std::move(point_position_MPtr);
	assert(point_position);
	assert(point_position.unique());

	success = true;
}

void TKCM::snap_grid_point(
	const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& point_position,
	const	Bifrost::Math::uint3& grid_scale,
			Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& grid_point
){
	grid_point = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float3>>();
	if (!point_position || point_position->empty()){
		return;
	}
	if (grid_scale.x == 0 && grid_scale.y == 0){ return; }
	if (grid_scale.y == 0 && grid_scale.z == 0){ return; }
	if (grid_scale.z == 0 && grid_scale.x == 0){ return; }

	Bifrost::Math::float3 bBoxMin, bBoxMax;
	bBoxMin = point_position->at(0);
	bBoxMax = point_position->at(0);
	for (uInt i = 1; i < point_position->size(); ++i){
		bBoxMin.x = point_position->at(i).x < bBoxMin.x ? point_position->at(i).x : bBoxMin.x;
		bBoxMin.y = point_position->at(i).y < bBoxMin.y ? point_position->at(i).y : bBoxMin.y;
		bBoxMin.z = point_position->at(i).z < bBoxMin.z ? point_position->at(i).z : bBoxMin.z;

		bBoxMax.x = bBoxMax.x < point_position->at(i).x ? point_position->at(i).x : bBoxMax.x;
		bBoxMax.y = bBoxMax.y < point_position->at(i).y ? point_position->at(i).y : bBoxMax.y;
		bBoxMax.z = bBoxMax.z < point_position->at(i).z ? point_position->at(i).z : bBoxMax.z;
	}
	
	float sX = grid_scale.x == 0 ? 0.0f : 1.0f / float(grid_scale.x);
	float sY = grid_scale.y == 0 ? 0.0f : 1.0f / float(grid_scale.y);
	float sZ = grid_scale.z == 0 ? 0.0f : 1.0f / float(grid_scale.z);
	int minX = int(std::round(bBoxMin.x * sX));
	int minY = int(std::round(bBoxMin.y * sY));
	int minZ = int(std::round(bBoxMin.z * sZ));
	int maxX = int(std::round(bBoxMax.x * sX));
	int maxY = int(std::round(bBoxMax.y * sY));
	int maxZ = int(std::round(bBoxMax.z * sZ));
	int xCount = int(std::round((maxX-minX))) + 1;
	int yCount = int(std::round((maxY-minY))) + 1;
	int zCount = int(std::round((maxZ-minZ))) + 1;
	int xyCount = xCount * yCount;
	grid_point->reserve(xCount * yCount * zCount);

	for (int i = 0; i < xCount; ++i){
		Bifrost::Math::float3 v;
		v.x = float(minX + i) * grid_scale.x;
		grid_point->push_back(v);
	}
	for (int i = 0; i < yCount; ++i){
		if (i == 0){
			for (int j = 0; j < xCount; ++j){
				grid_point->at(j).y = float(minY + i) * grid_scale.y;
			}
		} else{
			for (int j = 0; j < xCount; ++j){
				Bifrost::Math::float3 v;
				v.x = grid_point->at(j).x;
				v.y = float(minY + i) * grid_scale.y;
				grid_point->push_back(v);
			}		
		}
	}
	for (int i = 0; i < zCount; ++i){
		if (i == 0){
			for (int j = 0; j < xyCount; ++j){
				grid_point->at(j).z = float(minZ + i) * grid_scale.z;
			}
		} else{
			for (int j = 0; j < xyCount; ++j){
				Bifrost::Math::float3 v;
				v.x = grid_point->at(j).x;
				v.y = grid_point->at(j).y;
				v.z = float(minZ + i) * grid_scale.z;
				grid_point->push_back(v);
			}
		}
	}
}
