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
	// 出力データの準備
	point_position = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float3>>();
	face_vertex = Amino::newMutablePtr<Amino::Array<uInt>>();
	face_offset = Amino::newMutablePtr<Amino::Array<uInt>>();
	point_dst_to_source = Amino::newMutablePtr<Amino::Array<uInt>>();
	face_dst_to_source = Amino::newMutablePtr<Amino::Array<uInt>>();
	face_vertex_dst_to_source = Amino::newMutablePtr<Amino::Array<uInt>>();
	success = false;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 入力データのチェック

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
	// マージ対象の頂点番号リストretargetMap[sauceの頂点数]を作成する
	// 値＝0～ ： マージ先の頂点番号(予め検索しておいた一定距離内の頂点リスト内から最も小さい頂点番号)をセットする
	// 値＝-1　： マージ処理の対象外やマージ先となる頂点
	Amino::Array<int> retargetMap( poiCount );
	#pragma omp parallel for
	for ( int poiNum = 0; poiNum < poiCount; ++poiNum){
		// マスクが有効な場合
		if (use_tag_data) {
			// 頂点が処理対象になっていなければ-1をセットして終了
			if (invert_tag == false) {
				if (point_tag_data_to_merge->at ( poiNum ) == false) { retargetMap[poiNum] = -1; continue; }
			} else {
				if (point_tag_data_to_merge->at ( poiNum ) == true) { retargetMap[poiNum] = -1; continue; }
			}
		}
		
		// 一定距離内に自身以外の頂点が見つからなかった場合は-1をセットして終了
		if (target_point_indices->at(poiNum)->size() <= 1){
			retargetMap[poiNum] = -1;
			continue;
		}

		// 一定距離内に含まれる頂点のリストから最も若い番号を取得する
		size_t targetID = poiCount + 1; // ダミー値をセット
		for (int i = 0; i < target_point_indices->at(poiNum)->size(); ++i){
			size_t poiID = target_point_indices->at(poiNum)->at(i);
			if (use_tag_data){
				// 頂点が処理対象でなければ終了
				if (invert_tag == false){
					if (point_tag_data_to_merge->at(poiID) == false){ continue; }
				} else{
					if (point_tag_data_to_merge->at(poiID) == true){ continue; }
				}
			}
			if (poiID < targetID){ targetID = poiID; }
		}
		// 頂点リスト内の最も若い番号が自身の番号より上の場合は終了
		if (poiNum <= targetID ){
			retargetMap[poiNum] = -1;
			continue;
		}

		// マージ先の頂点番号を登録する
		retargetMap[poiNum] = int(targetID);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// マージによって発生する頂点番号のズレを事前に算出して[idOffsetMap]にセットしておく
	Amino::Array<int> idOffsetMap( poiCount );
	int offset = 0;
	for ( int i = 0; i < poiCount; ++i ){
		idOffsetMap[i] = offset;
		offset += retargetMap[i] == -1 ? 0 : 1;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// [retargetMap]と[idOffsetMap]を使って頂点マージ後の新しいface_vertexとpolygon_sizeデータを作成する
	// merge face vertex
	Amino::Array<uInt> mergedPolyPoiIDs(verCount);
	// ポリゴンを構成する頂点リストを走査する
	#pragma omp parallel for
	for ( int vertexID = 0; vertexID < verCount; ++vertexID){
		int poiID = source_face_vertex->at( vertexID );
		// 登録されている頂点番号がマージによって別の番号に置き換わると判定されているか確認
		if ( retargetMap[poiID] == -1 ){
			// 変更がない場合は番号のオフセットだけ適応する
			mergedPolyPoiIDs[vertexID] = poiID - idOffsetMap[poiID];
		} else{
			// マージ先の頂点も別のマージ先が設定されていた場合は末端まで辿り最終的なマージ先の頂点番号を見つける
			while ( true ){
				if ( retargetMap[retargetMap[poiID]] == -1 ){ break; }
				poiID = retargetMap[poiID];
			}
			mergedPolyPoiIDs[vertexID] = retargetMap[poiID] - idOffsetMap[retargetMap[poiID]];
		}
	}
		
	// 頂点のマージによってポリゴンの画数が変化する場合があるため、face_vertexとpolygon_sizeのリストを作りなおす
	Amino::Array<Amino::Array<uInt>> restPolyPoiIDList( polyCount );
	Amino::Array<Amino::Array<bool>> vertexPropCond ( polyCount ); // step-8用のリスト
	// 頂点マージ処理後の頂点リストをポリゴン毎に走査する
	#pragma omp parallel for
	for ( int polyNum = 0; polyNum < polyCount; ++polyNum){
		// ポリゴンを構成する頂点を走査し、新しい頂点リストを作成する
		for ( uInt thisVerID = source_face_offset->at ( polyNum ); thisVerID < source_face_offset->at ( polyNum + 1); ++thisVerID){
			uInt poiID = mergedPolyPoiIDs.at( thisVerID );
			// 同じ頂点番号がポリゴンの構成リストに登録されているか確認する
			bool find = false;
			for (uInt checkVerID = source_face_offset->at ( polyNum ); checkVerID < source_face_offset->at ( polyNum + 1 ); ++checkVerID) {
				if (thisVerID == checkVerID) { continue; }
				if (mergedPolyPoiIDs[checkVerID] == poiID) {
					// 同じ頂点番号があった場合、頂点番号の編集が行われる前の元の頂点番号を比較する
					if (source_face_vertex->at ( checkVerID ) < source_face_vertex->at ( thisVerID )) {
						// 元の番号が自身より若いものがあった場合は、競合が存在したと記録する
						find = true;
						break;
					}				
				}
			}
			if (find == false) {
				// 同じ頂点番号が見つからなかった場合は出力に含める頂点として記録する
				restPolyPoiIDList[polyNum].push_back ( poiID );
			}
			vertexPropCond[polyNum].push_back ( find ); // true=削除となるバーテックス
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// [restPolyPoiIDList]と[vertexPropCond]からメッシュトポロジーデータを新規作成する
	// プロパティ―転送処理で使用するIDリストも同時に作成する（コメント"Prop"の箇所）
	
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
			continue; // ポリゴンを構成する頂点数が不足している場合はスキップ
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
