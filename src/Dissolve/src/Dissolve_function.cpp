#include "Dissolve_nodedef.h"

void TKCM::dissolve_core(
	const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
	const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
	const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_face,
	const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_side,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_face,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_side,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacecy_index,
	const	TKCM::ComponentType& component_id_type,
	const	Amino::Ptr<Amino::Array<Amino::long_t>>& component_id_to_dissolve,
			Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
			Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
			Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
			Amino::MutablePtr<Amino::Array<uInt>>& point_dst_to_source,
			Amino::MutablePtr<Amino::Array<uInt>>& face_dst_to_source,
			Amino::MutablePtr<Amino::Array<uInt>>& face_vertex_dst_to_source,
			bool& success
){
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
	if (!source_point_position || source_point_position->empty() ||
		!source_face_vertex || source_face_vertex->empty() ||
		!source_face_offset || source_face_offset->empty() ||
		!face_ver_adjacent_edge_face || face_ver_adjacent_edge_face->empty() ||
		!face_ver_adjacent_edge_side || face_ver_adjacent_edge_side->empty() ||
		!point_face_adjacecy_face || point_face_adjacecy_face->empty() ||
		!point_face_adjacecy_side || point_face_adjacecy_side->empty() ||
		!point_face_adjacecy_index || point_face_adjacecy_index->empty() ||
		!component_id_to_dissolve || component_id_to_dissolve->empty()
		){
		return;
	}
	size_t poiCount = source_point_position->size();
	size_t polyCount = source_face_offset->size() - 1;
	size_t verCount = source_face_vertex->size();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// コンポーネントの削除処理で使用するリストを準備する
	
	// 頂点の削除による頂点番号のズレを記録しておくリスト (0=最終出力に含める頂点で頂点番号はそのまま、1~番号のずれる数, -2=削除予定)
	Amino::Array<int> poiCondition ( poiCount, 0 );
	// ポリゴンフェースの処理内容を記録しておくリスト　( false=削除処理の対象外(最終出力に含める）、true=削除予定)
	Amino::Array<bool> polyCondition ( polyCount, false );
	// edgeCondition[polyNum][PolyPoiSize] -1=削除処理の対象外　-2=削除予定のエッジ＋隣接ポリゴンなし 0~削除予定のエッジ＋隣接ポリゴンID
	Amino::Array<Amino::Array<int>> edgeCondition ( polyCount );
	// ポリゴンを構成する頂点番号のリスト(編集を行いやすくするためにBifrostのルールとは異なり２次元配列のリスト)
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

	// プロパティ―転送が有効の場合に使用する変数を準備
	Amino::Array<bool> propPolyCondition; // メッシュのクリーンアップ処理を行う前のpolyConditionを保持しておくための変数
	Amino::Array<Amino::Array<int>> origPolyPoiList; // メッシュのクリーンアップ処理を行う前のpolyPoiListを保持しておくための変数
//	if (transferNormalProp || transferUVsProp || transferOtherProp){
		origPolyPoiList = polyPoiList;
//	}
	
	// 削除予定エッジで隣接するポリゴンを記録する変数を準備　(true＝溶接処理が必要なポリゴン)
	Amino::Array<bool> processPoly ( polyCount, false );

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// コンポーネントの削除予定リストを作成する
	Amino::Array<bool> deleteVertexID(verCount, false);
	switch (component_id_type) {
		case TKCM::ComponentType::Point: {			
			// 削除対象の頂点を構成するバーテックスをtrueにする
			for (size_t i = 0; i < component_id_to_dissolve->size(); ++i){
				size_t poiNum = component_id_to_dissolve->at(i);
				// 頂点を構成するバーテックスを走査する			
				for (uInt j = point_face_adjacecy_index->at(poiNum); j < point_face_adjacecy_index->at(poiNum + 1); ++j){
					uInt polyNum = point_face_adjacecy_face->at(j);
					uInt localID = point_face_adjacecy_side->at(j);
					uInt vertexID = source_face_offset->at(polyNum) + localID;
					deleteVertexID[vertexID] = true;
				}
			}

			// ポリゴンを構成する全てのバーテックスが削除対象の場合はpolyConditionをtrueにする
			#pragma omp parallel for
			for (size_t polyNum = 0; polyNum < polyCount; ++polyNum) {
				bool delPoly = true;
				for (uInt verNum = source_face_offset->at ( polyNum ); verNum < source_face_offset->at ( polyNum + 1 ); ++verNum) {
					// ポリゴンを構成する全てのバーテックスが削除予定の場合はポリゴン自体を削除予定にする
					if (deleteVertexID[verNum] == false) { delPoly = false; break; }
				}
				// ↑でポリゴンが削除予定となった場合、ポリゴン内に連続するボーダーエッジが無いかチェックする
				if (delPoly) {
					bool border = false;
					for (uInt verNum1 = source_face_offset->at ( polyNum ); verNum1 < source_face_offset->at ( polyNum + 1 ) + 1; ++verNum1) {
						// ↑の+1と↓の処理は、最後のエッジの検査後にもう一度だけ最初のエッジを検査するための対応
						int verNum = verNum1;
						if (verNum1 == source_face_offset->at ( polyNum + 1 )) { verNum = source_face_offset->at ( polyNum ); }

						if (polyCount < face_ver_adjacent_edge_face->at(verNum)) {
							// ボーダーエッジだった場合
							if (border) {
								// 1つ前のエッジもボーダーだった場合はポリゴンの削除予定を解除して終了
								delPoly = false;
								break;
							} else {
								// 次回の処理のために、このエッジがボーダーであったと記しておく
								border = true;
							}
						} else {
							// 次回の処理のために、このエッジがボーダーでは無かったと記しておく
							border = false;
						}
					}
				}
				// ↑の検査の結果を受けて、ポリゴンを削除予定にする
				if (delPoly) { polyCondition[polyNum] = true; }
			}
			break;
		}
		case TKCM::ComponentType::Half_Edge:{
			for (int i = 0; i < component_id_to_dissolve->size(); ++i){
				if (component_id_to_dissolve->at(i) < verCount){
					deleteVertexID[component_id_to_dissolve->at(i)] = true;
				}
			}
			break;
		}
		case TKCM::ComponentType::Face: {
			// 削除対象のポリゴンを構成するバーテックスをtrueにする
			#pragma omp parallel for
			for (size_t i = 0; i < component_id_to_dissolve->size(); ++i){
				size_t polyNum = component_id_to_dissolve->at(i);
				int polySize = source_face_offset->at(polyNum + 1) - source_face_offset->at(polyNum);
				int polyFirstVertex = source_face_offset->at ( polyNum );
				for (int j = 0; j < polySize; ++j) {
					deleteVertexID[polyFirstVertex + j] = true;
				}
			}
			break;
		}
	}
	// メッシュのクリーンアップ処理が有効の場合にこの時点のpolyConditionを複製しておく
//	if (transferNormalProp || transferUVsProp || transferOtherProp){
		propPolyCondition = polyCondition;
//	}

	// 各ポリゴンの各エッジごとに状態を記録していく
	#pragma omp parallel for
	for (size_t polyNum = 0; polyNum < polyCount; ++polyNum) {
		int polySize = int(polyPoiList[polyNum].size ());
		edgeCondition[polyNum].resize ( polySize, -1 );
		for (int i = 0; i < polySize; ++i) {
			// ポリゴンを構成する頂点番号
			int poiID = polyPoiList[polyNum][i];
			// ポリゴン内の次点の頂点番号を取得
			int nextPoiID = TKCM::ArrayValueNextElement<int>( polyPoiList[polyNum], poiID );
			// ハーフエッジ番号(=バーテックス番号)を取得
			int halfEdgeNum = source_face_offset->at ( polyNum ) + i;

			// エッジで隣接するポリゴン番号とハーフエッジ番号を取得する
			int adjPolyNum = -1;
			int adjHalfEdgeNum;
			if (face_ver_adjacent_edge_face->at(halfEdgeNum) < polyCount) {
				adjPolyNum = face_ver_adjacent_edge_face->at(halfEdgeNum);
				adjHalfEdgeNum = source_face_offset->at ( adjPolyNum ) + face_ver_adjacent_edge_side->at(halfEdgeNum);
			}

			if (deleteVertexID[halfEdgeNum] == true) {
				// ハーフエッジが削除対象の場合
				if (adjPolyNum < 0) {
					// ハーフエッジで隣接するポリゴンが無い場合
					edgeCondition[polyNum][i] = -2;
				} else {
					// ハーフエッジで隣接するポリゴンがある場合
					edgeCondition[polyNum][i] = adjPolyNum;
					processPoly[polyNum] = true;
				}
			} else {
				// ハーフエッジが削除対象ではない場合
				if (adjPolyNum < 0) {
					// ハーフエッジで隣接するポリゴンが無い場合
					edgeCondition[polyNum][i] = -1;
				} else {
					// ハーフエッジで隣接するポリゴンがある場合
					if (deleteVertexID[adjHalfEdgeNum] == true) {
						// 反対側のハーフエッジが削除対象の場合
						edgeCondition[polyNum][i] = adjPolyNum;
						processPoly[polyNum] = true;
					} else {
						// 反対側のハーフエッジも削除対象外の場合
						edgeCondition[polyNum][i] = -1;
					}
				}
			}
		}
	}

	// ポリゴンポイントリストを順にマージする
	for (uInt polyNum = 0; polyNum < polyCount; ++polyNum) {
		if (polyCondition[polyNum] == true) { continue; } // 現在のポリゴンが削除対象の場合はスキップ
		if (processPoly[polyNum] == false) { continue; } // 現在のポリゴン内に削除するエッジが存在しなければスキップ

		// 頂点を順にマージしていくための新しいポイントリスト
		Amino::Array<int> newPolyPoiList;
		// マージ処理を行ったポリゴン番号を記録するリスト
		Amino::Array<int> processedPolyNum;

		int startPoiNum = polyPoiList[polyNum][0];
		int thisPoiNum = polyPoiList[polyNum][0];
		int thisPolyNum = polyNum;
		int thisLocalID = 0;
				
		while (true) {
			// エッジに隣接するポリゴン情報を取得する
			int adjacentPolyNum = edgeCondition[thisPolyNum][thisLocalID];
			if (adjacentPolyNum == -1 || adjacentPolyNum == -2 || polyCondition[adjacentPolyNum] == true) {
				// エッジが削除対象ではない、削除するエッジがボーダーエッジのため隣接ポリゴンが存在しない、隣接ポリゴンは削除予定のいずれか場合

				if (adjacentPolyNum == -1 || adjacentPolyNum == -2) {
					// エッジが削除対象ではない、削除するエッジの隣接ポリゴンが存在しない場合
					// 現在の頂点番号をリストに格納する
					newPolyPoiList.push_back ( thisPoiNum );
				}
				// 頂点番号/ローカル番号を次点へ進める
				thisPoiNum = TKCM::ArrayValueNextElement<int>( polyPoiList[thisPolyNum], thisPoiNum );
				thisLocalID = edgeCondition[thisPolyNum].size () - 1 == thisLocalID ? 0 : thisLocalID + 1;
			} else {
				// 隣接ポリゴンがマージ対象の場合

				// 現在のポリゴン番号を記録しておく
				processedPolyNum.push_back ( thisPolyNum );

				// ポリゴン番号を隣接ポリゴン番号に更新する
				thisPolyNum = adjacentPolyNum;
				// 隣接ポリゴン内で現在の頂点番号のローカルIDを取得する
				for (thisLocalID = 0; thisLocalID < polyPoiList[thisPolyNum].size (); ++thisLocalID) {
					if (polyPoiList[thisPolyNum][thisLocalID] == thisPoiNum) { break; }
				}
			}

			// 次の処理対象がスタートに戻ったら終了
			if (thisPolyNum == polyNum && startPoiNum == thisPoiNum) { break; }	
		}

		if (newPolyPoiList.empty () || newPolyPoiList.size () == 0) { newPolyPoiList.resize ( 1 ); }
		polyPoiList[polyNum] = newPolyPoiList;
					
		// 今回マージ処理を行ったポリゴンを削除予定にする (マージのベースにしたポリゴンはスキップする)
		for (int i = 0; i < processedPolyNum.size (); ++i) {
			if (processedPolyNum[i] == polyNum) { continue; }
			polyCondition[processedPolyNum[i]] = true;
		}
	}
				
	// ポリゴンポイントリストを整理
	for (uInt polyNum = 0; polyNum < polyCount; ++polyNum) {
		if (polyCondition[polyNum] == true) { continue; } // ポリゴンが削除処理の対象に設定されている場合はスキップ
		if (polyPoiList[polyNum].size () < 3) { continue; }

		Amino::Array<int> restPoiNum; // ポリゴンポイント番号の一時保存先
		restPoiNum.reserve ( polyPoiList[polyNum].size () );
		for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
			int poiID = polyPoiList[polyNum][i];
			if (poiCondition[poiID] < 0) { continue; } // 削除予定の頂点はカウントしない
			// 同じ頂点番号が複数使用予定の場合はスキップ
			if (2 <= restPoiNum.size () && poiID == restPoiNum[restPoiNum.size () - 2]) {
				restPoiNum.pop_back ();
				continue;
			}
			restPoiNum.push_back ( poiID ); // 頂点番号を一時保存リストに追加
		}
		// 頂点番号リストの最初と最後で同じ頂点番号が複数使用されていないか確認する
		// 頂点番号リストが[0,6,4,5,7,6]の場合は[6,4,5,7]となる
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
		
		if (restPoiNum.empty () || restPoiNum.size () == 0) { restPoiNum.resize ( 1 ); }
		polyPoiList[polyNum] = std::move(restPoiNum);
	}
				
	// 頂点は全て削除対象外にする
	poiCondition.assign ( poiCondition.size (), 0 );

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// メッシュのストラクトデータを出力する
	// プロパティ―転送処理で使用するIDリストも同時に作成する（コメント"Prop"の箇所）
	point_position->reserve(poiCount);
	face_vertex->reserve(verCount);
	face_offset->reserve(polyCount);
	point_dst_to_source->reserve(poiCount);
	face_dst_to_source->reserve(polyCount);

	int offset = 0;
	for (uInt i = 0; i < poiCount; ++i) {
		if (poiCondition[i] < 0) {
			offset++;
		} else {
			poiCondition[i] = offset;
		}
	}

	// ポリゴンの画数のリスト/頂点番号リスト
	face_offset->push_back ( 0 );
	for (uInt polyNum = 0; polyNum < polyCount; ++polyNum) {
		// ポリゴンが削除予定に設定されている場合はスキップ
		if (polyCondition[polyNum] == true) { 
			polyPoiList[polyNum].resize ( 0 );
			continue; 
		} 

		int polySize = 0;// ポリゴンの画数
		for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
			int poiID = polyPoiList[polyNum][i];
			if (poiCondition[poiID] < 0) { continue; } // 削除予定の頂点はカウントしない
			polySize++; // ポリゴンの画数をインクリメント
		}

		// ポリゴンを構成する頂点数が不足している場合はスキップ
		if (polySize <= 2) { 
			polyPoiList[polyNum].resize ( 0 ); 
			continue; 
		}

		// ポリゴンの画数のリストに追加（加算）
		face_offset->push_back (face_offset->back () + polySize );
		// ポイント番号を整理してポリゴンポイントリストに追加
		for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
			int poiNum = polyPoiList[polyNum][i];
			if (poiCondition[poiNum] < 0) { continue; }
			face_vertex->push_back ( poiNum - poiCondition[poiNum] );
		}
		face_dst_to_source->push_back(polyNum); // prop
	}

	// 頂点の位置データ
	for (uInt i = 0; i < poiCount; ++i) {
		if (poiCondition[i] < 0) { continue; }
		point_position->push_back ( source_point_position->at ( i ) );
		point_dst_to_source->push_back(i); // prop
	}	

	point_position->shrink_to_fit();
	face_vertex->shrink_to_fit();
	face_offset->shrink_to_fit();
	point_dst_to_source->shrink_to_fit();
	face_dst_to_source->shrink_to_fit();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// バーテックスプロパティ―の転送処理で使用するIDリストを作成する
	face_vertex_dst_to_source->reserve(verCount);

	// Dissolveの処理を再度実行する
	for (uInt polyNum = 0; polyNum < polyCount; ++polyNum){
		if (propPolyCondition[polyNum] == true){ continue; } // 現在のポリゴンが削除対象の場合はスキップ
		if (polyPoiList[polyNum].size() < 3){ continue; } // ポリゴンが出力対象ではなかった場合はスキップ

		// 現在のポリゴン内に削除するエッジが存在しなければオリジナルのバーテックスを全て転送リストに登録する
		if (processPoly[polyNum] == false){
			for (uInt verNum = source_face_offset->at(polyNum); verNum < source_face_offset->at(polyNum + 1); ++verNum){
				face_vertex_dst_to_source->push_back(verNum);
			}
			continue;
		}

		// マージ処理を行ったポリゴン番号を記録するリスト
		Amino::Array<int> processedPolyNum;

		int startPoiNum = origPolyPoiList[polyNum][0];
		int thisPoiNum = origPolyPoiList[polyNum][0];
		int thisPolyNum = polyNum;
		int thisLocalID = 0;
		int nThisPoi = polyPoiList[polyNum][0];
		int nLocalID = 0;

		while (true){
			// polyPoiListの頂点番号と現在行っているDissolve処理の頂点番号が一致した場合
			if (nThisPoi == thisPoiNum){
				if (0 <= poiCondition[thisPoiNum]){
					// 頂点が出力に含まれる場合はこのバーテックスをプロパティ―転送対象として記録する
					int verNum = source_face_offset->at(thisPolyNum) + thisLocalID;
					face_vertex_dst_to_source->push_back(verNum);
				}
				// polyPoiListの頂点番号を次に進める
				nLocalID++;
				nThisPoi = polyPoiList[polyNum][nLocalID];
			}

			// エッジに隣接するポリゴン情報を取得する
			int adjacentPolyNum = edgeCondition[thisPolyNum][thisLocalID];
			// エッジに隣接するポリゴンが削除予定か確認する
			if (adjacentPolyNum == -1 || adjacentPolyNum == -2 || propPolyCondition[adjacentPolyNum] == true){
				// エッジが削除対象ではない || 削除するエッジがボーダーエッジのため隣接ポリゴンが存在しない || 隣接ポリゴンは削除予定のいずれか場合
				// = このエッジは出力に含める可能性がある

				// 頂点番号/ローカル番号を次点へ進める
				thisPoiNum = TKCM::ArrayValueNextElement<int>(origPolyPoiList[thisPolyNum], thisPoiNum);
				thisLocalID = origPolyPoiList[thisPolyNum].size() - 1 == thisLocalID ? 0 : thisLocalID + 1;
			} else{
				// 隣接ポリゴンがマージ対象の場合 = このエッジは出力に含まない

				// 現在のポリゴン番号を記録しておく
				processedPolyNum.push_back(thisPolyNum);
				// ポリゴン番号を隣接ポリゴン番号に更新する
				thisPolyNum = adjacentPolyNum;
				// 隣接ポリゴン内で現在の頂点番号のローカルIDを取得する
				for (thisLocalID = 0; thisLocalID < origPolyPoiList[thisPolyNum].size(); ++thisLocalID){
					if (origPolyPoiList[thisPolyNum][thisLocalID] == thisPoiNum){ break; }
				}
			}

			// 次の処理対象の頂点がスタートに戻ったら終了
			if (thisPolyNum == polyNum && startPoiNum == thisPoiNum){ break; }
		}
	}
	face_vertex_dst_to_source->shrink_to_fit();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	success = true;
}
