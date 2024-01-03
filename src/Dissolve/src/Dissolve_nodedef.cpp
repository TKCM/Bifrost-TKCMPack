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
	step-1 : meshからトポロジーデータを取得する
	step-2 : コンポーネントの削除処理で使用するリストを準備する
	step-3 : 削除予定リストを作成する
	step-4 : 削除予定リストを元にstep-2で準備したリストを更新する
		Delete = コンポーネント要素を削除したリストにする（メッシュに穴が開く）
		Dissolve = 元の立体を保持するようにコンポーネント要素を削除したリストにする（削除したコンポーネントを共有するポリゴンはマージ）
	step-5 : メッシュのクリーンナップ処理
		overlapPolygon = 全く同じ頂点リストで生成するポリゴンを削除するようにリストを更新する
		degenerate = 面積が０のポリゴンを削除するようにリストを更新する
		inlinePoint = 立体の構成で使用していない直線状の頂点を削除するようにリストを更新する
		unusedPoint = ポリゴン構成で使用されていない頂点を削除するようにリストを更新する
	step-6 : リストを元にメッシュトポロジーデータを新規作成する
	step-7 : メッシュデータを作成する
	step-8 : メッシュの主要なプロパティ―群を継承する（FaceEdgeやユーザーデータタイプのプロパティ―は省く）
	step-9 : メッシュデータを出力する
	*/

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
	if (transferNormalProp || transferUVsProp || transferOtherProp){
		origPolyPoiList = polyPoiList;
	}
	
	// 削除予定エッジで隣接するポリゴンがあればtrue ＝ 溶接処理が必要なポリゴン
	Amino::Array<bool> processPoly ( polyCount, false );

	///////////////////////////////////////////////////////////////////////////////////////////////////// step-3
	
	// コンポーネントの削除リストを作成する
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
	
	// コンポーネントタイプごとに削除処理を行う
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

					// restDeleteIDにdeleteIDの情報を複製した後、deleteIDの配列をバーテックス数でリセットする
					Amino::Array<bool> restDeleteID( deleteID.size() );
					restDeleteID = deleteID;
					
					deleteID.assign ( verCount, false );
					
					// リセットしたdeleteIDに対し、削除対象の頂点を構成するバーテックスをtrueにする
					#pragma omp parallel for
					for (size_t poiNum = 0; poiNum < poiCount; ++poiNum) {
						if (restDeleteID[poiNum] == false) { continue; }

						// 頂点を構成するバーテックスをtrueにする			
						for (uInt i = point_face_adjacecy_index->at ( poiNum ); i < point_face_adjacecy_index->at ( poiNum + 1 ); ++i) {
							uInt polyNum = point_face_adjacecy_face->at(i);
							uInt localID = point_face_adjacecy_side->at(i);
							uInt vertexID = source_face_offset->at ( polyNum ) + localID;
							deleteID[vertexID] = true;
						}
						valid = true;
					}

					// ポリゴンを構成する全てのバーテックスが削除対象の場合はpolyConditionをtrueにする
					#pragma omp parallel for
					for (size_t polyNum = 0; polyNum < polyCount; ++polyNum) {
						bool delPoly = true;
						for (uInt verNum = source_face_offset->at ( polyNum ); verNum < source_face_offset->at ( polyNum + 1 ); ++verNum) {
							// ポリゴンを構成する全てのバーテックスが削除予定の場合はポリゴン自体を削除予定にする
							if (deleteID[verNum] == false) { delPoly = false; break; }
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
*/
				case TKCM::ComponentType::Half_Edge:
					deleteVertexID.resize(verCount, false);
					for (int i = 0; i < point_to_dissolve->size(); ++i){
						deleteVertexID[point_to_dissolve->at(i)] = true;
					}
					break;

/*				case TKCM::ComponentType::Face: {
					if (haveAdjProp == false || deleteID.size () != polyCount) { break; }

					// restDeleteIDにdeleteIDの情報を複製した後、deleteIDの配列をバーテックス数でリセットする
					Amino::Array<bool> restDeleteID ( deleteID.size () );
					restDeleteID = deleteID;
					
					deleteID.assign ( verCount, false );

					// リセットしたdeleteIDに対し、削除対象のポリゴンを構成するバーテックスをtrueにする
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
			// メッシュのクリーンアップ処理が有効の場合にこの時点のpolyConditionを複製しておく
			if (transferNormalProp || transferUVsProp || transferOtherProp){
				propPolyCondition.resize(polyCount);
				propPolyCondition = propPolyCondition;
			}

			if(valid = true){
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
						}
						else {
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
					polyPoiList[polyNum].resize ( newPolyPoiList.size () );
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
						// 削除予定の頂点はカウントしない
						if (poiCondition[poiID] < 0) { continue; }
						// 同じ頂点番号が複数使用予定の場合はスキップ
						if (2 <= restPoiNum.size () && poiID == restPoiNum[restPoiNum.size () - 2]) {
							restPoiNum.pop_back ();
							continue;
						}
						// 頂点番号を一時保存リストに追加
						restPoiNum.push_back ( poiID );
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
					restPoiNum.shrink_to_fit ();
					if (restPoiNum.empty () || restPoiNum.size () == 0) { restPoiNum.resize ( 1 ); }
					polyPoiList[polyNum].resize ( restPoiNum.size () );
					polyPoiList[polyNum] = restPoiNum;
				}
				
				// 頂点は全て削除対象外にする
				poiCondition.assign ( poiCondition.size (), 0 );
			}
		}
	}


	// メインの削除処理はここまで
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// ここからはメッシュのクリーンナップ処理

	///////////////////////////////////////////////////////////////////////////////////////////////////// step-5
	// 全く同じ頂点リストのポリゴンを削除（表裏の場合は削除しない）
/*	if (overlapPolygon){
		// ポリゴンを構成する頂点の番号を合計したリストを作成する
		Amino::Array<int> polyPoiNumSumList ( polyCount );
		BIF::parallel_for ( tbb::blocked_range<size_t> ( 0, polyCount ), [&]( const tbb::blocked_range<size_t>& range ) {
			for (int polyNum = int(range.begin ()); polyNum != range.end (); ++polyNum) {
				if (polyCondition[polyNum] == true) {
					// 別のポリゴンへマージ予定orポリゴン削除予定に設定されている場合は不正値をセットしてからスキップ
					polyPoiNumSumList[polyNum] = -1;
					continue;
				}

				polyPoiNumSumList[polyNum] = 0;
				for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
					polyPoiNumSumList[polyNum] += polyPoiList[polyNum][i];
				}
			}
		} );

		// 自身より小さいポリゴン番号を対象にして、完全に重なっているポリゴンを探す
		BIF::parallel_for ( tbb::blocked_range<size_t> ( 0, polyCount ), [&]( const tbb::blocked_range<size_t>& range ) {
			for (int polyNum = int(range.begin ()); polyNum != range.end (); ++polyNum) {
				int polyNSum = polyPoiNumSumList[polyNum]; // ポリゴンを構成する頂点番号の合計値を取得
				if (polyNSum <= 2) { continue; } // 値が不正値の場合はスキップ

				for (int i = 0; i < polyNum; ++i) {
					// 0番から順にpolyNumが同じポリゴンを探す 自分自身まで辿り着いたら終了
					if (polyNSum == polyPoiNumSumList[i]) {
						if (polyPoiList[polyNum].size () != polyPoiList[i].size ()) { continue; } // polyNumが同じでもポリゴンの画数が異なる場合はスキップ

						// ポリゴンを構成する頂点番号とその並び順が全て同じかチェックする
						bool same = true;
						for (int j = 0; j < polyPoiList[polyNum].size (); ++j) {
							int thisP = polyPoiList[polyNum][j];
							int nextP = TKCM::Math::GetNextVal ( polyPoiList[polyNum], thisP );

							if (nextP != TKCM::Math::GetNextVal ( polyPoiList[i], thisP )) {
								same = false;
								break;
							}
						}
						if (same == false) { continue; } // 頂点番号リストが異なる場合はスキップ

						// 頂点番号リストが合致した場合
						// 後の処理でこのポリゴンは削除対象であると判定されるように、mergeTargetFaceIdにtrueをセットする
						polyCondition[polyNum] = true;
					}
				}
			}
		} );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// 面積が０のポリゴンを削除
	if (degenerate) {
		BIF::parallel_for ( tbb::blocked_range<size_t> ( 0, polyCount ), [&]( const tbb::blocked_range<size_t>& range ) {
			for (int polyNum = int(range.begin ()); polyNum != range.end (); ++polyNum) {
				if (polyCondition[polyNum] == true) { continue; } // 別のポリゴンへマージ予定orポリゴン削除予定に設定されている場合はスキップ

				int polySize = int ( polyPoiList[polyNum].size () ); // ポリゴンの画数を取得
				if (polySize <= 2) { polyCondition[polyNum] = true; continue; }

				const Bifrost::Math::float3& p0 = source_point_position->at ( polyPoiList[polyNum][0] ); // 最初の頂点の位置を取得
				// サブトライアングルの面積を順に算出していく
				bool zero = true;
				for (int j = 0; j < polySize - 2; ++j) {
					const Bifrost::Math::float3& p1 = source_point_position->at ( polyPoiList[polyNum][j + 1] );
					const Bifrost::Math::float3& p2 = source_point_position->at ( polyPoiList[polyNum][j + 2] );

					float area = TKCM::Math::GetTriangleArea ( p0, p1, p2 ); // 面積
					if (0.0001f < area) {
						// 面積が０ではないサブトライアングルが１つでもあればループ終了
						zero = false;
						break;
					}
				}
				if (zero) {
					// zero = trueのままの場合はゼロ面積ポリゴンと判定し、削除処理の対象とするためにtrueをセットする
					polyCondition[polyNum] = true;
				}
			}
			} );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	// 不要な頂点を削除
	if (inlinePoint || unusedPoint) {
		// 頂点で構成しているポリゴンIDを纏めたリスト
		Amino::Array< Amino::Array < int >> poiPolyNum ( poiCount );
		for (uInt polyNum = 0; polyNum < polyCount; ++polyNum) {
			if (polyCondition[polyNum] == true) { continue; } // ポリゴンが削除処理の対象に設定されている場合はスキップ

			// ポリゴンで使用する頂点番号を整理する
			Amino::Array < int > restPoiIDlist ( 0 );
			for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
				int poiID = polyPoiList[polyNum][i];
				if (poiCondition[poiID] < 0) { continue; } // 削除予定の頂点はリストに追加しない
				restPoiIDlist.push_back ( poiID );
			}
			// ポリゴンの画数が不正値の場合はスキップ
			if (restPoiIDlist.size() < 3) { continue; } 
			// 頂点ごとに使用しているポリゴン番号を記録していく
			for (int i = 0; i < restPoiIDlist.size(); ++i) {
				int poiID = restPoiIDlist[i]; // ポリゴンを構成する頂点番号
				poiPolyNum[poiID].push_back ( polyNum ); // 頂点番号の配列にポリゴン番号を記録する
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
					if (3 <= poiPolyNum[poiID].size ()) { continue; } // ３以上のポリゴンで使用されている場合はスキップ

					switch (poiPolyNum[poiID].size ()) {
						case 2:
						{
							// 頂点で構成されている２つのポリゴンが隣接しているか確認する
							int polyNum0 = poiPolyNum[poiID][0];
							int nextPoiID0 = TKCM::Math::GetNextVal ( polyPoiList, polyNum0, poiID, poiCondition );
							int prevPoiID0 = TKCM::Math::GetPreviousVal ( polyPoiList, polyNum0, poiID, poiCondition );
							if (nextPoiID0 == prevPoiID0 || poiID == nextPoiID0 || poiID == prevPoiID0) { break; }
							int polyNum1 = poiPolyNum[poiID][1];
							int nextPoiID1 = TKCM::Math::GetNextVal ( polyPoiList, polyNum1, poiID, poiCondition );
							int prevPoiID1 = TKCM::Math::GetPreviousVal ( polyPoiList, polyNum1, poiID, poiCondition );
							if (nextPoiID1 == prevPoiID1 || poiID == nextPoiID1 || poiID == prevPoiID1) { break; }

							if (nextPoiID0 != prevPoiID1 && prevPoiID0 != nextPoiID1) { break; }// 隣接していない場合はスキップ
							// 隣接している場合はcase=1に進む
						}
						case 1:
						{
							// 直線状の頂点（=形状を構成するために使用していない頂点）か確認する
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
	// 頂点番号の変更リストを作成
	int offset = 0;
	for (uInt i = 0; i < poiCount; ++i) {
		if (poiCondition[i] < 0) {
			offset++;
		}
		else {
			poiCondition[i] = offset;
		}
	}

	Amino::Array<uInt> newPolySize; // ポリゴンの画数のリスト (配列要素は0スタートの加算式)
	Amino::Array<uInt> newPolyPoiID; // ポリゴンを構成する頂点番号リスト
	Amino::Array<Bifrost::Math::float3> newPoiPos; // 頂点の位置データ
	face_offset->reserve ( polyCount );
	face_vertex->reserve ( verCount );
	point_position->reserve ( poiCount );

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

			// 削除予定の頂点はカウントしない
			if (poiCondition[poiID] < 0) { continue; }

			// ポリゴンの画数をインクリメント
			polySize++;
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
	}

	// 頂点の位置データ
	for (uInt i = 0; i < poiCount; ++i) {
		if (poiCondition[i] < 0) { continue; }
		point_position->push_back ( source_point_position->at ( i ) );
	}	

	point_position->shrink_to_fit();
	face_vertex->shrink_to_fit();
	face_offset->shrink_to_fit();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////// step-7

	//////////////////////////////////////////////////////////////////////////////////////////////////////////// step-8
	// プロパティ―を引き継ぐ

/*	if (transferNormalProp || transferUVsProp || transferOtherProp){
		// 頂点あるいはポリゴンの削除によるバーテックスの削除予定状態を調べておく
		if (mode == TKCM::Delete::Mode::Delete) {
			// 削除予定のバーテックスをリストにする (false=削除しない、true=削除予定)
			Amino::Array<bool> verCondition ( verCount, false );

			for (uInt polyNum = 0; polyNum < polyCount; ++polyNum) {
				if (polyCondition.at ( polyNum ) == true || polyPoiList[polyNum].size () < 3) {
					// ポリゴンが出力対象ではなかった場合
					for (uInt verNum = source_face_offset->at ( polyNum ); verNum < source_face_offset->at ( polyNum + 1 ); ++verNum) {
						// ポリゴンを構成する全てのバーテックスを削除予定にする
						verCondition[verNum] = true;
					}
				} else {
					// ポリゴンがresuktに出力されている場合
					for (uInt verNum = source_face_offset->at ( polyNum ); verNum < source_face_offset->at ( polyNum + 1 ); ++verNum) {
						int polySize = 0;// ポリゴンの画数のカウンター
						for (size_t i = 0; i < polyPoiList[polyNum].size (); ++i) {
							int poiID = polyPoiList[polyNum][i];
							// 頂点が削除予定の場合はバーテックスを削除予定にしてカウントアップをスキップ
							if (poiCondition[poiID] < 0) {
								verCondition[source_face_offset->at ( polyNum ) + i] = true;
								continue;
							}
							polySize++;// ポリゴンの画数をインクリメント
						}
						// ポリゴンの画数が不整値の場合
						if (polySize < 3) {
							// ポリゴンを構成する全てのバーテックスを削除予定にする
							for (uInt verNum = source_face_offset->at ( polyNum ); verNum < source_face_offset->at ( polyNum + 1 ); ++verNum) {
								verCondition[verNum] = true;
							}
						}
					}
				}
			}
			if (transferNormalProp) { // 法線
				TKCM::Delete::Mesh::CreateNormalProp ( *resultBifMesh, *source, poiCondition, verCondition, polyCondition );
			}
			if (transferUVsProp) { // UV
				TKCM::Delete::Mesh::CreateUvProp ( *resultBifMesh, *source, verCondition );
			}
			if (transferOtherProp) {
				TKCM::Delete::Mesh::CreatePointComponentProp ( *resultBifMesh, *source, poiCondition ); // point_componentタイプ
				TKCM::Delete::Mesh::CreateFaceComponentProp ( *resultBifMesh, *source, polyCondition ); // face_componentタイプ
				TKCM::Delete::Mesh::CreateVertexComponentProp ( *resultBifMesh, *source, verCondition ); // face_vertex_componentタイプ

				// 未使用（=全ての値がfalseのコンポーネントタグ）のプロパティ―を削除する
				if (unusedComponentTag) {
					TKCM::Delete::Mesh::EraseUnusedComponentTag ( *resultBifMesh );
				}
			}
		}else{ // mode == TKCM::Delete::Mode::Dissolve
			// オリジナルから転送するバーテックスIDを順に格納するリスト
			Amino::Array<uInt> vertexOrder;
			vertexOrder.reserve ( verCount );
			// Dissolveの処理を再度実行する
			for (uInt polyNum = 0; polyNum < polyCount; ++polyNum) {
				if (propPolyCondition[polyNum] == true) { continue; } // 現在のポリゴンが削除対象の場合はスキップ
				if (polyPoiList[polyNum].size () < 3) { continue; } // ポリゴンが出力対象ではなかった場合はスキップ
				
				// 現在のポリゴン内に削除するエッジが存在しなければオリジナルのバーテックスを全て転送リストに登録する
				if (processPoly[polyNum] == false) { 
					for (uInt verNum = source_face_offset->at ( polyNum ); verNum < source_face_offset->at ( polyNum + 1 ); ++verNum) {
						vertexOrder.push_back ( verNum );
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

				while (true) {
					// polyPoiListの頂点番号と現在行っているDissolve処理の頂点番号が一致した場合
					if (nThisPoi == thisPoiNum ) {
						if (0 <= poiCondition[thisPoiNum]) {
							// 頂点が出力に含まれる場合はこのバーテックスをプロパティ―転送対象として記録する
							int verNum = source_face_offset->at ( thisPolyNum ) + thisLocalID;
							vertexOrder.push_back ( verNum );
						}
						// polyPoiListの頂点番号を次に進める
						nLocalID++;
						nThisPoi = polyPoiList[polyNum][nLocalID];
					}

					// エッジに隣接するポリゴン情報を取得する
					int adjacentPolyNum = edgeCondition[thisPolyNum][thisLocalID];
					// エッジに隣接するポリゴンが削除予定か確認する
					if (adjacentPolyNum == -1 || adjacentPolyNum == -2 || propPolyCondition[adjacentPolyNum] == true) {
						// エッジが削除対象ではない || 削除するエッジがボーダーエッジのため隣接ポリゴンが存在しない || 隣接ポリゴンは削除予定のいずれか場合
						// = このエッジは出力に含める可能性がある

						// 頂点番号/ローカル番号を次点へ進める
						thisPoiNum = TKCM::Math::GetNextVal ( origPolyPoiList[thisPolyNum], thisPoiNum );
						thisLocalID = origPolyPoiList[thisPolyNum].size () - 1 == thisLocalID ? 0 : thisLocalID + 1;
					} else {
						// 隣接ポリゴンがマージ対象の場合 = このエッジは出力に含まない

						// 現在のポリゴン番号を記録しておく
						processedPolyNum.push_back ( thisPolyNum );
						// ポリゴン番号を隣接ポリゴン番号に更新する
						thisPolyNum = adjacentPolyNum;
						// 隣接ポリゴン内で現在の頂点番号のローカルIDを取得する
						for (thisLocalID = 0; thisLocalID < origPolyPoiList[thisPolyNum].size (); ++thisLocalID) {
							if (origPolyPoiList[thisPolyNum][thisLocalID] == thisPoiNum) { break; }
						}
					}
					
					// 次の処理対象の頂点がスタートに戻ったら終了
					if (thisPolyNum == polyNum && startPoiNum == thisPoiNum) { break; }
				}
				// 今回マージ処理を行ったポリゴンを以降のループ処理で未使用にする
				for (int i = 0; i < processedPolyNum.size (); ++i) {
					if (processedPolyNum[i] == polyNum) { continue; }
					propPolyCondition[processedPolyNum[i]] = true;
				}
			}
			if (transferNormalProp) { // 法線
				TKCM::Dissolve::Mesh::CreateNormalProp ( *resultBifMesh, *source, poiCondition, vertexOrder, polyCondition );
			}
			if (transferUVsProp) { // UV
				TKCM::Dissolve::Mesh::CreateUvProp ( *resultBifMesh, *source, vertexOrder );
			}
			if (transferOtherProp) {
				TKCM::Delete::Mesh::CreatePointComponentProp ( *resultBifMesh, *source, poiCondition ); // point_componentタイプ
				TKCM::Delete::Mesh::CreateFaceComponentProp ( *resultBifMesh, *source, polyCondition ); // face_componentタイプ
				TKCM::Dissolve::Mesh::CreateVertexComponentProp ( *resultBifMesh, *source, vertexOrder ); // face_vertex_componentタイプ
				
				// 未使用（=全ての値がfalseのコンポーネントタグ）のプロパティ―を削除する
				if (unusedComponentTag) {
					TKCM::Delete::Mesh::EraseUnusedComponentTag ( *resultBifMesh );
				}
			}
		}
	}
*/
	success = true;
}
