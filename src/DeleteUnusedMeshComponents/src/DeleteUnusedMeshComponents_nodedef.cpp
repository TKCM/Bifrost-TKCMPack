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
		!source_face_offset || source_face_offset->empty()
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

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 全く同じ頂点リストのポリゴンを削除（表裏の場合は削除しない）
	if (overlapPolygon){
		// ポリゴンを構成する頂点の番号を合計したリストを作成する
		Amino::Array<int> polyPoiNumSumList ( polyCount );
		#pragma omp parallel for
		for (int polyNum = 0; polyNum < polyCount; ++polyNum) {
			polyPoiNumSumList[polyNum] = 0;
			for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
				polyPoiNumSumList[polyNum] += polyPoiList[polyNum][i];
			}
		}

		// 自身より小さいポリゴン番号を対象にして、完全に重なっているポリゴンを探す
		#pragma omp parallel for
		for (int polyNum = 0; polyNum < polyCount; ++polyNum) {
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
						int nextP = TKCM::ArrayValueNextElement<int>( polyPoiList[polyNum], thisP );

						if (nextP != TKCM::ArrayValueNextElement<int>( polyPoiList[i], thisP )) {
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
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 面積が０のポリゴンを削除
	if (degenerate) {
		#pragma omp parallel for
		for (int polyNum = 0; polyNum < polyCount; ++polyNum) {
			if (polyCondition[polyNum] == true) { continue; } // 別のポリゴンへマージ予定orポリゴン削除予定に設定されている場合はスキップ

			int polySize = int ( polyPoiList[polyNum].size () ); // ポリゴンの画数を取得
			if (polySize <= 2) { polyCondition[polyNum] = true; continue; }

			const Bifrost::Math::float3& p0 = source_point_position->at ( polyPoiList[polyNum][0] ); // 最初の頂点の位置を取得
			// サブトライアングルの面積を順に算出していく
			bool zero = true;
			for (int j = 0; j < polySize - 2; ++j) {
				const Bifrost::Math::float3& p1 = source_point_position->at ( polyPoiList[polyNum][j + 1] );
				const Bifrost::Math::float3& p2 = source_point_position->at ( polyPoiList[polyNum][j + 2] );

				float area = TKCM::GetTriangleArea ( p0, p1, p2 ); // 面積
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
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
			#pragma omp parallel for
			for (size_t poiID = 0; poiID < poiCount; ++poiID) {
				if (0 == poiPolyNum[poiID].size ()) { poiCondition[poiID] = -2; }
			}
		}

		if (inlinePoint) {
			#pragma omp parallel for
			for (int poiID = 0; poiID < poiCount; ++poiID) {
				if (poiCondition[poiID] < 0) { continue; }
				if (3 <= poiPolyNum[poiID].size ()) { continue; } // ３以上のポリゴンで使用されている場合はスキップ

				switch (poiPolyNum[poiID].size ()) {
					case 2:
					{
						// 頂点で構成されている２つのポリゴンが隣接しているか確認する
						int polyNum0 = poiPolyNum[poiID][0];
						int nextPoiID0 = TKCM::Array2DValueNextElement( polyPoiList, polyNum0, poiID, poiCondition );
						int prevPoiID0 = TKCM::Array2DValuePreviousElement( polyPoiList, polyNum0, poiID, poiCondition );
						if (nextPoiID0 == prevPoiID0 || poiID == nextPoiID0 || poiID == prevPoiID0) { break; }
						int polyNum1 = poiPolyNum[poiID][1];
						int nextPoiID1 = TKCM::Array2DValueNextElement( polyPoiList, polyNum1, poiID, poiCondition );
						int prevPoiID1 = TKCM::Array2DValuePreviousElement( polyPoiList, polyNum1, poiID, poiCondition );
						if (nextPoiID1 == prevPoiID1 || poiID == nextPoiID1 || poiID == prevPoiID1) { break; }

						if (nextPoiID0 != prevPoiID1 && prevPoiID0 != nextPoiID1) { break; }// 隣接していない場合はスキップ
						// 隣接している場合はcase=1に進む
					}
					case 1:
					{
						// 直線状の頂点（=形状を構成するために使用していない頂点）か確認する
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
	// メッシュのストラクトデータを出力する
	// プロパティ―転送処理で使用するIDリストも同時に作成する（コメント"Prop"の箇所）
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

	// ポリゴンの画数のリスト/頂点番号リスト
	face_offset->push_back ( 0 );
	for (uInt polyNum = 0; polyNum < polyCount; ++polyNum) {
		// ポリゴンが削除予定に設定されている場合はスキップ
		if (polyCondition[polyNum] == true){ continue; }

		int polySize = 0;// ポリゴンの画数
		for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
			int poiID = polyPoiList[polyNum][i];
			if (poiCondition[poiID] < 0) { continue; } // 削除予定の頂点はカウントしない
			polySize++; // ポリゴンの画数をインクリメント
		}
		// ポリゴンを構成する頂点数が不足している場合はスキップ
		if (polySize <= 2){ continue; }

		// ポリゴンの画数のリストに追加（加算）
		face_offset->push_back (face_offset->back () + polySize );
		// ポイント番号を整理してポリゴンポイントリストに追加
		for (int i = 0; i < polyPoiList[polyNum].size (); ++i) {
			int poiNum = polyPoiList[polyNum][i];
			if (poiCondition[poiNum] < 0) { continue; }
			face_vertex->push_back ( poiNum - poiCondition[poiNum] );
			face_vertex_dst_to_source->push_back(source_face_offset->at(polyNum) + i); // prop
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
	face_vertex_dst_to_source->shrink_to_fit();
	face_dst_to_source->shrink_to_fit();	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	success = true;
}
