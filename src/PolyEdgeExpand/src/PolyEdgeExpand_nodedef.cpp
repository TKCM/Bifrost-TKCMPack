#include "PolyEdgeExpand_nodedef.h"

void TKCM::poly_edge_expand_core(
	const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
		const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_face_normal,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_face,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_side,
		const	Amino::Ptr<Amino::Array<uInt>>& half_edge_to_expand,
		const	bool& apply_to_adjacent_half_edge,
		const	float& offset,
		const	uInt& division,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
				Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
				Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::uint4>>&	point_dst_to_source,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::float4>>& point_weight_to_source,
				Amino::MutablePtr<Amino::Array<uInt>>& face_dst_to_source,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::uint4>>&	face_vertex_dst_to_source,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::float4>>& face_vertex_weight_to_source,
				Amino::MutablePtr<Amino::Array<bool>>& new_face_tag_data,
				bool& success
){
	// 出力データの準備
	point_position = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float3>>();
	face_vertex = Amino::newMutablePtr<Amino::Array<uInt>>();
	face_offset = Amino::newMutablePtr<Amino::Array<uInt>>();
	point_dst_to_source = Amino::newMutablePtr<Amino::Array<Bifrost::Math::uint4>>();
	point_weight_to_source = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float4>>();
	face_dst_to_source = Amino::newMutablePtr<Amino::Array<uInt>>();
	face_vertex_dst_to_source = Amino::newMutablePtr<Amino::Array<Bifrost::Math::uint4>>();
	face_vertex_weight_to_source = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float4>>();
	new_face_tag_data = Amino::newMutablePtr<Amino::Array<bool>>();
	success = false;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 入力データのチェック
	if (!source_point_position || source_point_position->empty () ||
		!source_face_offset || source_face_offset->empty() ||
		!source_face_vertex || source_face_vertex->empty () ||
		!face_ver_adjacent_edge_face || face_ver_adjacent_edge_face->empty () ||
		!face_ver_adjacent_edge_side || face_ver_adjacent_edge_side->empty() ||
		!source_face_normal || source_face_normal->empty () || 
		!half_edge_to_expand || half_edge_to_expand->empty()
	){
		return;
	}
	if (half_edge_to_expand->size() == 0){ return; }

	uInt origPoiCount = source_point_position->size ();
	uInt origPolyCount = source_face_offset->size () - 1;
	uInt origVerCount = source_face_vertex->size ();
	float _offset = std::max(0.01f, offset);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// データ準備

	// trueと入力されたハーフエッジの反対側のハーフエッジもtrueにしたリスト
	Amino::Array<bool> edgeCond(origVerCount, false);
	#pragma omp parallel for
	for (size_t i = 0; i < half_edge_to_expand->size(); ++i){
		uInt verID = half_edge_to_expand->at(i);
		if (edgeCond[verID]==true){
			continue;
		}
		edgeCond[verID] = true;
		if (apply_to_adjacent_half_edge){
			uInt adjFaceID = face_ver_adjacent_edge_face->at(verID);
			if (adjFaceID < origPolyCount){
				uInt adjVerID = source_face_offset->at(adjFaceID) + face_ver_adjacent_edge_side->at(verID);
				edgeCond[adjVerID] = true;
			}
		}
	}

	// 処理対象となるエッジの始点・終点のバーテックス番号
	Amino::Array<uInt> edgeStartVerID, edgeEndVerID;
	// 処理対象となるエッジで構成されるフェースの番号
	Amino::Array<uInt> faceIDs;
	edgeStartVerID.reserve ( origVerCount );
	edgeEndVerID.reserve ( origVerCount );
	faceIDs.reserve ( origVerCount );
	for (size_t polyNum = 0; polyNum < origPolyCount; ++polyNum) {
		for (uInt verNum = source_face_offset->at ( polyNum ); verNum < source_face_offset->at ( polyNum + 1 ); ++verNum) {
			if (edgeCond[verNum]) {
				edgeStartVerID.push_back ( verNum );
				edgeEndVerID.push_back ( TKCM::GetNextFaceVertexIndex( polyNum, verNum, source_face_offset) );
				faceIDs.push_back ( polyNum );
			}
		}
	}

	//////////////////////////////////////////// 押し出すエッジ間での接続状態を調べておく
	// エッジの始点に接続しているハーフエッジ番号（edgeStartVerIDリストのID）を記録する
	Amino::Array<int> connectEdgeID_s ( edgeStartVerID.size (), -1 );
	// エッジの終点に接続しているハーフエッジ番号（edgeStartVerIDリストのID）を記録する
	Amino::Array<int> connectEdgeID_e ( edgeStartVerID.size (), -1 );
	// 接続が無かったエッジの終点を元に新規作成した頂点番号を取得するためのオフセットを記録する　
	// 頂点リストのデータ格納順は、[ソースの頂点群、エッジ始点を元にした新規頂点群、接続の無い終点を元にした新規頂点、サブディブ、、、]となる
	// ここに記録する値は「エッジ始点を元にした新規頂点群」の末尾からの番号オフセット
	Amino::Array<int> newEndPointIDOffset ( edgeStartVerID.size (), -1 );

	#pragma omp parallel for
	for (size_t i = 0; i < edgeEndVerID.size(); ++i) {
		for (size_t j = 0; j < edgeStartVerID.size (); ++j) {
			if (edgeEndVerID[i] == edgeStartVerID[j]) {
				connectEdgeID_e[i] = j;
				connectEdgeID_s[j] = i;
				break;
			}
		}
	}
	int offsetSum = 0;
	for (size_t i = 0; i < connectEdgeID_e.size (); ++i) {
		if (0 <= connectEdgeID_e[i]) { continue; }
		newEndPointIDOffset[i] = offsetSum;
		offsetSum++;
	}
	int nonConnectedEndPointSum = offsetSum;	
	
	// 新規頂点の生成位置の計算で使用する向き情報を算出する
	Amino::Array<Bifrost::Math::float3> pushDir_s, pushDir_e;
	pushDir_s.resize ( edgeStartVerID.size () );
	pushDir_e.resize ( edgeStartVerID.size () );
	#pragma omp parallel for
	for (size_t i = 0; i < edgeStartVerID.size(); ++i) {
		// エッジの始点・終点に接続するエッジの向きと平行なベクトルを記録する
		uInt pVerID = TKCM::GetPreviousFaceVertexIndex( faceIDs[i], edgeStartVerID[i], source_face_offset);
		uInt pPoiID = source_face_vertex->at(pVerID);
		uInt nVerID = TKCM::GetNextFaceVertexIndex( faceIDs[i], edgeEndVerID[i], source_face_offset);
		uInt nPoiID = source_face_vertex->at(nVerID);
		const Bifrost::Math::float3& esPos = source_point_position->at(source_face_vertex->at(edgeStartVerID[i]));
		const Bifrost::Math::float3& eePos = source_point_position->at(source_face_vertex->at(edgeEndVerID[i]));
		const Bifrost::Math::float3& pPos = source_point_position->at(pPoiID);
		const Bifrost::Math::float3& nPos = source_point_position->at(nPoiID);
		Bifrost::Math::float3 edgePVec = TKCM::Sub(esPos, pPos);
		Bifrost::Math::float3 edgeNVec = TKCM::Sub(eePos, nPos);
		Bifrost::Math::float3 edgeVec = TKCM::Sub(eePos, esPos);

		Bifrost::Math::float3 plane_noraml_s, plane_noraml_e;
		if (TKCM::AlmostParallel ( edgePVec, edgeVec )) {
			plane_noraml_s = Bifrost::Math::float3{ 0,1,0 };
		} else {
			plane_noraml_s = TKCM::Cross ( edgePVec, edgeVec );
		}
		if (TKCM::AlmostParallel ( edgeNVec, edgeVec )) {
			plane_noraml_e = Bifrost::Math::float3{ 0,1,0 };
		} else {
			plane_noraml_e = TKCM::Cross ( edgeNVec, edgeVec );
		}
		plane_noraml_s = TKCM::Normal ( plane_noraml_s );
		plane_noraml_e = TKCM::Normal ( plane_noraml_e );
		const Bifrost::Math::float3& faceNormal = source_face_normal->at ( faceIDs[i] );
		if (TKCM::Dot ( faceNormal, plane_noraml_s ) < 0.0f) { plane_noraml_s = TKCM::Mul(plane_noraml_s ,-1.0f); }
		if (TKCM::Dot ( faceNormal, plane_noraml_e ) < 0.0f) { plane_noraml_e = TKCM::Mul(plane_noraml_e ,-1.0f); }

		pushDir_s[i] = TKCM::Cross ( plane_noraml_s, edgeVec );
		pushDir_s[i] = TKCM::Normal ( pushDir_s[i] );
		pushDir_e[i] = TKCM::Cross ( plane_noraml_e, edgeVec );
		pushDir_e[i] = TKCM::Normal ( pushDir_e[i] );
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 出力データの下準備

	// 頂点
	point_position->resize( origPoiCount + (edgeStartVerID.size () + nonConnectedEndPointSum) * (division + 1) );
	point_dst_to_source->resize(origPoiCount + (edgeStartVerID.size() + nonConnectedEndPointSum) * (division + 1)); // prop
	point_weight_to_source->resize(origPoiCount + (edgeStartVerID.size() + nonConnectedEndPointSum) * (division + 1)); // prop
	#pragma omp parallel for
	for (size_t i = 0; i < origPoiCount; ++i){
		point_position->at(i) = source_point_position->at(i);
		point_dst_to_source->at(i).x = i;
		point_weight_to_source->at(i).x = 1.0f;
	}

	// ポリゴンリスト
	// ソースメッシュのデータをコピーし、追加分のポリゴンのメモリを確保する
	face_vertex->resize( origVerCount + (edgeStartVerID.size () * 4) * (division + 1) );
	face_vertex_dst_to_source->resize(origVerCount + (edgeStartVerID.size() * 4) * (division + 1)); // prop
	face_vertex_weight_to_source->resize(origVerCount + (edgeStartVerID.size() * 4) * (division + 1)); // prop
	#pragma omp parallel for
	for (size_t i = 0; i < origVerCount; ++i){
		face_vertex->at(i) = source_face_vertex->at(i);
		face_vertex_dst_to_source->at(i).x = i; // prop
		face_vertex_weight_to_source->at(i).x = 1.0f; // prop
	}

	face_offset->resize( origPolyCount + (edgeStartVerID.size () * (division + 1)) + 1 );
	new_face_tag_data->resize(origPolyCount + (edgeStartVerID.size() * (division + 1)), false); // tag
	face_dst_to_source->resize(origPolyCount + (edgeStartVerID.size() * (division + 1))); // prop
	#pragma omp parallel for
	for (size_t i = 0; i < origPolyCount+1; ++i){
		face_offset->at(i) = source_face_offset->at(i);
	}
	#pragma omp parallel for
	for (size_t i = 0; i < origPolyCount; ++i){
		face_dst_to_source->at(i) = i;
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 押し出しに必要な頂点を追加する
	#pragma omp parallel for
	for (size_t i = 0; i < edgeStartVerID.size(); ++i) {
		// 押し出す方向
		Bifrost::Math::float3 dir;
		float reDistance = _offset;
		int prevVerID = TKCM::GetPreviousFaceVertexIndex(faceIDs[i], edgeStartVerID[i], source_face_offset);
		int prevPoiID = source_face_vertex->at(prevVerID);
		int esPoiID = source_face_vertex->at(edgeStartVerID[i]);
		int eePoiID = source_face_vertex->at(edgeEndVerID[i]);
		if (0 <= connectEdgeID_s[i]) { // エッジの始点に別の処理対象のエッジが隣接する場合
			dir = TKCM::Add(pushDir_s[i], pushDir_e[connectEdgeID_s[i]]);
		} else { // エッジの始点に別の処理対象のエッジが隣接しない場合
			const Bifrost::Math::float3& pPoiPos = source_point_position->at(prevPoiID);
			const Bifrost::Math::float3& ePoiPos = source_point_position->at(esPoiID);
			dir = TKCM::Sub(pPoiPos, ePoiPos);
			reDistance = std::min ( TKCM::Length ( dir ) - 0.01f, _offset );
		}
				
		// 隣接するエッジとのpushDirの平均を算出する（押し出しの方向）
		dir = TKCM::Normal ( dir );
		float rad = TKCM::UnitAngleTo ( dir, pushDir_s[i] );
				
		// このエッジのpushDirの方向から↑の新しい方向への角度を取得し、_offsetを元に斜辺の長さを算出する（押し出し距離）
		float hypotenuse = reDistance / std::cos ( rad );

		const Bifrost::Math::float3& poiPos = source_point_position->at(esPoiID);
		Bifrost::Math::float3& offset = TKCM::Mul(dir, hypotenuse);
		Bifrost::Math::float3 setPoiPos = TKCM::Add(poiPos, offset);
		point_position->at(origPoiCount + i) = setPoiPos;

		// prop
		point_dst_to_source->at(origPoiCount + i).x = esPoiID;
		point_dst_to_source->at(origPoiCount + i).y = eePoiID;
		point_dst_to_source->at(origPoiCount + i).z = prevPoiID;
		point_weight_to_source->at(origPoiCount + i) = 
			TKCM::GetTriBarycentric(
				setPoiPos,
				source_point_position->at(esPoiID), 
				source_point_position->at(eePoiID), 
				source_point_position->at(prevPoiID)
			);
	}
	#pragma omp parallel for
	for (size_t i = 0; i < edgeEndVerID.size(); ++i) {
		// エッジの終点に隣接する別のエッジが無い場合は、押し出す元の終点を起点にした新規頂点を登録する
		if (connectEdgeID_e[i] == -1) {
			int nextVerID = TKCM::GetNextFaceVertexIndex(faceIDs[i], edgeEndVerID[i], source_face_offset);
			int nextPoiID = source_face_vertex->at(nextVerID);
			int esPoiID = source_face_vertex->at(edgeStartVerID[i]);
			int eePoiID = source_face_vertex->at(edgeEndVerID[i]);
			// 終点に隣接するエッジ（処理対象のエッジではない）の向き
			const Bifrost::Math::float3& nPoiPos = source_point_position->at(nextPoiID);
			const Bifrost::Math::float3& ePoiPos = source_point_position->at(eePoiID);
			Bifrost::Math::float3 dir = TKCM::Sub(nPoiPos, ePoiPos);
			float reDistance = std::min ( TKCM::Length ( dir ) - 0.01f, _offset );
					
			dir = TKCM::Normal ( dir );
					
			// このエッジのpushDirの方向から↑の新しい方向への角度を取得し、_offsetを元に斜辺の長さを算出する（押し出し距離）
			float rad = TKCM::UnitAngleTo ( dir, pushDir_e[i] );
			float hypotenuse = reDistance / std::cos ( rad );

			const Bifrost::Math::float3& poiPos = source_point_position->at(eePoiID);
			Bifrost::Math::float3 offset = TKCM::Mul(dir, hypotenuse);
			int setID = origPoiCount + edgeEndVerID.size() + newEndPointIDOffset[i];
			Bifrost::Math::float3 setPoiPos = TKCM::Add(poiPos, offset);
			point_position->at(setID) = setPoiPos;
			
			// prop
			point_dst_to_source->at(setID).x = esPoiID;
			point_dst_to_source->at(setID).y = eePoiID;
			point_dst_to_source->at(setID).z = nextPoiID;
			point_weight_to_source->at(setID) =
				TKCM::GetTriBarycentric(
					setPoiPos,
					source_point_position->at(esPoiID),
					source_point_position->at(eePoiID),
					source_point_position->at(nextPoiID)
				);
		}
	}

	// 側面の分割頂点を追加する
	int sidePointCount = edgeStartVerID.size () + nonConnectedEndPointSum;
	if (0 < division ) {
		float t = 1.0f / (division + 1);
			
		#pragma omp parallel for
		for (size_t ite = 0; ite < division * edgeStartVerID.size(); ++ite) {
			uInt i = ite / edgeStartVerID.size ();
			uInt j = ite - edgeStartVerID.size () * i;
			uInt div = i + 1;

			const Bifrost::Math::float3& p0 = point_position->at(source_face_vertex->at ( edgeStartVerID[j] ));
			const Bifrost::Math::float3& p1 = point_position->at(origPoiCount + j);
			int setID = origPoiCount + (sidePointCount * div) + j;
			Bifrost::Math::float3 setPoiPos = TKCM::Lerp(p0, p1, t * (i + 1));
			point_position->at(setID) = setPoiPos;
			
			// prop
			point_dst_to_source->at(setID) = point_dst_to_source->at(origPoiCount + j);
			point_weight_to_source->at(setID) =
				TKCM::GetTriBarycentric(
					setPoiPos,
					source_point_position->at(point_dst_to_source->at(setID).x),
					source_point_position->at(point_dst_to_source->at(setID).y),
					source_point_position->at(point_dst_to_source->at(setID).z)
				);

			if (connectEdgeID_e[j] == -1) {
				const Bifrost::Math::float3& p0 = point_position->at(source_face_vertex->at ( edgeEndVerID[j] ));
				const Bifrost::Math::float3& p1 = point_position->at(origPoiCount + edgeStartVerID.size () + newEndPointIDOffset[j]);
				int setID = origPoiCount + (sidePointCount * div) + edgeStartVerID.size() + newEndPointIDOffset[j];
				Bifrost::Math::float3 setPoiPos = TKCM::Lerp(p0, p1, t * (i + 1));
				point_position->at(setID) = setPoiPos;

				// prop
				point_dst_to_source->at(setID) = point_dst_to_source->at(origPoiCount + edgeStartVerID.size() + newEndPointIDOffset[j]);
				point_weight_to_source->at(setID) =
					TKCM::GetTriBarycentric(
						setPoiPos,
						source_point_position->at(point_dst_to_source->at(setID).x),
						source_point_position->at(point_dst_to_source->at(setID).y),
						source_point_position->at(point_dst_to_source->at(setID).z)
					);
			}
		}
	}

	// ポリゴンの画数リストを作成する
	#pragma omp parallel for
	for (size_t i = 0; i < edgeStartVerID.size() * (division + 1); ++i){
		face_offset->at(origPolyCount + i + 1) = origVerCount + ((i + 1) * 4);
		new_face_tag_data->at(origPolyCount + i) = true; // tag
		face_dst_to_source->at(origPolyCount + i) = origPolyCount + 1; // prop 新規ポリゴンはスースメッシュから引き継ぐタグが無いため、範囲外の値をセットしておく
	}
	
	// ポリゴンを構成する頂点番号リストを作る
	size_t sideVerCount = edgeStartVerID.size () * 4;
	#pragma omp parallel for
	for (size_t ite = 0; ite < edgeStartVerID.size() * (division + 1); ++ite) {
		uInt i = ite / (division + 1); //for (size_t i = 0; i < edgeStartPoiID.size (); ++i) {
		uInt j = ite - (division + 1) * i; //for (uInt j = 0; j < divisions + 1; ++j) {						
		size_t i4 = i * 4;
		uInt div = j + 1;
		size_t verOffset = sideVerCount * j;

		uInt p0, p1, p2, p3;
		if (j == 0) {
			p0 = source_face_vertex->at(edgeStartVerID[i]);
			p1 = source_face_vertex->at(edgeEndVerID[i]);
		} else {
			p0 = origPoiCount + (sidePointCount * j) + i;
			if (0 <= connectEdgeID_e[i]) {
				p1 = origPoiCount + (sidePointCount * j) + connectEdgeID_e[i];
			} else {
				p1 = origPoiCount + (sidePointCount * j) + edgeStartVerID.size () + newEndPointIDOffset[i];
			}
		}
		face_vertex->at(origVerCount + verOffset + i4) = p0;
		face_vertex->at(origVerCount + verOffset + i4 + 1) = p1;
		
		if (j == division) {
			if (0 <= connectEdgeID_e[i]) {
				// エッジの終点が他の押し出しエッジと共有ではなかった場合
				// 終点用に準備しておいた頂点番号を記録する
				p2 = origPoiCount + connectEdgeID_e[i];
			} else {
				// エッジの終点が他の押し出しエッジと共有だった場合は、共有の頂点の番号を記録する
				p2 = origPoiCount + connectEdgeID_e.size () + newEndPointIDOffset[i];
			}
			p3 = origPoiCount + i;
		} else {
			if (0 <= connectEdgeID_e[i]) {
				p2 = origPoiCount + (sidePointCount * div) + connectEdgeID_e[i];
			} else {
				p2 = origPoiCount + (sidePointCount * div) + edgeStartVerID.size () + newEndPointIDOffset[i];
			}
			p3 = origPoiCount + (sidePointCount * div) + i;
		}
		face_vertex->at(origVerCount + verOffset + i4 + 2) = p2;
		face_vertex->at(origVerCount + verOffset + i4 + 3) = p3;

		// prop
		face_vertex_dst_to_source->at(origVerCount + verOffset + i4) = TKCM::ToVertexIndies(faceIDs[i], point_dst_to_source->at(p0), source_face_offset, source_face_vertex);
		face_vertex_weight_to_source->at(origVerCount + verOffset + i4) = point_weight_to_source->at(p0);

		face_vertex_dst_to_source->at(origVerCount + verOffset + i4 + 1) = TKCM::ToVertexIndies(faceIDs[i], point_dst_to_source->at(p1), source_face_offset, source_face_vertex);
		face_vertex_weight_to_source->at(origVerCount + verOffset + i4 + 1) = point_weight_to_source->at(p1);

		face_vertex_dst_to_source->at(origVerCount + verOffset + i4 + 2) = TKCM::ToVertexIndies(faceIDs[i], point_dst_to_source->at(p2), source_face_offset, source_face_vertex);
		face_vertex_weight_to_source->at(origVerCount + verOffset + i4 + 2) = point_weight_to_source->at(p2);

		face_vertex_dst_to_source->at(origVerCount + verOffset + i4 + 3) = TKCM::ToVertexIndies(faceIDs[i], point_dst_to_source->at(p3), source_face_offset, source_face_vertex);
		face_vertex_weight_to_source->at(origVerCount + verOffset + i4 + 3) = point_weight_to_source->at(p3);
	}

	// 元のポリゴンを構成する頂点番号を新規追加した頂点に変更する
	#pragma omp parallel for
	for (size_t i = 0; i < edgeStartVerID.size(); ++i) {
		int p0 = origPoiCount + i;
		face_vertex->at(edgeStartVerID[i]) = p0;
		face_vertex_dst_to_source->at(edgeStartVerID[i]) = TKCM::ToVertexIndies(faceIDs[i], point_dst_to_source->at(p0), source_face_offset, source_face_vertex);
		face_vertex_weight_to_source->at(edgeStartVerID[i]) = point_weight_to_source->at(p0);
		if (connectEdgeID_e[i] == -1) {
			int p1 = origPoiCount + edgeEndVerID.size() + newEndPointIDOffset[i];
			face_vertex->at(edgeEndVerID[i]) = p1;
			face_vertex_dst_to_source->at(edgeEndVerID[i]) = TKCM::ToVertexIndies(faceIDs[i], point_dst_to_source->at(p1), source_face_offset, source_face_vertex);
			face_vertex_weight_to_source->at(edgeEndVerID[i]) = point_weight_to_source->at(p1);
		}
	}
		
	success = true;
}

隣り合うエッジを連続して処理する際に、それぞれのハーフエッジの新規頂点の位置が微妙に誤差がある