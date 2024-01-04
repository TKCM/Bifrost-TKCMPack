#include "functions.h"

namespace TKCM {
namespace Delete {
namespace Mesh {
	void EraseUnusedComponentTag ( Bifrost::Object& polygonMesh ) {
		// boolタイプのプロパティ―を取得し、全ての値がfalseの場合はプロパティ―を削除する

		// 
		Amino::Array<Amino::String> propNames = BifrostExp::getGeoPropertyKeysByTarget ( polygonMesh, "face_component" );
		UInt32 compCount = TKCM::Mesh::GetPolygonCount ( polygonMesh );
		for (size_t i = 0; i < propNames.size (); ++i) {
			auto boolProp = BifrostExp::getDataGeoPropValues<Amino::bool_t> ( polygonMesh, propNames[i] );
			if (!boolProp || boolProp->empty () || compCount != boolProp->size ()) { continue; }

			bool allFalse = true;
			for (size_t j = 0; j < boolProp->size (); ++j) {
				if (boolProp->at ( j )) { allFalse = false; break; }
			}

			if (allFalse) { polygonMesh.eraseProperty ( propNames[i] ); }
		}

		propNames = BifrostExp::getGeoPropertyKeysByTarget ( polygonMesh, "face_vertex_component" );
		compCount = TKCM::Mesh::GetVertexCount ( polygonMesh );
		for (size_t i = 0; i < propNames.size (); ++i) {
			auto boolProp = BifrostExp::getDataGeoPropValues<Amino::bool_t> ( polygonMesh, propNames[i] );
			if (!boolProp || boolProp->empty () || compCount != boolProp->size ()) { continue; }

			bool allFalse = true;
			for (size_t j = 0; j < boolProp->size (); ++j) {
				if (boolProp->at ( j )) { allFalse = false; break; }
			}

			if (allFalse) { polygonMesh.eraseProperty ( propNames[i] ); }
		}

		propNames = BifrostExp::getGeoPropertyKeysByTarget ( polygonMesh, "point_component" );
		compCount = TKCM::Mesh::GetPointCount ( polygonMesh );
		for (size_t i = 0; i < propNames.size (); ++i) {
			auto boolProp = BifrostExp::getDataGeoPropValues<Amino::bool_t> ( polygonMesh, propNames[i] );
			if (!boolProp || boolProp->empty () || compCount != boolProp->size ()) { continue; }

			bool allFalse = true;
			for (size_t j = 0; j < boolProp->size (); ++j) {
				if (boolProp->at ( j )) { allFalse = false; break; }
			}

			if (allFalse) { polygonMesh.eraseProperty ( propNames[i] ); }
		}
	}

	void CreateNormalProp ( 
		Bifrost::Object& target, 
		const Bifrost::Object& source, 
		const Amino::Array<int>& poiCondition, 
		const Amino::Array<bool>& verCondition,
		const Amino::Array<bool>& faceCondition ) 
	{
		////////////////////////////////////////////////////////////////////
		// フェース法線をターゲットメッシュの新規プロパティ―として転送
		auto faceNormalPtr = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3> ( source, "face_normal" );
		if (faceNormalPtr && target.hasProperty ( "face_normal" ) == false) {
			TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::float3> ( target, source, "face_normal", faceCondition );
		}

		////////////////////////////////////////////////////////////////////
		// 頂点法線をターゲットメッシュの新規プロパティ―として転送
		auto pointNormalPtr = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3> ( source, "point_normal" );
		if (pointNormalPtr && target.hasProperty ( "point_normal" ) == false) {
			TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::float3> ( target, source, "point_normal", poiCondition );
		}

		////////////////////////////////////////////////////////////////////
		// バーテックス法線をターゲットメッシュの新規プロパティ―として転送
		auto verNormalPtr = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3> ( source, "face_vertex_normal" );
		auto verNormalIdPtr = BifrostExp::getRangeGeoPropIndices<UInt32> ( source, "face_vertex_normal_index" );
		if (!verNormalPtr || verNormalPtr->empty () || !verNormalIdPtr || verNormalIdPtr->empty ()) { return; }
		if (target.hasProperty ( "face_vertex_normal" ) || target.hasProperty ( "face_vertex_normal_index" )) { return; }
			
		UInt32 normalValCount = UInt32 ( verNormalPtr->size () );
		UInt32 verCount = UInt32 ( verNormalIdPtr->size () );
		Amino::Array<int> useCount ( normalValCount, 0 );
		for (UInt32 verNum = 0; verNum < verCount; ++verNum) {
			if (verCondition[verNum] == true) { continue; } // バーテックスが削除予定の場合はスキップ

			UInt32 normalNum = verNormalIdPtr->at ( verNum ); // バーテックスで使用している法線データの番号を取得
			useCount[normalNum] ++; // 法線データの使用回数をカウントアップする
		}

		// 法線データの番号の変化数を算出しておく
		Amino::Array<int> normalOffset ( normalValCount, 0 );
		int offset = 0;
		for (size_t j = 0; j < useCount.size (); ++j) {
			if (useCount[j] == 0) {
				offset++;
			} else {
				normalOffset[j] = offset;
			}
		}

		{ // 法線のIDリストを新規作成してターゲットメッシュのプロパティ―にセットする
			Amino::Array<UInt32> newVerNormalIDs;
			newVerNormalIDs.reserve ( verCount );
			for (UInt32 verNum = 0; verNum < verCount; ++verNum) {
				if (verCondition[verNum] == true) { continue; } // バーテックスが削除予定の場合はスキップ

				UInt32 normalNum = verNormalIdPtr->at ( verNum ); // バーテックスで使用している法線データの番号を取得
				newVerNormalIDs.push_back ( normalNum - (normalOffset[normalNum]) );
			}
			newVerNormalIDs.shrink_to_fit ();
			Amino::Ptr<Amino::Array<UInt32>> newProp = Amino::newClassPtr<Amino::Array<UInt32>> ( newVerNormalIDs );
			Amino::MutablePtr<Bifrost::Object> rangeObj = Bifrost::createObject ();
			BifrostExp::populateRangeGeoProp<Amino::Array<UInt32>> ( newProp, "face_vertex_component", *rangeObj );
			target.setProperty ( "face_vertex_normal_index", std::move ( rangeObj ) );
		}

		{ // 法線データを新規作成してターゲットメッシュのプロパティ―にセットする 
			Amino::Array<Bifrost::Math::float3> newVerNormal;
			newVerNormal.reserve ( normalValCount );
			for (UInt32 normalNum = 0; normalNum < normalValCount; ++normalNum) {
				if (0 == useCount[normalNum]) { continue; } // 法線データが未使用の場合はスキップ

				newVerNormal.push_back ( verNormalPtr->at ( normalNum ) );
			}
			newVerNormal.shrink_to_fit ();
			auto newProp = Amino::newClassPtr<Amino::Array<Bifrost::Math::float3>> ( newVerNormal );
			Amino::MutablePtr<Bifrost::Object> propObj = Bifrost::createObject ();
			BifrostExp::populateDataGeoProp<Bifrost::Math::float3> ( Bifrost::Math::float3{ 0,1,0 }, std::move ( newProp ), "face_vertex_normal_index", *propObj );
			target.setProperty ( "face_vertex_normal", std::move ( propObj ) );
		}
	}

	void CreateUvProp (
		Bifrost::Object& target,
		const Bifrost::Object& source,
		const Amino::Array<bool>& verCondition )
	{
		auto polySizePtr = BifrostExp::getDataGeoPropValues<BifGeoIndex> ( source, "face_offset" );
		auto polyPoiListPtr = BifrostExp::getDataGeoPropValues<BifGeoIndex> ( source, "face_vertex" );
		auto poiPosPtr = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3> ( source, "point_position" );
		UInt32 verCount = UInt32(polyPoiListPtr->size ());

		for (int i = 0; i < 5; ++i) {
			// "face_vertex_uv*"と"face_vertex_uv*_index"のプロパティ―名を準備する
			Amino::String uvName = "face_vertex_uv";
			if(0<i){ uvName += std::to_string ( i ).c_str (); }
			Amino::String uvIdName = uvName + "_index";

			// ソースとターゲットメッシュにUVプロパティーが存在するか確認
			auto verUvPtr = BifrostExp::getDataGeoPropValues<Bifrost::Math::float2> ( source, uvName );
			auto verUvIdPtr = BifrostExp::getRangeGeoPropIndices<UInt32> ( source, uvIdName );
			if (!verUvPtr || verUvPtr->empty () || !verUvIdPtr || verUvIdPtr->empty ()) { continue; }
			if(target.hasProperty( uvName ) || target.hasProperty ( uvIdName )) { continue; }

			UInt32 uvValCount = UInt32(verUvPtr->size ());
			// UV座標データの使用回数をカウントする (0=未使用、1~=使用回数)
			Amino::Array<int> useCount ( uvValCount, 0 );
			for (UInt32 verNum = 0; verNum < verCount; ++verNum) {
				if (verCondition[verNum] == true) { continue; } // バーテックスが削除予定の場合はスキップ
				
				UInt32 uvNum = verUvIdPtr->at ( verNum ); // バーテックスで使用しているUV座標データの番号を取得
				useCount[uvNum] ++; // UV座標データの使用回数をカウントアップする
			}

			// UV座標データの番号の変化数を算出しておく
			Amino::Array<int> uvOffset ( uvValCount, 0 );
			int offset = 0;
			for (size_t j = 0; j < useCount.size (); ++j) {
				if (useCount[j] == 0) {
					offset++;
				} else {
					uvOffset[j] = offset;
				}
			}

			{ // UVのIDリストを新規作成してターゲットメッシュのプロパティ―にセットする
				Amino::Array<UInt32> newUvIDs;
				newUvIDs.reserve ( verCount );
				for (UInt32 verNum = 0; verNum < verCount; ++verNum) {
					if (verCondition[verNum] == true) { continue; } // バーテックスが削除予定の場合はスキップ

					UInt32 uvNum = verUvIdPtr->at ( verNum ); // バーテックスで使用しているUV座標データの番号を取得
					newUvIDs.push_back ( uvNum - (uvOffset[uvNum]) ); // UV座標データの番号の変化数を適応してからリストに追加する
				}
				newUvIDs.shrink_to_fit ();
				Amino::Ptr<Amino::Array<UInt32>> newProp = Amino::newClassPtr<Amino::Array<UInt32>> ( newUvIDs );
				Amino::MutablePtr<Bifrost::Object> rangeObj = Bifrost::createObject ();
				BifrostExp::populateRangeGeoProp<Amino::Array<UInt32>> ( newProp, "face_vertex_component", *rangeObj );
				target.setProperty ( uvIdName, std::move ( rangeObj ) );
			}

			{ // UV座標データを新規作成してターゲットメッシュのプロパティ―にセットする
				Amino::Array<Bifrost::Math::float2> newUv;
				newUv.reserve ( uvValCount );
				for (UInt32 uvNum = 0; uvNum < uvValCount; ++uvNum) {
					if (0 == useCount[uvNum]) { continue; } // UV座標データが未使用の場合はスキップ
					
					newUv.push_back ( verUvPtr->at ( uvNum ) );					
				}
				newUv.shrink_to_fit ();
				auto newProp = Amino::newClassPtr<Amino::Array<Bifrost::Math::float2>> ( newUv );
				Amino::MutablePtr<Bifrost::Object> propObj = Bifrost::createObject ();
				BifrostExp::populateDataGeoProp<Bifrost::Math::float2> ( Bifrost::Math::float2{ 0,0 }, std::move ( newProp ), uvIdName, *propObj );
				target.setProperty ( uvName, std::move ( propObj ) );
			}
		}
	}

	void CreatePointComponentProp ( Bifrost::Object& target, const Bifrost::Object& source, const Amino::Array<int>& poiCondition ) {
		const Amino::Array<Amino::String> pointPropNames = BifrostExp::getGeoPropertyKeysByTarget ( source, "point_component" );
		for (int i = 0; i < pointPropNames.size (); ++i) {
			// ターゲットに既に同名のプロパティ―が存在する場合はスキップ
			if (target.hasProperty ( pointPropNames[i] )) { continue; }
			// 法線プロパティ―は転送しない
			if (pointPropNames[i] == "point_normal") { continue; }
			
			auto propPtr_char = BifrostExp::getDataGeoPropValues<Amino::char_t> ( source, pointPropNames[i] );
			if (propPtr_char) { TKCM::Delete::Mesh::CreatePointCompPropCore<Amino::char_t> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_short = BifrostExp::getDataGeoPropValues<Amino::short_t> ( source, pointPropNames[i] );
			if (propPtr_short) { TKCM::Delete::Mesh::CreatePointCompPropCore<Amino::short_t> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_int = BifrostExp::getDataGeoPropValues<Amino::int_t> ( source, pointPropNames[i] );
			if (propPtr_int) { TKCM::Delete::Mesh::CreatePointCompPropCore<Amino::int_t> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_long = BifrostExp::getDataGeoPropValues<Amino::long_t> ( source, pointPropNames[i] );
			if (propPtr_long) { TKCM::Delete::Mesh::CreatePointCompPropCore<Amino::long_t> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_uchar = BifrostExp::getDataGeoPropValues<Amino::uchar_t> ( source, pointPropNames[i] );
			if (propPtr_uchar) { TKCM::Delete::Mesh::CreatePointCompPropCore<Amino::uchar_t> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_ushort = BifrostExp::getDataGeoPropValues<Amino::ushort_t> ( source, pointPropNames[i] );
			if (propPtr_ushort) { TKCM::Delete::Mesh::CreatePointCompPropCore<Amino::ushort_t> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_uint = BifrostExp::getDataGeoPropValues<Amino::uint_t> ( source, pointPropNames[i] );
			if (propPtr_uint) { TKCM::Delete::Mesh::CreatePointCompPropCore<Amino::uint_t> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_ulong = BifrostExp::getDataGeoPropValues<Amino::ulong_t> ( source, pointPropNames[i] );
			if (propPtr_ulong) { TKCM::Delete::Mesh::CreatePointCompPropCore<Amino::ulong_t> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_bool = BifrostExp::getDataGeoPropValues<Amino::bool_t> ( source, pointPropNames[i] );
			if (propPtr_bool) { TKCM::Delete::Mesh::CreatePointCompPropCore<Amino::bool_t> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_String = BifrostExp::getDataGeoPropValues<Amino::String> ( source, pointPropNames[i] );
			if (propPtr_String) { TKCM::Delete::Mesh::CreatePointCompPropCore<Amino::String> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_float = BifrostExp::getDataGeoPropValues<Amino::float_t> ( source, pointPropNames[i] );
			if (propPtr_float) { TKCM::Delete::Mesh::CreatePointCompPropCore<Amino::float_t> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_double = BifrostExp::getDataGeoPropValues<Amino::double_t> ( source, pointPropNames[i] );
			if (propPtr_double) { TKCM::Delete::Mesh::CreatePointCompPropCore<Amino::double_t> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_int2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int2> ( source, pointPropNames[i] );
			if (propPtr_int2) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::int2> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_int2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int2x2> ( source, pointPropNames[i] );
			if (propPtr_int2x2) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::int2x2> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_int2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int2x3> ( source, pointPropNames[i] );
			if (propPtr_int2x3) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::int2x3> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_int3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int3> ( source, pointPropNames[i] );
			if (propPtr_int3) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::int3> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_float2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float2> ( source, pointPropNames[i] );
			if (propPtr_float2) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::float2> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_float4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float4> ( source, pointPropNames[i] );
			if (propPtr_float4) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::float4> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_float3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3> ( source, pointPropNames[i] );
			if (propPtr_float3) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::float3> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_float2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float2x2> ( source, pointPropNames[i] );
			if (propPtr_float2x2) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::float2x2> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_float2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float2x3> ( source, pointPropNames[i] );
			if (propPtr_float2x3) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::float2x3> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_float3x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3x3> ( source, pointPropNames[i] );
			if (propPtr_float3x3) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::float3x3> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_float4x4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float4x4> ( source, pointPropNames[i] );
			if (propPtr_float4x4) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::float4x4> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_double2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double2> ( source, pointPropNames[i] );
			if (propPtr_double2) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::double2> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_double3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double3> ( source, pointPropNames[i] );
			if (propPtr_double3) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::double3> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_double4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double4> ( source, pointPropNames[i] );
			if (propPtr_double4) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::double4> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_double2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double2x2> ( source, pointPropNames[i] );
			if (propPtr_double2x2) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::double2x2> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_double2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double2x3> ( source, pointPropNames[i] );
			if (propPtr_double2x3) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::double2x3> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_double3x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double3x3> ( source, pointPropNames[i] );
			if (propPtr_double3x3) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::double3x3> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_double4x4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double4x4> ( source, pointPropNames[i] );
			if (propPtr_double4x4) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::double4x4> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_short2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short2> ( source, pointPropNames[i] );
			if (propPtr_short2) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::short2> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_short2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short2x2> ( source, pointPropNames[i] );
			if (propPtr_short2x2) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::short2x2> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_short2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short2x3> ( source, pointPropNames[i] );
			if (propPtr_short2x3) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::short2x3> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_short3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short3> ( source, pointPropNames[i] );
			if (propPtr_short3) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::short3> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_uchar3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::uchar3> ( source, pointPropNames[i] );
			if (propPtr_uchar3) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::uchar3> ( target, source, pointPropNames[i], poiCondition ); continue; }

			auto propPtr_uchar4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::uchar4> ( source, pointPropNames[i] );
			if (propPtr_uchar4) { TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::uchar4> ( target, source, pointPropNames[i], poiCondition ); continue; }

		}
	}

	void CreateFaceComponentProp ( Bifrost::Object& target, const Bifrost::Object& source, const Amino::Array<bool>& faceCondition ) {
		const Amino::Array<Amino::String> facePropNames = BifrostExp::getGeoPropertyKeysByTarget ( source, "face_component" );
		for (int i = 0; i < facePropNames.size (); ++i) {
			// ターゲットに既に同名のプロパティ―が存在する場合はスキップ
			if (target.hasProperty ( facePropNames[i] )) { continue; }
			// 法線プロパティ―は転送しない
			if (facePropNames[i] == "face_normal") { continue; }

			auto propPtr_char = BifrostExp::getDataGeoPropValues<Amino::char_t> ( source, facePropNames[i] );
			if (propPtr_char) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Amino::char_t> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_short = BifrostExp::getDataGeoPropValues<Amino::short_t> ( source, facePropNames[i] );
			if (propPtr_short) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Amino::short_t> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_int = BifrostExp::getDataGeoPropValues<Amino::int_t> ( source, facePropNames[i] );
			if (propPtr_int) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Amino::int_t> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_long = BifrostExp::getDataGeoPropValues<Amino::long_t> ( source, facePropNames[i] );
			if (propPtr_long) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Amino::long_t> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_uchar = BifrostExp::getDataGeoPropValues<Amino::uchar_t> ( source, facePropNames[i] );
			if (propPtr_uchar) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Amino::uchar_t> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_ushort = BifrostExp::getDataGeoPropValues<Amino::ushort_t> ( source, facePropNames[i] );
			if (propPtr_ushort) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Amino::ushort_t> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_uint = BifrostExp::getDataGeoPropValues<Amino::uint_t> ( source, facePropNames[i] );
			if (propPtr_uint) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Amino::uint_t> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_ulong = BifrostExp::getDataGeoPropValues<Amino::ulong_t> ( source, facePropNames[i] );
			if (propPtr_ulong) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Amino::ulong_t> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_bool = BifrostExp::getDataGeoPropValues<Amino::bool_t> ( source, facePropNames[i] );
			if (propPtr_bool) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Amino::bool_t> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_String = BifrostExp::getDataGeoPropValues<Amino::String> ( source, facePropNames[i] );
			if (propPtr_String) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Amino::String> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_float = BifrostExp::getDataGeoPropValues<Amino::float_t> ( source, facePropNames[i] );
			if (propPtr_float) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Amino::float_t> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_double = BifrostExp::getDataGeoPropValues<Amino::double_t> ( source, facePropNames[i] );
			if (propPtr_double) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Amino::double_t> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_int2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int2> ( source, facePropNames[i] );
			if (propPtr_int2) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::int2> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_int2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int2x2> ( source, facePropNames[i] );
			if (propPtr_int2x2) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::int2x2> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_int2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int2x3> ( source, facePropNames[i] );
			if (propPtr_int2x3) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::int2x3> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_int3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int3> ( source, facePropNames[i] );
			if (propPtr_int3) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::int3> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_float2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float2> ( source, facePropNames[i] );
			if (propPtr_float2) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::float2> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_float4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float4> ( source, facePropNames[i] );
			if (propPtr_float4) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::float4> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_float3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3> ( source, facePropNames[i] );
			if (propPtr_float3) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::float3> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_float2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float2x2> ( source, facePropNames[i] );
			if (propPtr_float2x2) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::float2x2> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_float2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float2x3> ( source, facePropNames[i] );
			if (propPtr_float2x3) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::float2x3> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_float3x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3x3> ( source, facePropNames[i] );
			if (propPtr_float3x3) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::float3x3> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_float4x4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float4x4> ( source, facePropNames[i] );
			if (propPtr_float4x4) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::float4x4> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_double2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double2> ( source, facePropNames[i] );
			if (propPtr_double2) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::double2> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_double3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double3> ( source, facePropNames[i] );
			if (propPtr_double3) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::double3> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_double4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double4> ( source, facePropNames[i] );
			if (propPtr_double4) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::double4> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_double2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double2x2> ( source, facePropNames[i] );
			if (propPtr_double2x2) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::double2x2> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_double2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double2x3> ( source, facePropNames[i] );
			if (propPtr_double2x3) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::double2x3> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_double3x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double3x3> ( source, facePropNames[i] );
			if (propPtr_double3x3) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::double3x3> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_double4x4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double4x4> ( source, facePropNames[i] );
			if (propPtr_double4x4) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::double4x4> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_short2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short2> ( source, facePropNames[i] );
			if (propPtr_short2) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::short2> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_short2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short2x2> ( source, facePropNames[i] );
			if (propPtr_short2x2) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::short2x2> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_short2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short2x3> ( source, facePropNames[i] );
			if (propPtr_short2x3) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::short2x3> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_short3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short3> ( source, facePropNames[i] );
			if (propPtr_short3) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::short3> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_uchar3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::uchar3> ( source, facePropNames[i] );
			if (propPtr_uchar3) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::uchar3> ( target, source, facePropNames[i], faceCondition ); continue; }

			auto propPtr_uchar4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::uchar4> ( source, facePropNames[i] );
			if (propPtr_uchar4) { TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::uchar4> ( target, source, facePropNames[i], faceCondition ); continue; }

		}
	}

	void CreateVertexComponentProp ( Bifrost::Object& target, const Bifrost::Object& source, const Amino::Array<bool>& verCondition ) {
		const Amino::Array<Amino::String> vertexPropNames = BifrostExp::getGeoPropertyKeysByTarget ( source, "face_vertex_component" );
		for (int i = 0; i < vertexPropNames.size (); ++i) {
			// ターゲットに既に同名のプロパティ―が存在する場合はスキップ
			if (target.hasProperty ( vertexPropNames[i] )) { continue; }
			// 法線プロパティ―は転送しない
			if (vertexPropNames[i] == "face_vertex_normal_index") { continue; }
			// UVプロパティ―は転送しない
			bool uvProp = false;
			for (int j = 0; j < 9; ++j) {
				// "face_vertex_uv*"と"face_vertex_uv*_index"のプロパティ―名を準備する "face_vertex_uv*_index"
				Amino::String uvName = "face_vertex_uv";
				if (0 < j) { uvName += std::to_string ( j ).c_str (); }
				Amino::String uvIdName = uvName + "_index";

				if (vertexPropNames[i] == uvIdName) { uvProp = true; break; }
			}
			if (uvProp) { continue; }

			auto propPtr_char = BifrostExp::getDataGeoPropValues<Amino::char_t> ( source, vertexPropNames[i] );
			if (propPtr_char) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Amino::char_t> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_short = BifrostExp::getDataGeoPropValues<Amino::short_t> ( source, vertexPropNames[i] );
			if (propPtr_short) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Amino::short_t> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_int = BifrostExp::getDataGeoPropValues<Amino::int_t> ( source, vertexPropNames[i] );
			if (propPtr_int) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Amino::int_t> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_long = BifrostExp::getDataGeoPropValues<Amino::long_t> ( source, vertexPropNames[i] );
			if (propPtr_long) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Amino::long_t> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_uchar = BifrostExp::getDataGeoPropValues<Amino::uchar_t> ( source, vertexPropNames[i] );
			if (propPtr_uchar) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Amino::uchar_t> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_ushort = BifrostExp::getDataGeoPropValues<Amino::ushort_t> ( source, vertexPropNames[i] );
			if (propPtr_ushort) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Amino::ushort_t> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_uint = BifrostExp::getDataGeoPropValues<Amino::uint_t> ( source, vertexPropNames[i] );
			if (propPtr_uint) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Amino::uint_t> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_ulong = BifrostExp::getDataGeoPropValues<Amino::ulong_t> ( source, vertexPropNames[i] );
			if (propPtr_ulong) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Amino::ulong_t> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_bool = BifrostExp::getDataGeoPropValues<Amino::bool_t> ( source, vertexPropNames[i] );
			if (propPtr_bool) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Amino::bool_t> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_String = BifrostExp::getDataGeoPropValues<Amino::String> ( source, vertexPropNames[i] );
			if (propPtr_String) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Amino::String> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_float = BifrostExp::getDataGeoPropValues<Amino::float_t> ( source, vertexPropNames[i] );
			if (propPtr_float) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Amino::float_t> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_double = BifrostExp::getDataGeoPropValues<Amino::double_t> ( source, vertexPropNames[i] );
			if (propPtr_double) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Amino::double_t> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_int2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int2> ( source, vertexPropNames[i] );
			if (propPtr_int2) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::int2> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_int2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int2x2> ( source, vertexPropNames[i] );
			if (propPtr_int2x2) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::int2x2> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_int2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int2x3> ( source, vertexPropNames[i] );
			if (propPtr_int2x3) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::int2x3> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_int3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int3> ( source, vertexPropNames[i] );
			if (propPtr_int3) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::int3> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_float2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float2> ( source, vertexPropNames[i] );
			if (propPtr_float2) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::float2> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_float4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float4> ( source, vertexPropNames[i] );
			if (propPtr_float4) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::float4> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_float3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3> ( source, vertexPropNames[i] );
			if (propPtr_float3) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::float3> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_float2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float2x2> ( source, vertexPropNames[i] );
			if (propPtr_float2x2) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::float2x2> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_float2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float2x3> ( source, vertexPropNames[i] );
			if (propPtr_float2x3) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::float2x3> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_float3x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3x3> ( source, vertexPropNames[i] );
			if (propPtr_float3x3) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::float3x3> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_float4x4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float4x4> ( source, vertexPropNames[i] );
			if (propPtr_float4x4) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::float4x4> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_double2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double2> ( source, vertexPropNames[i] );
			if (propPtr_double2) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::double2> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_double3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double3> ( source, vertexPropNames[i] );
			if (propPtr_double3) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::double3> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_double4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double4> ( source, vertexPropNames[i] );
			if (propPtr_double4) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::double4> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_double2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double2x2> ( source, vertexPropNames[i] );
			if (propPtr_double2x2) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::double2x2> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_double2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double2x3> ( source, vertexPropNames[i] );
			if (propPtr_double2x3) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::double2x3> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_double3x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double3x3> ( source, vertexPropNames[i] );
			if (propPtr_double3x3) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::double3x3> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_double4x4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double4x4> ( source, vertexPropNames[i] );
			if (propPtr_double4x4) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::double4x4> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_short2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short2> ( source, vertexPropNames[i] );
			if (propPtr_short2) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::short2> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_short2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short2x2> ( source, vertexPropNames[i] );
			if (propPtr_short2x2) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::short2x2> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_short2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short2x3> ( source, vertexPropNames[i] );
			if (propPtr_short2x3) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::short2x3> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_short3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short3> ( source, vertexPropNames[i] );
			if (propPtr_short3) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::short3> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_uchar3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::uchar3> ( source, vertexPropNames[i] );
			if (propPtr_uchar3) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::uchar3> ( target, source, vertexPropNames[i], verCondition ); continue; }

			auto propPtr_uchar4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::uchar4> ( source, vertexPropNames[i] );
			if (propPtr_uchar4) { TKCM::Delete::Mesh::CreateVertexCompPropCore<Bifrost::Math::uchar4> ( target, source, vertexPropNames[i], verCondition ); continue; }
		}
	}

	template<class T>
	void CreatePointCompPropCore (
		Bifrost::Object& target,
		const Bifrost::Object& source,
		const Amino::String& propName,
		const Amino::Array<int>& poiCondition )
	{
		if (target.hasProperty ( propName )) { return; }

		auto propPtr = BifrostExp::getDataGeoPropValues<T> ( source, propName );
		size_t count = std::min ( propPtr->size (), poiCondition.size () );

		Amino::Array<T> newPropVal;
		newPropVal.reserve ( count );
		for (int i = 0; i < count; ++i) {
			if (poiCondition[i] < 0) { continue; }
			newPropVal.push_back ( propPtr->at ( i ) );
		}
		newPropVal.shrink_to_fit ();

		auto newProp = Amino::newClassPtr<Amino::Array<T>> ( newPropVal );
		Amino::MutablePtr<Bifrost::Object> propObj = Bifrost::createObject ();
		BifrostExp::populateDataGeoProp<T> ( T (), std::move ( newProp ), "point_component", *propObj );

		target.setProperty ( propName, std::move ( propObj ) );
	}

	template<class T>
	void CreateFaceCompPropCore (
		Bifrost::Object& target,
		const Bifrost::Object& source,
		const Amino::String& propName,
		const Amino::Array<bool>& faceCondition )
	{
		if (target.hasProperty ( propName )) { return; }

		auto propPtr = BifrostExp::getDataGeoPropValues<T> ( source, propName );
		size_t count = std::min ( propPtr->size (), faceCondition.size () );

		Amino::Array<T> newPropVal;
		newPropVal.reserve ( count );
		for (int i = 0; i < count; ++i) {
			if (faceCondition[i] == true) { continue; }
			newPropVal.push_back ( propPtr->at ( i ) );
		}
		newPropVal.shrink_to_fit ();

		auto newProp = Amino::newClassPtr<Amino::Array<T>> ( newPropVal );
		Amino::MutablePtr<Bifrost::Object> propObj = Bifrost::createObject ();
		BifrostExp::populateDataGeoProp<T> ( T (), std::move ( newProp ), "face_component", *propObj );

		target.setProperty ( propName, std::move ( propObj ) );
	}

	template<class T>
	void CreateVertexCompPropCore (
		Bifrost::Object& target,
		const Bifrost::Object& source,
		const Amino::String& propName,
		const Amino::Array<bool>& verCondition )
	{
		if (target.hasProperty ( propName )) { return; }

		auto propPtr = BifrostExp::getDataGeoPropValues<T> ( source, propName );
		size_t count = std::min ( propPtr->size (), verCondition.size () );

		Amino::Array<T> newPropVal;
		newPropVal.reserve ( count );
		for (int i = 0; i < count; ++i) {
			if (verCondition[i] == true) { continue; }
			newPropVal.push_back ( propPtr->at ( i ) );
		}
		newPropVal.shrink_to_fit ();

		auto newProp = Amino::newClassPtr<Amino::Array<T>> ( newPropVal );
		Amino::MutablePtr<Bifrost::Object> propObj = Bifrost::createObject ();
		BifrostExp::populateDataGeoProp<T> ( T (), std::move ( newProp ), "face_vertex_component", *propObj );

		target.setProperty ( propName, std::move ( propObj ) );
	}
}}}

namespace TKCM {
namespace Dissolve {
namespace Mesh {
	void CreateNormalProp (
		Bifrost::Object& target,
		const Bifrost::Object& source,
		const Amino::Array<int>& poiCondition,
		const Amino::Array<UInt32>& verOrder,
		const Amino::Array<bool>& faceCondition )
	{
		////////////////////////////////////////////////////////////////////
		// フェース法線をターゲットメッシュの新規プロパティ―として転送
		auto faceNormalPtr = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3> ( source, "face_normal" );
		if (faceNormalPtr && target.hasProperty ( "face_normal" ) == false) {
			TKCM::Delete::Mesh::CreateFaceCompPropCore<Bifrost::Math::float3> ( target, source, "face_normal", faceCondition );
		}

		////////////////////////////////////////////////////////////////////
		// 頂点法線を転送をターゲットメッシュの新規プロパティ―として転送
		auto pointNormalPtr = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3> ( source, "point_normal" );
		if (pointNormalPtr && target.hasProperty ( "point_normal" ) == false) {
			TKCM::Delete::Mesh::CreatePointCompPropCore<Bifrost::Math::float3> ( target, source, "point_normal", poiCondition );
		}

		////////////////////////////////////////////////////////////////////
		// バーテックス法線をターゲットメッシュの新規プロパティ―として転送
		auto verNormalPtr = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3> ( source, "face_vertex_normal" );
		auto verNormalIdPtr = BifrostExp::getRangeGeoPropIndices<UInt32> ( source, "face_vertex_normal_index" );
		if (!verNormalPtr || verNormalPtr->empty () || !verNormalIdPtr || verNormalIdPtr->empty ()) { return; }
		if (target.hasProperty ( "face_vertex_normal" ) || target.hasProperty ( "face_vertex_normal_index" )) { return; }

		UInt32 normalValCount = UInt32 ( verNormalPtr->size () );
		UInt32 verCount = UInt32 ( verOrder.size () );
		Amino::Array<UInt32> useCount ( normalValCount, 0 );
		for (UInt32 i = 0; i < verCount; ++i) {
			UInt32 normalValID = verNormalIdPtr->at ( verOrder[i] );
			useCount[normalValID] ++; // 法線データの使用回数をカウントアップする
		}

		// 法線データの番号の変化数を算出しておく
		Amino::Array<UInt32> normalOffset ( normalValCount, 0 );
		UInt32 offset = 0;
		for (size_t j = 0; j < useCount.size (); ++j) {
			if (useCount[j] == 0) {
				offset++;
			} else {
				normalOffset[j] = offset;
			}
		}

		{ // 法線のIDリストを新規作成してターゲットメッシュのプロパティ―にセットする
			Amino::Array<UInt32> newVerNormalIDs;
			newVerNormalIDs.reserve ( verCount );
			for (UInt32 verNum = 0; verNum < verCount; ++verNum) {
				UInt32 normalNum = verNormalIdPtr->at ( verOrder[verNum] ); // バーテックスで使用している法線データの番号を取得
				newVerNormalIDs.push_back ( normalNum - (normalOffset[normalNum]) );
			}
			newVerNormalIDs.shrink_to_fit ();
			Amino::Ptr<Amino::Array<UInt32>> newProp = Amino::newClassPtr<Amino::Array<UInt32>> ( newVerNormalIDs );
			Amino::MutablePtr<Bifrost::Object> rangeObj = Bifrost::createObject ();
			BifrostExp::populateRangeGeoProp<Amino::Array<UInt32>> ( newProp, "face_vertex_component", *rangeObj );
			target.setProperty ( "face_vertex_normal_index", std::move ( rangeObj ) );
		}

		{ // 法線データを新規作成してターゲットメッシュのプロパティ―にセットする
			Amino::Array<Bifrost::Math::float3> newVerNormal;
			newVerNormal.reserve ( normalValCount );
			for (UInt32 normalNum = 0; normalNum < normalValCount; ++normalNum) {
				if (0 == useCount[normalNum]) { continue; } // 法線データが未使用の場合はスキップ

				newVerNormal.push_back ( verNormalPtr->at ( normalNum ) );
			}
			newVerNormal.shrink_to_fit ();
			auto newProp = Amino::newClassPtr<Amino::Array<Bifrost::Math::float3>> ( newVerNormal );
			Amino::MutablePtr<Bifrost::Object> propObj = Bifrost::createObject ();
			BifrostExp::populateDataGeoProp<Bifrost::Math::float3> ( Bifrost::Math::float3{ 0,1,0 }, std::move ( newProp ), "face_vertex_normal_index", *propObj );
			target.setProperty ( "face_vertex_normal", std::move ( propObj ) );
		}
	}

	void CreateUvProp (
		Bifrost::Object& target,
		const Bifrost::Object& source,
		const Amino::Array<UInt32>& verOrder )
	{
		auto polySizePtr = BifrostExp::getDataGeoPropValues<BifGeoIndex> ( source, "face_offset" );
		auto polyPoiListPtr = BifrostExp::getDataGeoPropValues<BifGeoIndex> ( source, "face_vertex" );
		auto poiPosPtr = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3> ( source, "point_position" );
		UInt32 verCount = UInt32 ( verOrder.size () );

		for (int i = 0; i < 5; ++i) {
			// "face_vertex_uv*"と"face_vertex_uv*_index"のプロパティ―名を準備する
			Amino::String uvName = "face_vertex_uv";
			if (0 < i) { uvName += std::to_string ( i ).c_str (); }
			Amino::String uvIdName = uvName + "_index";

			// ソースとターゲットメッシュにUVプロパティーが存在するか確認
			auto verUvPtr = BifrostExp::getDataGeoPropValues<Bifrost::Math::float2> ( source, uvName );
			auto verUvIdPtr = BifrostExp::getRangeGeoPropIndices<UInt32> ( source, uvIdName );
			if (!verUvPtr || verUvPtr->empty () || !verUvIdPtr || verUvIdPtr->empty ()) { continue; }
			if (target.hasProperty ( uvName ) || target.hasProperty ( uvIdName )) { continue; }

			UInt32 uvValCount = UInt32 ( verUvPtr->size () );
			// UV座標データの使用回数をカウントする (0=未使用、1~=使用回数)
			Amino::Array<int> useCount ( uvValCount, 0 );
			for (UInt32 verNum = 0; verNum < verCount; ++verNum) {
				UInt32 uvNum = verUvIdPtr->at ( verOrder[verNum] ); // バーテックスで使用しているUV座標データの番号を取得
				useCount[uvNum] ++; // UV座標データの使用回数をカウントアップする
			}

			// UV座標データの番号の変化数を算出しておく
			Amino::Array<int> uvOffset ( uvValCount, 0 );
			int offset = 0;
			for (size_t j = 0; j < useCount.size (); ++j) {
				if (useCount[j] == 0) {
					offset++;
				} else {
					uvOffset[j] = offset;
				}
			}

			{ // UVのIDリストを新規作成してターゲットメッシュのプロパティ―にセットする
				Amino::Array<UInt32> newUvIDs;
				newUvIDs.reserve ( verCount );
				for (UInt32 verNum = 0; verNum < verCount; ++verNum) {
					UInt32 uvNum = verUvIdPtr->at ( verOrder[verNum] ); // バーテックスで使用しているUV座標データの番号を取得
					newUvIDs.push_back ( uvNum - (uvOffset[uvNum]) ); // UV座標データの番号の変化数を適応してからリストに追加する
				}
				newUvIDs.shrink_to_fit ();
				Amino::Ptr<Amino::Array<UInt32>> newProp = Amino::newClassPtr<Amino::Array<UInt32>> ( newUvIDs );
				Amino::MutablePtr<Bifrost::Object> rangeObj = Bifrost::createObject ();
				BifrostExp::populateRangeGeoProp<Amino::Array<UInt32>> ( newProp, "face_vertex_component", *rangeObj );
				target.setProperty ( uvIdName, std::move ( rangeObj ) );
			}

			{ // UV座標データを新規作成してターゲットメッシュのプロパティ―にセットする
				Amino::Array<Bifrost::Math::float2> newUv;
				newUv.reserve ( uvValCount );
				for (UInt32 uvNum = 0; uvNum < uvValCount; ++uvNum) {
					if (0 == useCount[uvNum]) { continue; } // UV座標データが未使用の場合はスキップ

					newUv.push_back ( verUvPtr->at ( uvNum ) );
				}
				newUv.shrink_to_fit ();
				auto newProp = Amino::newClassPtr<Amino::Array<Bifrost::Math::float2>> ( newUv );
				Amino::MutablePtr<Bifrost::Object> propObj = Bifrost::createObject ();
				BifrostExp::populateDataGeoProp<Bifrost::Math::float2> ( Bifrost::Math::float2{ 0,0 }, std::move ( newProp ), uvIdName, *propObj );
				target.setProperty ( uvName, std::move ( propObj ) );
			}
		}
	}

	void CreateVertexComponentProp ( 
		Bifrost::Object& target, 
		const Bifrost::Object& source, 
		const Amino::Array<UInt32>& verOrder 
	) {
		const Amino::Array<Amino::String> vertexPropNames = BifrostExp::getGeoPropertyKeysByTarget ( source, "face_vertex_component" );
		for (int i = 0; i < vertexPropNames.size (); ++i) {
			if (target.hasProperty ( vertexPropNames[i] )) { continue; }

			auto propPtr_char = BifrostExp::getDataGeoPropValues<Amino::char_t> ( source, vertexPropNames[i] );
			if (propPtr_char) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Amino::char_t> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_short = BifrostExp::getDataGeoPropValues<Amino::short_t> ( source, vertexPropNames[i] );
			if (propPtr_short) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Amino::short_t> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_int = BifrostExp::getDataGeoPropValues<Amino::int_t> ( source, vertexPropNames[i] );
			if (propPtr_int) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Amino::int_t> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_long = BifrostExp::getDataGeoPropValues<Amino::long_t> ( source, vertexPropNames[i] );
			if (propPtr_long) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Amino::long_t> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_uchar = BifrostExp::getDataGeoPropValues<Amino::uchar_t> ( source, vertexPropNames[i] );
			if (propPtr_uchar) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Amino::uchar_t> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_ushort = BifrostExp::getDataGeoPropValues<Amino::ushort_t> ( source, vertexPropNames[i] );
			if (propPtr_ushort) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Amino::ushort_t> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_uint = BifrostExp::getDataGeoPropValues<Amino::uint_t> ( source, vertexPropNames[i] );
			if (propPtr_uint) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Amino::uint_t> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_ulong = BifrostExp::getDataGeoPropValues<Amino::ulong_t> ( source, vertexPropNames[i] );
			if (propPtr_ulong) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Amino::ulong_t> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_bool = BifrostExp::getDataGeoPropValues<Amino::bool_t> ( source, vertexPropNames[i] );
			if (propPtr_bool) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Amino::bool_t> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_String = BifrostExp::getDataGeoPropValues<Amino::String> ( source, vertexPropNames[i] );
			if (propPtr_String) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Amino::String> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_float = BifrostExp::getDataGeoPropValues<Amino::float_t> ( source, vertexPropNames[i] );
			if (propPtr_float) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Amino::float_t> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_double = BifrostExp::getDataGeoPropValues<Amino::double_t> ( source, vertexPropNames[i] );
			if (propPtr_double) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Amino::double_t> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_int2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int2> ( source, vertexPropNames[i] );
			if (propPtr_int2) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::int2> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_int2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int2x2> ( source, vertexPropNames[i] );
			if (propPtr_int2x2) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::int2x2> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_int2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int2x3> ( source, vertexPropNames[i] );
			if (propPtr_int2x3) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::int2x3> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_int3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::int3> ( source, vertexPropNames[i] );
			if (propPtr_int3) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::int3> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_float2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float2> ( source, vertexPropNames[i] );
			if (propPtr_float2) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::float2> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_float4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float4> ( source, vertexPropNames[i] );
			if (propPtr_float4) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::float4> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_float3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3> ( source, vertexPropNames[i] );
			if (propPtr_float3) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::float3> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_float2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float2x2> ( source, vertexPropNames[i] );
			if (propPtr_float2x2) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::float2x2> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_float2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float2x3> ( source, vertexPropNames[i] );
			if (propPtr_float2x3) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::float2x3> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_float3x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float3x3> ( source, vertexPropNames[i] );
			if (propPtr_float3x3) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::float3x3> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_float4x4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::float4x4> ( source, vertexPropNames[i] );
			if (propPtr_float4x4) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::float4x4> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_double2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double2> ( source, vertexPropNames[i] );
			if (propPtr_double2) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::double2> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_double3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double3> ( source, vertexPropNames[i] );
			if (propPtr_double3) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::double3> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_double4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double4> ( source, vertexPropNames[i] );
			if (propPtr_double4) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::double4> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_double2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double2x2> ( source, vertexPropNames[i] );
			if (propPtr_double2x2) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::double2x2> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_double2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double2x3> ( source, vertexPropNames[i] );
			if (propPtr_double2x3) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::double2x3> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_double3x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double3x3> ( source, vertexPropNames[i] );
			if (propPtr_double3x3) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::double3x3> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_double4x4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::double4x4> ( source, vertexPropNames[i] );
			if (propPtr_double4x4) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::double4x4> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_short2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short2> ( source, vertexPropNames[i] );
			if (propPtr_short2) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::short2> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_short2x2 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short2x2> ( source, vertexPropNames[i] );
			if (propPtr_short2x2) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::short2x2> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_short2x3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short2x3> ( source, vertexPropNames[i] );
			if (propPtr_short2x3) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::short2x3> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_short3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::short3> ( source, vertexPropNames[i] );
			if (propPtr_short3) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::short3> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_uchar3 = BifrostExp::getDataGeoPropValues<Bifrost::Math::uchar3> ( source, vertexPropNames[i] );
			if (propPtr_uchar3) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::uchar3> ( target, source, vertexPropNames[i], verOrder ); continue; }

			auto propPtr_uchar4 = BifrostExp::getDataGeoPropValues<Bifrost::Math::uchar4> ( source, vertexPropNames[i] );
			if (propPtr_uchar4) { TKCM::Dissolve::Mesh::CreateVertexCompPropCore<Bifrost::Math::uchar4> ( target, source, vertexPropNames[i], verOrder ); continue; }
		}
	}

	template<class T>
	void CreateVertexCompPropCore (
		Bifrost::Object& target,
		const Bifrost::Object& source,
		const Amino::String& propName,
		const Amino::Array<UInt32>& verOrder )
	{
		if (target.hasProperty ( propName )) { return; }

		auto propPtr = BifrostExp::getDataGeoPropValues<T> ( source, propName );
		size_t count = std::min ( propPtr->size (), verOrder.size () );

		Amino::Array<T> newPropVal ( count );
		BIF::parallel_for ( tbb::blocked_range<size_t> ( 0, count ), [&]( const tbb::blocked_range<size_t>& range ) {
			for (size_t i = range.begin (); i != range.end (); ++i) {
				newPropVal[i] = propPtr->at ( verOrder[i] );
			}
		} );

		auto newProp = Amino::newClassPtr<Amino::Array<T>> ( newPropVal );
		Amino::MutablePtr<Bifrost::Object> propObj = Bifrost::createObject ();
		BifrostExp::populateDataGeoProp<T> ( T (), std::move ( newProp ), "face_vertex_component", *propObj );

		target.setProperty ( propName, std::move ( propObj ) );
	}
}}}