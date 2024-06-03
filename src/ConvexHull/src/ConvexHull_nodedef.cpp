#include "ConvexHull_nodedef.h"
#include "ConvexHull_functions.h"

void TKCM::VHACD_core(
	const Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
	const Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
	const Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
	const uint32_t max_convex_hulls,
	const uint32_t resolution,
	const double   minimum_volume_percent_error_allowed,
	const uint32_t max_recursion_depth,
	const bool     shrink_wrap,
	const uint32_t fill_mode,
	const uint32_t max_num_vertices_per_cnvex_hull,
	const bool     async_ACD,
	const uint32_t min_edge_length,
	const bool     find_best_plane,
		  Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
		  Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
		  Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
		  bool& success
){
	// 出力データの準備
	point_position = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float3>>();
	face_vertex = Amino::newMutablePtr<Amino::Array<uInt>>();
	face_offset = Amino::newMutablePtr<Amino::Array<uInt>>();
	success = false;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 入力データのチェック
	if (!source_point_position || source_point_position->empty() ||
		!source_face_vertex || source_face_vertex->empty() ||
		!source_face_offset || source_face_offset->empty()
	){ 
		return; 
	}

	const uint32_t poiCount = source_point_position->size();
	const uint32_t polyCount = source_face_offset->size() - 1;
	const uint32_t verCount = source_face_vertex->size();

	if (poiCount < 4 || polyCount < 2 || verCount < 6){ return; }

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// V-HACD

	VHACD::IVHACD::Parameters params;
	params.m_maxConvexHulls = max_convex_hulls;
	params.m_resolution = resolution;
	params.m_minimumVolumePercentErrorAllowed = minimum_volume_percent_error_allowed;
	params.m_maxRecursionDepth = max_recursion_depth;
	params.m_shrinkWrap = shrink_wrap;
	params.m_maxNumVerticesPerCH = max_num_vertices_per_cnvex_hull;
	params.m_asyncACD = async_ACD;
	params.m_minEdgeLength = min_edge_length;
	params.m_findBestPlane = find_best_plane;
	params.m_fillMode = static_cast<VHACD::FillMode>(fill_mode);

	VHACD::IVHACD* vhacd = VHACD::CreateVHACD();
	bool ret = vhacd->Compute(
		(const float* const)&source_point_position->front().x,
		poiCount,
		(const uint32_t* const)&source_face_vertex->front(),
		polyCount,
		params
	);	

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 結果をメッシュトポロジーデータとして出力
	if (ret){
		point_position->reserve(max_convex_hulls * max_num_vertices_per_cnvex_hull);
		face_vertex->reserve(max_convex_hulls * max_num_vertices_per_cnvex_hull * 3);
		face_offset->reserve(max_convex_hulls * max_num_vertices_per_cnvex_hull);
		face_offset->push_back(0);

		uint32_t nConvexHulls = vhacd->GetNConvexHulls();
		VHACD::IVHACD::ConvexHull cHullCell;
		// ConvexHullのデータを順にコンバートしていく
		for (uint32_t cHull_i = 0; cHull_i < nConvexHulls; ++cHull_i){
			int startID = point_position->size();
			
			vhacd->GetConvexHull(cHull_i, cHullCell);
			
			// 頂点位置
			const std::vector<VHACD::Vertex> &poiPos = cHullCell.m_points;
			for (int i = 0; i < poiPos.size(); ++i){ 
				point_position->push_back(
					Bifrost::Math::float3{ float(poiPos[i].mX),float(poiPos[i].mY),float(poiPos[i].mZ) }
				);
			}

			// 頂点/ポリゴンリスト
			const std::vector<VHACD::Triangle> &tri = cHullCell.m_triangles;
			for (int i = 0; i < tri.size(); ++i){
				face_vertex->push_back(startID + tri[i].mI0);
				face_vertex->push_back(startID + tri[i].mI1);
				face_vertex->push_back(startID + tri[i].mI2);

				face_offset->push_back(face_offset->back() + 3);
			}
		}

		point_position->shrink_to_fit();
		face_vertex->shrink_to_fit();
		face_offset->shrink_to_fit();
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	vhacd->Clean();
	vhacd->Release();
	success = ret;
}

void TKCM::quick_hull_core(
	const Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
		  Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
		  Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
		  Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
		  bool& success
){
	// 出力データの準備
	point_position = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float3>>();
	face_vertex = Amino::newMutablePtr<Amino::Array<uInt>>();
	face_offset = Amino::newMutablePtr<Amino::Array<uInt>>();
	success = false;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 入力データのチェック
	if (!source_point_position || source_point_position->empty()){ return; }
	const uint32_t poiCount = source_point_position->size();

	if (poiCount < 4 ){ return; }

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// QuickHull

	quickhull::QuickHull<float> qh;
	quickhull::ConvexHull<float> hull = qh.getConvexHull(&source_point_position->front().x, poiCount, false, false);
	const std::vector<size_t>& indexBuffer = hull.getIndexBuffer();
	const quickhull::VertexDataSource<float>& vertexBuffer = hull.getVertexBuffer();

	face_vertex->assign(indexBuffer.begin(), indexBuffer.end());

	face_offset->resize(indexBuffer.size() / 3 + 1 );
	face_offset->at(0) = 0;
	#pragma omp parallel for
	for (int i = 0; i < indexBuffer.size() / 3; ++i){
		face_offset->at(i + 1) = (i+1) * 3;
	}
	
	point_position->resize(vertexBuffer.size());
	#pragma omp parallel for
	for (int i = 0; i < vertexBuffer.size(); ++i){
		point_position->at(i) =	Bifrost::Math::float3{ vertexBuffer[i].x, vertexBuffer[i].y, vertexBuffer[i].z };
	}	

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	success = true;
}

void TKCM::plane_quick_hull_core(
	const Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
	const int plane_axis,
		  Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
		  Amino::MutablePtr<Amino::Array<uInt>>& orig_point_index,
		  bool& success
){
	// 出力データの準備
	point_position = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float3>>();
	orig_point_index = Amino::newMutablePtr<Amino::Array<uInt>>();
	success = false;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// 入力データのチェック
	if (!source_point_position || source_point_position->empty()){ return; }
	const uint32_t poiCount = source_point_position->size();

	if (poiCount < 4){ return; }

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// QuickHull-plane
	Amino::Array<uInt> point_index = TKCM::convexHullPlane(*source_point_position, plane_axis);

	point_position->reserve(point_index.size());
	orig_point_index->reserve(point_index.size());
	for (int i = 0; i < point_index.size(); ++i){
		Bifrost::Math::float3 p = source_point_position->at(point_index.at(i));
		switch (plane_axis){
			case 0: p.y = 0.0f; break; // x,z
			case 1: p.x = 0.0f; break; // y,z
			case 2: p.z = 0.0f; break; // x,y
		}
		point_position->push_back(p);
		orig_point_index->push_back(point_index[i]);
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	success = true;
}