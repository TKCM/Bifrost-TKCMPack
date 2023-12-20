#include "InstantMeshes_nodedef.h"

void TKCM::InstantMeshes (
	const Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
	const Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
	const Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
	const int& remeshAs,
	const int& type,
	const unsigned int& vertexCount,
	const unsigned int& faceCount,
	const float& scaleVal,
	const bool& pureQuad,
	const unsigned int& smoothIter,
	const bool& deterministic,
	const bool& extrinsic,
	const bool& alignToBoundaries,
	Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
	Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
	Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
	bool& success
) {
	point_position = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float3>>();
	face_vertex = Amino::newMutablePtr<Amino::Array<uInt>>();
	face_offset = Amino::newMutablePtr<Amino::Array<uInt>>();

	// viewer.h:72
	Float scale = -1;
	int face_count = -1;
	int vertex_count = -1;
	int knn_points = 10;
	Float creaseAngle = -1;

	// viewer.cpp:57
	int rosy, posy;
	switch (remeshAs) {
		case 0:// triangles_6x6:
			rosy = 6;
			posy = 3;
			break;
		case 1:// quads_2x4:
			rosy = 2;
			posy = 4;
			break;
		case 2:// quads_4x4:
			rosy = 4;
			posy = 4;
			break;
	}

	switch (type) {
		case 0:// vertex_count:
			vertex_count = int( vertexCount );
			break;
		case 1:// face_count:
			face_count = int( faceCount );
			break;
		case 2:// face_scale:
			scale = int( scaleVal );
			break;
	}

	// batch.cpp:49
	try {
		MatrixXu F;
		MatrixXf N;
		MatrixXf V;
		VectorXf A;
		std::set<uint32_t> crease_in, crease_out;
		BVH* bvh = nullptr;
		AdjacencyMatrix adj = nullptr;

		// convert Bifrost mesh structure to InstantMeshes structure data
		TKCM::BifMeshToInstantMeshesStructure( F, N, V, source_point_position, source_face_vertex, source_face_offset);

		bool pointcloud = F.size() == 0;
		MeshStats stats = compute_mesh_stats ( F, V, deterministic );
		if (pointcloud) {
			bvh = new BVH ( &F, &V, &N, stats.mAABB );
			bvh->build ();
			adj = generate_adjacency_matrix_pointcloud ( V, N, bvh, stats, knn_points, deterministic );
			A.resize ( V.cols () );
			A.setConstant ( 1.0f );
		}

		if (scale < 0 && vertex_count < 0 && face_count < 0) {
			vertex_count = V.cols () / 16;
		}

		if (scale > 0) {
			Float face_area = posy == 4 ? (scale * scale) : (std::sqrt ( 3.f ) / 4.f * scale * scale);
			face_count = stats.mSurfaceArea / face_area;
			vertex_count = posy == 4 ? face_count : (face_count / 2);
		}
		else if (face_count > 0) {
			Float face_area = stats.mSurfaceArea / face_count;
			vertex_count = posy == 4 ? face_count : (face_count / 2);
			scale = posy == 4 ? std::sqrt ( face_area ) : (2 * std::sqrt ( face_area * std::sqrt ( 1.f / 3.f ) ));
		}
		else if (vertex_count > 0) {
			face_count = posy == 4 ? vertex_count : (vertex_count * 2);
			Float face_area = stats.mSurfaceArea / face_count;
			scale = posy == 4 ? std::sqrt ( face_area ) : (2 * std::sqrt ( face_area * std::sqrt ( 1.f / 3.f ) ));
		}

		MultiResolutionHierarchy mRes;

		if (!pointcloud) {
			// Subdivide the mesh if necessary
			VectorXu V2E, E2E;
			VectorXb boundary, nonManifold;
			if (stats.mMaximumEdgeLength * 2 > scale || stats.mMaximumEdgeLength > stats.mAverageEdgeLength * 2) {
				build_dedge ( F, V, V2E, E2E, boundary, nonManifold );
				subdivide ( F, V, V2E, E2E, boundary, nonManifold, std::min ( scale / 2, (Float)stats.mAverageEdgeLength * 2 ), deterministic );
			}

			// Compute a directed edge data structure
			build_dedge ( F, V, V2E, E2E, boundary, nonManifold );

			// Compute adjacency matrix
			adj = generate_adjacency_matrix_uniform ( F, V2E, E2E, nonManifold );

			// Compute vertex/crease normals
			if (creaseAngle >= 0)
				generate_crease_normals ( F, V, V2E, E2E, boundary, nonManifold, creaseAngle, N, crease_in );
			else
				generate_smooth_normals ( F, V, V2E, E2E, nonManifold, N );

			// Compute dual vertex areas
			compute_dual_vertex_areas ( F, V, V2E, E2E, nonManifold, A );

			mRes.setE2E ( std::move ( E2E ) );
		}

		// Build multi-resolution hierarrchy
		mRes.setAdj ( std::move ( adj ) );
		mRes.setF ( std::move ( F ) );
		mRes.setV ( std::move ( V ) );
		mRes.setA ( std::move ( A ) );
		mRes.setN ( std::move ( N ) );
		mRes.setScale ( scale );
		mRes.build ( deterministic );
		mRes.resetSolution ();

		if ( alignToBoundaries && !pointcloud) {
			mRes.clearConstraints ();
			for (uint32_t i = 0; i < 3 * mRes.F ().cols (); ++i) {
				if (mRes.E2E ()[i] == INVALID) {
					uint32_t i0 = mRes.F ()(i % 3, i / 3);
					uint32_t i1 = mRes.F ()((i + 1) % 3, i / 3);
					Vector3f p0 = mRes.V ().col ( i0 ), p1 = mRes.V ().col ( i1 );
					Vector3f edge = p1 - p0;
					if (edge.squaredNorm () > 0) {
						edge.normalize ();
						mRes.CO ().col ( i0 ) = p0;
						mRes.CO ().col ( i1 ) = p1;
						mRes.CQ ().col ( i0 ) = mRes.CQ ().col ( i1 ) = edge;
						mRes.CQw ()[i0] = mRes.CQw ()[i1] = mRes.COw ()[i0] = mRes.COw ()[i1] = 1.0f;
					}
				}
			}
			mRes.propagateConstraints ( rosy, posy );
		}

		bvh = new BVH(&mRes.F(), &mRes.V(), &mRes.N(), stats.mAABB);
		bvh->build ();

		// batch.cpp:170
		Optimizer optimizer ( mRes, false );
		optimizer.setRoSy(rosy);
		optimizer.setPoSy ( posy );
		optimizer.setExtrinsic ( extrinsic );

		optimizer.optimizeOrientations ( -1 );
		optimizer.notify ();
		optimizer.wait ();
		std::map<uint32_t, uint32_t> sing;
		compute_orientation_singularities ( mRes, sing, extrinsic, rosy );

		optimizer.optimizePositions ( -1 );
		optimizer.notify ();
		optimizer.wait ();

		optimizer.shutdown ();

		MatrixXf O_extr, N_extr, Nf_extr;
		std::vector<std::vector<TaggedLink>> adj_extr;
		extract_graph ( optimizer.mRes, extrinsic, optimizer.rosy (), optimizer.posy (), adj_extr, O_extr, N_extr, crease_in, crease_out, deterministic );

		MatrixXu F_extr;
		extract_faces ( adj_extr, O_extr, N_extr, Nf_extr, F_extr, posy, optimizer.mRes.scale (), crease_out, true, pureQuad, bvh, int ( smoothIter ) );

		// convert InstantMeshes structure to InstantMeshes Bifrost mesh structure data
		TKCM::InstantMeshesToBifMeshStructure(point_position, face_vertex, face_offset, F_extr, Nf_extr, O_extr);

		if (bvh) {
			delete bvh;
		}
	} catch (const std::exception& e) {
		success = false;
		return;
	}
	
	success = true;
}