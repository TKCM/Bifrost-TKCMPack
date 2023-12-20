#include "InstantMeshes_functions.h"

namespace TKCM{
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void BifMeshToInstantMeshesStructure(
		MatrixXu& F,
		MatrixXf& N,
		MatrixXf& V,
		const Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& point_position,
		const Amino::Ptr<Amino::Array<uInt>>& face_vertex,
		const Amino::Ptr<Amino::Array<uInt>>& face_offset )
	{
		if (!point_position || point_position->empty()){ return; }
		if (!face_vertex || face_vertex->empty()){ return; }
		if (!face_offset || face_offset->empty()){ return; }

		// point position
		int pointCount = point_position->size ();
		V.resize ( 3, pointCount );
		#pragma omp parallel for
		for (int i = 0; i < pointCount; ++i) {
			V ( 0, i ) = point_position->at(i).x;
			V ( 1, i ) = point_position->at(i).y;
			V ( 2, i ) = point_position->at(i).z;
		}

		// point normal
		N.resize ( 3, pointCount );
		
		// polygon size
		int polygonCount = face_offset->size () - 1;
		
		// packed polygon point indices
		std::vector<unsigned int> ids;
		ids.reserve ( size_t(3 * polygonCount) );
		for (int i = 0; i < polygonCount; ++i) {
			int polygonSize = face_offset->at ( i + 1 ) - face_offset->at ( i );
			int triangleCount = polygonSize - 2;
			int firstPointID = face_offset->at ( i );
			for (int j = 0; j < triangleCount; ++j) {
				ids.push_back ( face_vertex->at ( firstPointID ) );
				ids.push_back ( face_vertex->at ( size_t(firstPointID + j + 1 )) );
				ids.push_back ( face_vertex->at ( size_t(firstPointID + j + 2 )) );
			}
		}
		F = Eigen::Map<MatrixXu> ( ids.data(), 3, ids.size() / 3 );
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void InstantMeshesToBifMeshStructure(
		Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
		Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
		Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
		const MatrixXu& F, 
		const MatrixXf& N, 
		const MatrixXf& V)
	{
		// polygon info
		int polySizeSum = 0;
		face_offset->push_back ( polySizeSum );

		for ( uint32_t f = 0; f < F.cols(); ++f ){
			if ( F.rows() == 4 && F( 2, f ) == F( 3, f ) ){ continue; }

			int nbF = 0;
			for ( uint32_t i = 0; i < F.rows(); ++i ){
				face_vertex->push_back ( F( i, f ) );
				nbF++;
			}
			if ( 3 <= nbF ){
				polySizeSum += nbF;
				face_offset->push_back ( polySizeSum );
			}
		}

		// point position
		point_position->resize(V.cols());
		#pragma omp parallel for
		for ( int i = 0; i < V.cols(); ++i ){
			point_position->at(i).x = V( 0, i );
			point_position->at(i).y = V( 1, i );
			point_position->at(i).z = V( 2, i );
		}

		std::map<uint32_t, std::pair<uint32_t, std::map<uint32_t, uint32_t>>> irregular;
		size_t nIrregular = 0;
		if ( F.rows() == 4 ){
			for ( uint32_t f = 0; f < F.cols(); ++f ){
				if ( F( 2, f ) == F( 3, f ) ){
					nIrregular++;
					auto& value = irregular[F( 2, f )];
					value.first = f;
					value.second[F( 0, f )] = F( 1, f );
				}
			}
		}
		for ( auto item : irregular ){
			auto face = item.second;
			uint32_t v = face.second.begin()->first, first = v, i = 0;
			int nbF = 0;
			while ( true ){
				face_vertex->push_back ( v );
				nbF++;
				v = face.second[v];
				i++;

				if ( v == first || i == face.second.size() ){ break; }
			}
			if ( 3 <= nbF ){
				polySizeSum += nbF;
				face_offset->push_back ( polySizeSum );
			}
			nbF = 0;
			while ( i != face.second.size() ){
				face_vertex->push_back ( v );
				nbF++;
				i++;
			}
			if ( 3 <= nbF ){
				polySizeSum += nbF;
				face_offset->push_back ( polySizeSum );
			}
		}
	}
}