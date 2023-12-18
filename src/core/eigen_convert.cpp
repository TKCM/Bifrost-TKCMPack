#include "eigen_convert.h"


namespace TKCM {
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Vector3f -- float3
	void Convert ( BIF_O Eigen::Vector3f &eigen, BIF_I Bifrost::Math::float3 &bif ) {
		eigen << bif.x, bif.y, bif.z;
	}

	void Convert ( BIF_O Bifrost::Math::float3 &bif, BIF_I Eigen::Vector3f &eigen ) {
		bif.x = eigen.x();
		bif.y = eigen.y();
		bif.z = eigen.z();
	}

	void Convert ( BIF_O std::vector<Eigen::Vector3f>& eigen, BIF_I Amino::Array<Bifrost::Math::float3>& bif ) {
		eigen.resize ( bif.size () );

		#pragma omp parallel for
		for (int i = 0; i < bif.size (); ++i) {
			eigen[i] << bif[i].x, bif[i].y, bif[i].z;
		}
	}

	void Convert ( BIF_O Amino::Array<Bifrost::Math::float3>& bif, BIF_I std::vector<Eigen::Vector3f>& eigen ) {
		bif.resize ( eigen.size () );

		#pragma omp parallel for
		for (int i = 0; i < bif.size (); ++i) {
			bif[i].x = eigen[i].x();
			bif[i].y = eigen[i].y();
			bif[i].z = eigen[i].z();
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Matrix4f -- float4x4
	void Convert ( BIF_O Eigen::Matrix4f &eigen, BIF_I Bifrost::Math::float4x4 &bif ) {
		eigen << bif.c0.x, bif.c1.x, bif.c2.x, bif.c3.x, bif.c0.y, bif.c1.y, bif.c2.y, bif.c3.y, bif.c0.z, bif.c1.z, bif.c2.z, bif.c3.z, bif.c0.w, bif.c1.w, bif.c2.w, bif.c3.w;
	}

	/*void Convert ( BIF_O Bifrost::Math::float4x4 &bif, BIF_I Eigen::Matrix4f &eigen ) {
		bif.c0 << eigen.col ( 0 ).x, eigen.col ( 0 ).y, eigen.col ( 0 ).z, eigen.col ( 0 ).w;
		bif.c1 << eigen.col ( 1 ).x, eigen.col ( 1 ).y, eigen.col ( 1 ).z, eigen.col ( 1 ).w;
		bif.c2 << eigen.col ( 2 ).x, eigen.col ( 2 ).y, eigen.col ( 2 ).z, eigen.col ( 2 ).w;
		bif.c3 << eigen.col ( 3 ).x, eigen.col ( 3 ).y, eigen.col ( 3 ).z, eigen.col ( 3 ).w;
	}*/

	void Convert ( BIF_O std::vector<Eigen::Matrix4f> &eigen, BIF_I Amino::Array<Bifrost::Math::float4x4> &bif ) {
		eigen.resize ( bif.size () );

		#pragma omp parallel for
		for (int i = 0; i < bif.size (); ++i) {
			eigen[i] << bif[i].c0.x, bif[i].c1.x, bif[i].c2.x, bif[i].c3.x, bif[i].c0.y, bif[i].c1.y, bif[i].c2.y, bif[i].c3.y, bif[i].c0.z, bif[i].c1.z, bif[i].c2.z, bif[i].c3.z, bif[i].c0.w, bif[i].c1.w, bif[i].c2.w, bif[i].c3.w;
		}
	}

	/*void Convert ( BIF_O Amino::Array<Bifrost::Math::float4x4> &bif, BIF_I std::vector<Eigen::Matrix4f> &eigen ) {
		bif.resize ( eigen.size () );

		#pragma omp parallel for
		for (int i = 0; i < eigen.size (); ++i) {
			bif[i].c0 << eigen[i].col ( 0 ).x, eigen[i].col ( 0 ).y, eigen[i].col ( 0 ).z, eigen[i].col ( 0 ).w;
			bif[i].c1 << eigen[i].col ( 1 ).x, eigen[i].col ( 1 ).y, eigen[i].col ( 1 ).z, eigen[i].col ( 1 ).w;
			bif[i].c2 << eigen[i].col ( 2 ).x, eigen[i].col ( 2 ).y, eigen[i].col ( 2 ).z, eigen[i].col ( 2 ).w;
			bif[i].c3 << eigen[i].col ( 3 ).x, eigen[i].col ( 3 ).y, eigen[i].col ( 3 ).z, eigen[i].col ( 3 ).w;
		}
	}*/

}