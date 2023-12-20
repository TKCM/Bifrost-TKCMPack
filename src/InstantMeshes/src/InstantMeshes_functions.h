#pragma once
#include <vector>

// Bifrost Amino
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>

// InstantMeshes
#include <src/batch.h>
#include <src/meshio.h>
#include <src/dedge.h>
#include <src/subdivide.h>
#include <src/meshstats.h>
#include <src/hierarchy.h> 
#include <src/field.h>
#include <src/normal.h>
#include <src/extract.h>
#include <src/bvh.h>
#include <src/smoothcurve.h>
#include <src/extract.h>

using uInt = unsigned int;

namespace TKCM{
	void BifMeshToInstantMeshesStructure(
		MatrixXu& F, 
		MatrixXf& N, 
		MatrixXf& V, 
		const Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& point_position, 
		const Amino::Ptr<Amino::Array<uInt>>& face_vertex, 
		const Amino::Ptr<Amino::Array<uInt>>& face_offset);

	void InstantMeshesToBifMeshStructure(
		Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
		Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
		Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
		const MatrixXu& F, 
		const MatrixXf& N, 
		const MatrixXf& V );

}