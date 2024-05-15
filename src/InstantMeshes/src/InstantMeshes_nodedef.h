#ifndef BIF_INSTANTMESHES_H
#define BIF_INSTANTMESHES_H

// Bifrost Amino
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>
#include <Amino/Core/Array.h>

// TKCM
#include "../../core/TKCM_core.h"
#include "InstantMeshes_functions.h"

namespace TKCM
{
	TKCM_DECL
	void InstantMeshes(
		const Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
		const Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
		const Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
		const int& remeshAs						AMINO_ANNOTATE("Amino::Port name=remesh_as"),
		const int& type							AMINO_ANNOTATE("Amino::Port name=target_type"),
		const unsigned int& vertexCount			AMINO_ANNOTATE("Amino::Port name=vertex_count"),
		const unsigned int& faceCount			AMINO_ANNOTATE("Amino::Port name=face_count"),
		const float& scaleVal					AMINO_ANNOTATE("Amino::Port name=face_scale"),
		const bool& pureQuad					AMINO_ANNOTATE("Amino::Port name=pure_quad"),
		const unsigned int& smoothIter			AMINO_ANNOTATE("Amino::Port name=smooth_iter"),
		const bool& deterministic				AMINO_ANNOTATE("Amino::Port name=deterministic"),
		const bool& extrinsic					AMINO_ANNOTATE("Amino::Port name=extrinsic"),
		const bool& alignToBoundaries			AMINO_ANNOTATE("Amino::Port name=align_to_boundaries"),
		Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
		Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
		Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
		bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Modeling_Toolbox::Internal::InstantMeshes_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);
}
#endif // BIF_INSTANTMESHES_H