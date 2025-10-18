#if 1930 <= _MSC_VER // Visual Studio 2022 (MSVC 14.30 以降)
#define _ALLOW_COMPILER_AND_STL_VERSION_MISMATCH // Clangのバージョンエラーが出るのを無視する
#endif

#ifndef BIF_TKCMCONVEXHULL_H
#define BIF_TKCMCONVEXHULL_H
#include <vector>
#include <algorithm>

// Bifrost Amino
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>
#include <Amino/Core/Array.h>

// VHACD
#define ENABLE_VHACD_IMPLEMENTATION 1
#define VHACD_DISABLE_THREADING 0
#include <oss/v-hacd/v-hacd-master/include/VHACD.h>

// QuickHull
#include <oss/quickhull/quickhull-master/QuickHull.hpp>
#include <oss/quickhull/quickhull-master/QuickHull.cpp>

// TKCM
#include "../../core/TKCM_core.h"

using uInt = unsigned int;

namespace TKCM
{
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	TKCM_DECL
	void VHACD_core(
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
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Geometry::Internal::VHACD_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);

	TKCM_DECL
	void quick_hull_core(
		const Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,			
			  Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
			  Amino::MutablePtr<Amino::Array<uInt>>& face_vertex,
			  Amino::MutablePtr<Amino::Array<uInt>>& face_offset,
			  bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Geometry::Internal::quick_hull_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);

	TKCM_DECL
	void plane_quick_hull_core(
		const Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& source_point_position,
		const int plane_axis,
			  Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& point_position,
			  Amino::MutablePtr<Amino::Array<uInt>>& orig_point_index,
			  bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Geometry::Internal::plane_quick_hull_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);
}
#endif // BIF_TKCMCONVEXHULL_H
