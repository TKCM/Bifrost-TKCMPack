#ifndef BIF_SNAP_H
#define BIF_SNAP_H
#include <limits>

// Bifrost Amino
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>
#include <Amino/Core/Array.h>
#include <Amino/Core/String.h>

// TKCM
#include "../../core/TKCM_core.h"
#include "../../core/TKCM_math.h"
#include "../../core/TKCM_geometry.h"
#include "../../core/TKCM_geometry.cpp"

namespace TKCM
{
	enum SnapTarget : uInt
	{
		Least_Target_Point_ID = 0,
		Closest_Target_Point = 1
	};

	TKCM_DECL
	void snap_grid(
				Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& point_position AMINO_ANNOTATE("Amino::InOut inName=in_point_position outName=out_point_position"),
		const	bool& use_tag_data,
		const	Amino::Ptr<Amino::Array<bool>>& point_tag_data_to_snap,
		const	bool& invert_tag,
		const	Bifrost::Math::uint3& grid_scale,
				Amino::MutablePtr<Amino::Array<bool>>& processed_point_tag_data,
				bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Modeling_Toolbox::Internal::snap_grid_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);

	TKCM_DECL
	void snap_geometry(
				Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& point_position AMINO_ANNOTATE("Amino::InOut inName=in_point_position outName=out_point_position"),
		const	bool& use_tag_data,
		const	Amino::Ptr<Amino::Array<bool>>& point_tag_data_to_snap,
		const	bool& invert_tag,
		const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& target_point_position,
		const	Amino::Ptr<Amino::Array<Amino::Ptr<Amino::Array<Amino::long_t>>>>& target_point_indices,
		const	bool& use_target_tag_data,
		const	Amino::Ptr<Amino::Array<bool>>& point_tag_data_of_target,
		const	bool& invert_target_tag,
		const	uInt& snap_using,
				Amino::MutablePtr<Amino::Array<bool>>& processed_point_tag_data,
				bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Modeling_Toolbox::Internal::snap_geometry_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);

	TKCM_DECL
	void snap_grid_point(
		const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& point_position,
		const	Bifrost::Math::uint3& grid_scale,
				Amino::MutablePtr<Amino::Array<Bifrost::Math::float3>>& grid_point
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=Modeling_Toolbox::Internal::snap_grid_point "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);
	
}

#endif // BIF_SNAP_H
