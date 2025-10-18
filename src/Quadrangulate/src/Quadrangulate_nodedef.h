#if 1930 <= _MSC_VER // Visual Studio 2022 (MSVC 14.30 以降)
#define _ALLOW_COMPILER_AND_STL_VERSION_MISMATCH // Clangのバージョンエラーが出るのを無視する
#endif

#ifndef BIF_QUADRANGULATE_H
#define BIF_QUADRANGULATE_H
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

// 四角形化するために削除するエッジを選別するノード
// delete_mesh_componentノードとセットで使用する
namespace TKCM {
	TKCM_DECL
	void quadrangulate_core(
		const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& point_position,
		const	Amino::Ptr<Amino::Array<uInt>>& face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& face_offset,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_face,
		const	Amino::Ptr<Amino::Array<uInt>>& face_ver_adjacent_edge_side,
				Amino::MutablePtr<Amino::Array<bool>>& half_edge_tag_data_to_dissolve
	) AMINO_ANNOTATE (
		"Amino::Node "
		"name=Modeling_Toolbox::Internal::quadrangulate_core "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] " 
	);
}
#endif // BIF_QUADRANGULATE_H
