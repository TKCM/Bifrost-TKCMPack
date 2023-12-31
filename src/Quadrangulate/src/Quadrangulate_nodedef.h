#ifndef BIF_QUADRANGULATE_H
#define BIF_QUADRANGULATE_H
// Bifrost Amino
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>
#include <Amino/Core/String.h>

// TKCM
#include "../../core/TKCM_core.h"
#include "../../core/TKCM_math.h"

// 四角形化するために削除するエッジを選別するノード
// delete_mesh_componentノードとセットで使用する
TKCM_DECL
void quadrangulate_core(
	const	Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& point_position,
	const	Amino::Ptr<Amino::Array<uInt>>& face_vertex,
	const	Amino::Ptr<Amino::Array<uInt>>& face_offset,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_face,
	const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_side,
			Amino::MutablePtr<Amino::Array<bool>>& half_edge_tag_data_to_delete
) AMINO_ANNOTATE (
	"Amino::Node "
	"name=TKCM::Modeling_Toolbox::Internal::quadrangulate_core "
	"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, false}] " 
);

#endif // BIF_QUADRANGULATE_H
