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

// �l�p�`�����邽�߂ɍ폜����G�b�W��I�ʂ���m�[�h
// delete_mesh_component�m�[�h�ƃZ�b�g�Ŏg�p����
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
