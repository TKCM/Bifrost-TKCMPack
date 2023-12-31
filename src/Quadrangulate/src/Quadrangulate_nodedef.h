#ifndef BIF_QUADRANGULATE_H
#define BIF_QUADRANGULATE_H
#include "../../common/TKCM_core.h"

// �l�p�`�����邽�߂ɍ폜����G�b�W��I�ʂ���m�[�h
// delete_mesh_component�m�[�h�ƃZ�b�g�Ŏg�p����
TKCM_DECL
void Quadrangulate (
	BIF_I Amino::Ptr<Bifrost::Object>& mesh			AMINO_ANNOTATE ( "Amino::Port name=mesh" ),
	BIF_O Amino::Ptr<Amino::Array<bool>>& deleteEdge	AMINO_ANNOTATE ( "Amino::Port name=delete_edge" )
) AMINO_ANNOTATE (
	"Amino::Node "
	"name=TKCM::Modeling_Toolbox::Internal::quadrangulate_select_edge "
	"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] " 
);

#endif // BIF_QUADRANGULATE_H
