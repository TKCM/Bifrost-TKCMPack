#ifndef BIF_QUADRANGULATE_H
#define BIF_QUADRANGULATE_H
#include "../../common/TKCM_core.h"

// 四角形化するために削除するエッジを選別するノード
// delete_mesh_componentノードとセットで使用する
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
