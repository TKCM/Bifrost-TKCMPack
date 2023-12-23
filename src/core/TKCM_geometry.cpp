#include "TKCM_geometry.h"

namespace TKCM {
	uInt GetNextFaceVertexIndex(
		const uInt& faceID,
		const uInt& vertexID,
		const Amino::Ptr<Amino::Array<uInt>>& face_offset)
	{
		if(face_offset->at(faceID+1) == vertexID+1){
			return face_offset->at(faceID);
		}else{
			return vertexID + 1;
		}
	}

}