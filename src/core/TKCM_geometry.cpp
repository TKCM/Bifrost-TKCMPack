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
	
	uInt GetPreviousFaceVertexIndex(
		const uInt& faceID,
		const uInt& vertexID,
		const Amino::Ptr<Amino::Array<uInt>>& face_offset)
	{
		if(face_offset->at(faceID) == vertexID){
			return face_offset->at(faceID+1) - 1;
		}else{
			return vertexID - 1;
		}
	}

	uInt GetEdgeLeftPolygonID(
		const Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
		const Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
		const Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_edge_face,
		const Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_edge_side,
		const Amino::Ptr<Amino::Array<uInt>>& point_face_adjacency_index,
		const uInt halfEdgeID
	){
		int poiID = source_face_vertex->at(halfEdgeID);
		for (uInt i = point_face_adjacency_index->at(poiID); i < point_face_adjacency_index->at(poiID + 1); ++i){
			uInt adjFaceID = point_face_adjacent_edge_face->at(i);
			if (source_face_offset->size() - 1 <= adjFaceID){
				return 4294967295;
			}

			uInt vertexID = source_face_offset->at(adjFaceID) + point_face_adjacent_edge_side->at(i);
			if (halfEdgeID == vertexID){
				return adjFaceID;
			}
		}
		return 4294967295;
	}

	uInt GetNextPointIndex(
		const uInt& faceID,
		const uInt& pointID,
		const Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
		const Amino::Ptr<Amino::Array<uInt>>& face_offset
	){
		if (face_offset->size() <= faceID){
			return -1;
		}

		for (uInt i = face_offset->at(faceID); i < face_offset->at(faceID + 1); ++i){
			if (source_face_vertex->at(i) == pointID){
				if (i == face_offset->at(faceID + 1) - 1){
					// pointIDがポリゴンを構成する最後のポイントだった場合
					return source_face_vertex->at(face_offset->at(faceID));
				} else{
					return source_face_vertex->at(i + 1);
				}
			}
		}
		return -1; // 無限値を返す
	}

	uInt GetPreviousPointIndex(
		const uInt& faceID,
		const uInt& pointID,
		const Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
		const Amino::Ptr<Amino::Array<uInt>>& face_offset
	){
		if (face_offset->size() <= faceID){ return -1; }

		for (uInt i = face_offset->at(faceID); i < face_offset->at(faceID + 1); ++i){
			if (source_face_vertex->at(i) == pointID){
				if (i == face_offset->at(faceID)){
					// pointIDがポリゴンを構成する最初のポイントだった場合
					return source_face_vertex->at(face_offset->at(faceID+1)-1);
				} else{
					return source_face_vertex->at(i + 1);
				}
			}
		}
		return -1; // 無限値を返す
	}

	uInt GetHalfEdge(
		const uInt& faceID,
		const uInt& pointID,
		const Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
		const Amino::Ptr<Amino::Array<uInt>>& face_offset
	){
		if (face_offset->size() <= faceID){ return -1; }

		for (uInt i = face_offset->at(faceID); i < face_offset->at(faceID + 1); ++i){
			if (source_face_vertex->at(i) == pointID){
				return i;
			}
		}
		return -1; // 無限値を返す
	}

	uInt GetLocalFaceVertexIndex(
		const uInt& faceID,
		const uInt& pointID,
		const Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
		const Amino::Ptr<Amino::Array<uInt>>& face_offset
	){
		if (face_offset->size() <= faceID){ return -1; }

		for (uInt i = face_offset->at(faceID); i < face_offset->at(faceID + 1); ++i){
			if (source_face_vertex->at(i) == pointID){
				return i - face_offset->at(faceID);
			}
		}
		return -1; // 無限値を返す
	}
}