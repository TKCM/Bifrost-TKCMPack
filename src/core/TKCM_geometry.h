#pragma once
#include <vector>
#include <sstream>
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>

using uInt = unsigned int;

namespace TKCM {
	uInt GetNextFaceVertexIndex(
		const uInt& faceID,
		const uInt& vertexID,
		const Amino::Ptr<Amino::Array<uInt>>& face_offset);
		
	uInt GetPreviousFaceVertexIndex(
		const uInt& faceID,
		const uInt& vertexID,
		const Amino::Ptr<Amino::Array<uInt>>& face_offset);

	uInt GetEdgeLeftPolygonID(
		const Amino::Ptr<Amino::Array<uInt>>& source_face_vertex,
		const Amino::Ptr<Amino::Array<uInt>>& source_face_offset,
		const Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_edge_face,
		const Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_edge_side,
		const Amino::Ptr<Amino::Array<uInt>>& point_face_adjacency_index,
		const uInt halfEdgeID);
}