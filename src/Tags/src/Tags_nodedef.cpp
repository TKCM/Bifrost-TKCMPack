#include "Tags_nodedef.h"
#include <string>

namespace TKCM
{
	TKCM::ComponentType ToComponentType(const Amino::String& component_type){
		std::string v("face_vertex_component");
		std::string p("point_component");
		std::string f("face_component");
		std::string e("half_edge_vertex_component");
		if (component_type == Amino::String(v.c_str(), v.size())){
			return TKCM::ComponentType::face_vertex_component;
		} else if (component_type == Amino::String(p.c_str(), p.size())){
			return TKCM::ComponentType::point_component;
		} else if (component_type == Amino::String(f.c_str(), f.size())){
			return TKCM::ComponentType::face_component;
		} else if (component_type == Amino::String(e.c_str(), e.size())){
			return TKCM::ComponentType::half_edge_vertex_component;
		} else{
			return TKCM::ComponentType::non;
		}
	}

	bool NewTagPtr(
				Amino::MutablePtr<Amino::Array<bool>>& newData,
		const	uInt& point_count,
		const	Amino::Ptr<Amino::Array<uInt>>& face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& face_offset,
		const	Amino::String& component_type
	){
		newData = Amino::newMutablePtr<Amino::Array<bool>>();
		if (!face_vertex || face_vertex->empty() || !face_offset || face_offset->empty()){ return false; }		

		uInt vertexCount = face_vertex->size();
		uInt polyCount = face_offset->size() - 1;

		TKCM::ComponentType comp_type = TKCM::ToComponentType(component_type);
		switch (comp_type){
			case TKCM::ComponentType::face_vertex_component:		newData->resize(vertexCount, false); break;
			case TKCM::ComponentType::point_component:				newData->resize(point_count, false); break;
			case TKCM::ComponentType::face_component:				newData->resize(polyCount, false); break;
			case TKCM::ComponentType::half_edge_vertex_component:	newData->resize(vertexCount, false); break;
			case TKCM::ComponentType::non:							return false;
		}
		return true;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	void point_tag_promote(
		const	uInt& point_count,
		const	Amino::Ptr<Amino::Array<uInt>>& face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& face_offset,
		const	Amino::Ptr<Amino::Array<bool>>& in_tag,
		const	Amino::String& out_component_type,
		const	uInt& logical_operation,
				Amino::MutablePtr<Amino::Array<bool>>& out_tag,
				bool& success
	){
		// �o�̓f�[�^����������
		if (TKCM::NewTagPtr(out_tag, point_count, face_vertex, face_offset, out_component_type) == false){ success = false; return; }

		// ���̓f�[�^�ɕs�����Ȃ����m�F����
		uInt vertexCount = face_vertex->size();
		uInt polyCount = face_offset->size() - 1;
		if (!in_tag || in_tag->empty() || in_tag->size() != point_count){ success = false; return; }

		TKCM::ComponentType out_type = TKCM::ToComponentType(out_component_type);
		switch (out_type){
			case TKCM::ComponentType::face_vertex_component:
			{
				#pragma omp parallel for
				for (size_t polyId = 0; polyId < polyCount; ++polyId){
					for (uInt verId = face_offset->at(polyId); verId < face_offset->at(polyId + 1); ++verId){
						uInt verPoiId = face_vertex->at(verId);
						out_tag->at(verId) = in_tag->at(verPoiId);
					}
				}
				break;
			}
			case TKCM::ComponentType::point_component:
			{
				#pragma omp parallel for
				for (size_t i = 0; i < in_tag->size(); ++i){
					out_tag->at(i) = in_tag->at(i);
				}
				break;
			}
			case TKCM::ComponentType::face_component:
			{
				#pragma omp parallel for
				for (size_t polyId = 0; polyId < polyCount; ++polyId){
					switch (logical_operation){
						case LogicalOperation::AND: // �|���S�����\������S���_��true�ł���΃^�O�t������
						{
							bool allTrue = true;
							for (uInt verId = face_offset->at(polyId); verId < face_offset->at(polyId + 1); ++verId){
								uInt verPoiId = face_vertex->at(verId);
								if (in_tag->at(verPoiId) == false){
									allTrue = false;
									break;
								}
							}
							out_tag->at(polyId) = allTrue;
							break;
						}
						case LogicalOperation::OR: // �|���S�����\�����钸�_�̂ǂꂩ�P�ł�true�ł���΃^�O�t������
						{
							for (uInt verId = face_offset->at(polyId); verId < face_offset->at(polyId + 1); ++verId){
								uInt verPoiId = face_vertex->at(verId);
								if (in_tag->at(verPoiId) == true){
									out_tag->at(polyId) = true;
									break;
								}
							}
							break;
						}
					}
				}
				break;
			}
			case TKCM::ComponentType::half_edge_vertex_component:
			{
				#pragma omp parallel for
				for (size_t polyId = 0; polyId < polyCount; ++polyId){
					for (uInt verId = face_offset->at(polyId); verId < face_offset->at(polyId + 1); ++verId){
						// �|���S�����\������n�[�t�G�b�W�̎n�_�E�I�_�̒��_�ԍ����擾����
						uInt edgeStartPoiId = face_vertex->at(verId);
						uInt nextFaceVerId = TKCM::GetNextFaceVertexIndex(polyId, verId, face_offset);
						uInt edgeEndPoiId = face_vertex->at(nextFaceVerId);
						switch (logical_operation){ // �G�b�W�̎n�_/�I�_�̗����̃|�C���g��true�ł���΃^�O�t������
							case LogicalOperation::AND:
								out_tag->at(verId) = in_tag->at(edgeStartPoiId) && in_tag->at(edgeEndPoiId);
								break;
							case LogicalOperation::OR: // �G�b�W�̎n�_/�I�_�̂ǂ��炩�̃|�C���g��true�ł���΃^�O�t������
								out_tag->at(verId) = in_tag->at(edgeStartPoiId) || in_tag->at(edgeEndPoiId);
								break;
						}
					}
				}
				break;
			}
		}
		success = true;
	}

	void vertex_tag_promote(
		const	uInt& point_count,
		const	Amino::Ptr<Amino::Array<uInt>>& face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& face_offset,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_face,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_side,
		const	Amino::Ptr<Amino::Array<uInt>>& point_face_adjacent_index,
		const	Amino::Ptr<Amino::Array<bool>>& in_tag,
		const	Amino::String& out_component_type,
		const	uInt& logical_operation,
				Amino::MutablePtr<Amino::Array<bool>>& out_tag,
				bool& success
	){
		// �o�̓f�[�^����������
		if (TKCM::NewTagPtr(out_tag, point_count, face_vertex, face_offset, out_component_type) == false){ success = false; return; }

		// ���̓f�[�^�ɕs�����Ȃ����m�F����
		uInt vertexCount = face_vertex->size();
		uInt polyCount = face_offset->size() - 1;
		if (!in_tag || in_tag->empty() || in_tag->size() != vertexCount ||
			!point_face_adjacent_face || point_face_adjacent_face->empty() ||
			!point_face_adjacent_side || point_face_adjacent_side->empty() ||
			!point_face_adjacent_index || point_face_adjacent_index->empty()
			){ success = false; return; }

		TKCM::ComponentType out_type = TKCM::ToComponentType(out_component_type);
		switch (out_type){
			case TKCM::ComponentType::face_vertex_component:
			{
				#pragma omp parallel for
				for (size_t i = 0; i < in_tag->size(); ++i){
					out_tag->at(i) = in_tag->at(i);
				}
				break;
			}
			case TKCM::ComponentType::point_component:
			{
				switch (logical_operation){
					case TKCM::LogicalOperation::AND: // �|�C���g���\������S�Ẵo�[�e�b�N�X��true�ł���΃^�O�t������
					{
						#pragma omp parallel for
						for (size_t poiID = 0; poiID < point_count; ++poiID){
							bool allTrue = true;
							for (int i = point_face_adjacent_index->at(poiID); i < point_face_adjacent_index->at(poiID + 1); ++i){
								uInt adjFaceID = point_face_adjacent_face->at(i);
								uInt adjVerID = face_offset->at(adjFaceID) + point_face_adjacent_side->at(i);
								if (in_tag->at(adjVerID) == false){
									allTrue = false;
									break;
								}
							}
							out_tag->at(poiID) = allTrue;
						}						
						break;
					}
					case TKCM::LogicalOperation::OR: // �|�C���g���\������o�[�e�b�N�X�̂ǂꂩ�P�ł�true�ł���΃^�O�t������
					{
						#pragma omp parallel for
						for (size_t verID = 0; verID < vertexCount; ++verID){
							if (in_tag->at(verID) == true){
								uInt poiID = face_vertex->at(verID);
								out_tag->at(poiID) = true;
							}
						}
						break;
					}
				}
				break;
			}
			case TKCM::ComponentType::face_component :
			{
				switch (logical_operation){
					case TKCM::LogicalOperation::AND: // �|���S�����\������S�o�[�e�b�N�X��true�ł���΃^�O�t������
					{
						#pragma omp parallel for
						for (size_t polyId = 0; polyId < polyCount; ++polyId){
							bool allTrue = true;
							for (uInt verId = face_offset->at(polyId); verId < face_offset->at(polyId + 1); ++verId){
								if (in_tag->at(verId) == false){
									allTrue = false;
									break;
								}
							}
							out_tag->at(polyId) = allTrue;
						}
						break;
					}
					case TKCM::LogicalOperation::OR: // �|���S�����\������o�[�e�b�N�X�̂ǂꂩ�P�ł�true�ł���΃^�O�t������
					{
						#pragma omp parallel for
						for (size_t polyId = 0; polyId < polyCount; ++polyId){
							for (uInt verId = face_offset->at(polyId); verId < face_offset->at(polyId + 1); ++verId){
								if (in_tag->at(verId) == true){
									out_tag->at(polyId) = true;
									break;
								}
							}
						}
						break;
					}
				}
				break;
			}
			case TKCM::ComponentType::half_edge_vertex_component:
			{
				#pragma omp parallel for
				for (size_t polyId = 0; polyId < polyCount; ++polyId){
					for (uInt verId = face_offset->at(polyId); verId < face_offset->at(polyId + 1); ++verId){
						// �n�[�t�G�b�W�̎n�_�E�I�_�̃o�[�e�b�N�X�ԍ����擾����
						uInt edgeStartVerId = verId;
						uInt edgeEndVerId = TKCM::GetNextFaceVertexIndex(polyId, edgeStartVerId, face_offset);
						
						switch (logical_operation){
							case TKCM::LogicalOperation::AND: // �G�b�W�̎n�_�ƏI�_�̗����̃o�[�e�b�N�X��true�ł���΃^�O�t������
								out_tag->at(verId) = in_tag->at(edgeStartVerId) && in_tag->at(edgeEndVerId);
								break;
							case TKCM::LogicalOperation::OR: // �G�b�W�̎n�_�ƏI�_�̂ǂ��炩�̃o�[�e�b�N�X��true�ł���΃^�O�t������
								out_tag->at(verId) = in_tag->at(edgeStartVerId) || in_tag->at(edgeEndVerId);
								break;
						}
					}
				}
				break;
			}
		}
		success = true;
	}

	void face_tag_promote(
		const	uInt& point_count,
		const	Amino::Ptr<Amino::Array<uInt>>& face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& face_offset,
		const	Amino::Ptr<Amino::Array<bool>>& in_tag,
		const	Amino::String& out_component_type,
				Amino::MutablePtr<Amino::Array<bool>>& out_tag,
				bool& success
	){
		// �o�̓f�[�^����������
		if (TKCM::NewTagPtr(out_tag, point_count, face_vertex, face_offset, out_component_type) == false ){ success = false; return; }
		
		// ���̓f�[�^�ɕs�����Ȃ����m�F����
		uInt vertexCount = face_vertex->size();
		uInt polyCount = face_offset->size() - 1;
		if (!in_tag || in_tag->empty() || in_tag->size() != polyCount){ success = false; return; }

		TKCM::ComponentType out_type = TKCM::ToComponentType(out_component_type);
		switch (out_type){
			case TKCM::ComponentType::face_vertex_component:
			case TKCM::ComponentType::half_edge_vertex_component:
			{
				#pragma omp parallel for
				for (size_t polyId = 0; polyId < polyCount; ++polyId){
					if (in_tag->at(polyId) == true){
						for (uInt verId = face_offset->at(polyId); verId < face_offset->at(polyId + 1); ++verId){
							out_tag->at(verId) = true;
						}
					}
				}
				break;
			}
			case TKCM::ComponentType::point_component:
			{
				#pragma omp parallel for
				for (size_t polyId = 0; polyId < polyCount; ++polyId){
					if (in_tag->at(polyId) == true){
						for (uInt verId = face_offset->at(polyId); verId < face_offset->at(polyId + 1); ++verId){
							uInt verPoiId = face_vertex->at(verId);
							out_tag->at(verPoiId) = true;
						}
					}
				}
				break;
			}
			case TKCM::ComponentType::face_component:
			{
				#pragma omp parallel for
				for (size_t i = 0; i < in_tag->size(); ++i){
					out_tag->at(i) = in_tag->at(i);
				}
				break;
			}
		}
		success = true;
	}

	void to_face_vertex(
				Amino::Ptr<Amino::Array<bool>>& tag_data,
		const	Amino::Ptr<Amino::Array<uInt>>& face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& face_offset,
				bool& success
	){
		Amino::MutablePtr<Amino::Array<bool>> rest = Amino::newMutablePtr<Amino::Array<bool>>(tag_data->size(), false);

		#pragma omp parallel for
		for (size_t polyId = 0; polyId < face_offset->size() - 1; ++polyId){
			for (uInt verId = face_offset->at(polyId); verId < face_offset->at(polyId + 1); ++verId){
				if (tag_data->at(verId) == true){
					uInt edgeEndVerId = TKCM::GetNextFaceVertexIndex(polyId, verId, face_offset);
					rest->at(edgeEndVerId) = true;
				}
			}
		}
		tag_data = std::move(rest);
		success = true;
	}

	void tag_mesh_border(
		const	uInt& point_count,
		const	Amino::Ptr<Amino::Array<uInt>>& face_vertex,
		const	Amino::Ptr<Amino::Array<uInt>>& face_offset,
		const	Amino::Ptr<Amino::Array<uInt>>& adjacent_edge_face,
		const	Amino::String& out_component_type,
				Amino::MutablePtr<Amino::Array<bool>>& tag_data,
				bool& success
	){
		tag_data = Amino::newMutablePtr<Amino::Array<bool>>();

		if (!face_vertex || face_vertex->empty() ||
			!face_offset || face_offset->empty() ||
			!adjacent_edge_face || adjacent_edge_face->empty())
		{
			success = false;
			return;
		}

		uInt vertexCount = int(face_vertex->size());
		uInt polyCount = int(face_offset->size()) - 1;

		std::string v("face_vertex_component");
		std::string p("point_component");
		std::string f("face_component");
		std::string e("half_edge_vertex_component");
		if (out_component_type == Amino::String(v.c_str(), v.size())){
			//////////////////////////////////////////////////////////// face_vertex
			tag_data->resize(vertexCount, false);

			#pragma omp parallel for
			for (size_t polyId = 0; polyId < polyCount; ++polyId){
				for (uInt verId = face_offset->at(polyId); verId < face_offset->at(polyId + 1); ++verId){
					uInt adjacentPolyId = adjacent_edge_face->at(verId);
					if (polyCount < adjacentPolyId){
						// �G�b�W�ŗאڂ���|���S���ԍ����K���l�ł͂Ȃ��i���אڂ����݂��Ȃ��j�ꍇ�̓^�O�t������
						tag_data->at(verId) = true;
						// �G�b�W�̃G���h�o�[�e�b�N�X���^�O�t������
						uInt edgeEndVerId = TKCM::GetNextFaceVertexIndex(polyId, verId, face_offset);
						tag_data->at(edgeEndVerId) = true;
					}
				}
			}
		} else if (out_component_type == Amino::String(p.c_str(), p.size())){
			//////////////////////////////////////////////////////////// point	
			tag_data->resize(point_count, false);

			#pragma omp parallel for 
			for (size_t polyId = 0; polyId < polyCount; ++polyId){
				for (uInt verId = face_offset->at(polyId); verId < face_offset->at(polyId + 1); ++verId){
					uInt adjacentPolyId = adjacent_edge_face->at(verId);
					if (polyCount < adjacentPolyId){
						// �G�b�W�ŗאڂ���|���S���ԍ����K���l�ł͂Ȃ��i���אڂ����݂��Ȃ��j�ꍇ�̓^�O�t������
						uInt poiId = face_vertex->at(verId);
						tag_data->at(poiId) = true;
					}
				}
			}
		} else if (out_component_type == Amino::String(f.c_str(), f.size())){
			//////////////////////////////////////////////////////////// face
			tag_data->resize(polyCount, false);

		#pragma omp parallel for 
			for (size_t polyId = 0; polyId < polyCount; ++polyId){
				for (uInt verId = face_offset->at(polyId); verId < face_offset->at(polyId + 1); ++verId){
					uInt adjacentPolyId = adjacent_edge_face->at(verId);
					if (polyCount < adjacentPolyId){
						// �G�b�W�ŗאڂ���|���S���ԍ����K���l�ł͂Ȃ��i���אڂ����݂��Ȃ��j�ꍇ���P�ł��������ꍇ�̓^�O�t������
						tag_data->at(polyId) = true;
						break;
					}
				}
			}
		}else if (out_component_type == Amino::String(e.c_str(), e.size())){
			//////////////////////////////////////////////////////////// half_edge_vertex
			tag_data->resize(vertexCount, false);

			#pragma omp parallel for
				for (size_t polyId = 0; polyId < polyCount; ++polyId){
					for (uInt verId = face_offset->at(polyId); verId < face_offset->at(polyId + 1); ++verId){
						uInt adjacentPolyId = adjacent_edge_face->at(verId);
						if (polyCount < adjacentPolyId){
							// �G�b�W�ŗאڂ���|���S���ԍ����K���l�ł͂Ȃ��i���אڂ����݂��Ȃ��j�ꍇ�̓^�O�t������
							tag_data->at(verId) = true;
						}
					}
				}
		} else{
			success = false;
			return;
		}
		success = true;
	}
}
