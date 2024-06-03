#ifndef BIF_TKCMCONVEXHULL_FNC_H
#define BIF_TKCMCONVEXHULL_FNC_H
#include <vector>
#include <algorithm>

// Bifrost Amino
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>
#include <Amino/Core/Array.h>

// TKCM
#include "../../core/TKCM_core.h"

using uInt = unsigned int;

namespace TKCM{
	bool isRightSide(const Bifrost::Math::float2& p, const Bifrost::Math::float2& p0, const Bifrost::Math::float2& p1){
		float z = ((p1.x - p0.x) * (p.y - p0.y)) - ((p.x - p0.x) * (p1.y - p0.y)); // cross product
		return 0.0f < z;
	}

	float distanceToLine(const Bifrost::Math::float2& p, const Bifrost::Math::float2& p0, const Bifrost::Math::float2& p1){
		return (abs(((p1.x - p0.x) * (p0.y - p.y)) - ((p0.x - p.x) * (p1.y - p0.y)))) / sqrt((pow((p1.x - p0.x), 2)) + (pow((p1.y - p0.y), 2)));
	}

	Amino::Array<uInt> upperHull(
		const int leftPoiID,
		const int rightPoiID,
		const Amino::Array<Bifrost::Math::float3>& uv,
		const int planeAxis,
		const Amino::Array<bool>& uvTag){
		Amino::Array<uInt> result;

		bool valid = false;
		for (int i = 0; i < uvTag.size(); ++i){
			if (uvTag[i]){
				valid = true;
				break;
			}
		}
		if (valid == false){ return result; }

		Amino::Array<bool> upper_uvs_tag(uv.size(), false);
		Bifrost::Math::float2 p, p0, p1;
		float maxDist = 0.0f;
		int furthestPoiID = -1;
		for (int i = 0; i < uv.size(); ++i){
			if (uvTag[i] == false){ continue; }
			switch (planeAxis){
				case 0: // x,z
					p = Bifrost::Math::float2{ uv[i].x, uv[i].z };
					p0 = Bifrost::Math::float2{ uv[leftPoiID].x, uv[leftPoiID].z };
					p1 = Bifrost::Math::float2{ uv[rightPoiID].x, uv[rightPoiID].z };
					break;
				case 1: // y,z
					p = Bifrost::Math::float2{ uv[i].y, uv[i].z };
					p0 = Bifrost::Math::float2{ uv[leftPoiID].y, uv[leftPoiID].z };
					p1 = Bifrost::Math::float2{ uv[rightPoiID].y, uv[rightPoiID].z };
					break;
				case 2: // x,y
					p = Bifrost::Math::float2{ uv[i].x, uv[i].y };
					p0 = Bifrost::Math::float2{ uv[leftPoiID].x, uv[leftPoiID].y };
					p1 = Bifrost::Math::float2{ uv[rightPoiID].x, uv[rightPoiID].y };
					break;
			}
			if (TKCM::isRightSide(p, p0, p1) == false){ continue; }

			upper_uvs_tag[i] = true;
			float dist = TKCM::distanceToLine(p, p0, p1);
			if (maxDist < dist){
				maxDist = dist;
				furthestPoiID = i;
			}
		}

		if (furthestPoiID != -1){
			result.push_back(furthestPoiID);
		}

		Amino::Array<uInt> region1 = TKCM::upperHull(leftPoiID, furthestPoiID, uv, planeAxis, upper_uvs_tag);
		Amino::Array<uInt> region3 = TKCM::upperHull(furthestPoiID, rightPoiID, uv, planeAxis, upper_uvs_tag);

		result.insert(result.end(), region1.begin(), region1.end());
		result.insert(result.end(), region3.begin(), region3.end());

		return result;
	}

	Amino::Array<uInt> GetSortedID(
		const Amino::Array<uInt>& uvID,
		const Amino::Array<Bifrost::Math::float3>& uv,
		const int& sortAxis,
		const bool& descendingOrder
	){
		if (uvID.size() < 2){ return uvID; }

		Amino::Array<uInt> result(uvID.size());
		std::vector<std::tuple<float, int>> val_id(uvID.size());
		#pragma omp parallel for if(1000 < cnt)
		for (int i = 0; i < uvID.size(); ++i){
			switch (sortAxis){
				case 0: std::get<0>(val_id[i]) = uv[uvID[i]].x; break; // x,z
				case 1: std::get<0>(val_id[i]) = uv[uvID[i]].y; break; // y,z
				case 2: std::get<0>(val_id[i]) = uv[uvID[i]].x; break; // x,y
			}
			std::get<1>(val_id[i]) = uvID[i];
		}
		std::sort(val_id.begin(), val_id.end());

		if (descendingOrder){
			#pragma omp parallel for if(1000 < cnt)
			for (int i = 0; i < uvID.size(); ++i){
				uInt ii = uvID.size() - i - 1;
				result[i] = std::get<1>(val_id[ii]);
			}
		} else{
			#pragma omp parallel for if(1000 < cnt)
			for (int i = 0; i < uvID.size(); ++i){
				result[i] = std::get<1>(val_id[i]);
			}
		}
		return result;
	}

	Amino::Array<uInt> convexHullPlane(const Amino::Array<Bifrost::Math::float3>& sourceUV, const int planeAxis){
		Amino::Array<uInt> result;

		if (sourceUV.empty()){ return result; }

		int leftPoiID = 0;
		int rightPoiID = 0;
		float minVal;
		switch (planeAxis){
			case 0: minVal = sourceUV.at(0).x; break; // x,z
			case 1: minVal = sourceUV.at(0).y; break; // y,z
			case 2: minVal = sourceUV.at(0).x; break; // x,y
			default: return result;
		}
		float maxVal = minVal;
		for (int i = 1; i < sourceUV.size(); ++i){
			float val;
			switch (planeAxis){
				case 0: val = sourceUV.at(i).x; break; // x,z
				case 1: val = sourceUV.at(i).y; break; // y,z
				case 2: val = sourceUV.at(i).x; break; // x,y
			}
			if (val < minVal){
				minVal = val;
				leftPoiID = i;
			}
			if (maxVal < val){
				maxVal = val;
				rightPoiID = i;
			}
		}
		if (leftPoiID == rightPoiID){ return result; }

		Amino::Array<bool> source_uvs_tag(sourceUV.size(), true);
		source_uvs_tag[leftPoiID] = false;
		source_uvs_tag[rightPoiID] = false;

		Amino::Array<uInt> allPointsUpperID = TKCM::upperHull(leftPoiID, rightPoiID, sourceUV, planeAxis, source_uvs_tag);
		Amino::Array<uInt> allPointsLowerID = TKCM::upperHull(rightPoiID, leftPoiID, sourceUV, planeAxis, source_uvs_tag);

		// ç∂âÒÇËÇ…êÆì⁄ÇµÇƒèoóÕÇ∑ÇÈ
		Amino::Array<uInt> dist_region1 = TKCM::GetSortedID(allPointsUpperID, sourceUV, planeAxis, false);
		Amino::Array<uInt> dist_region3 = TKCM::GetSortedID(allPointsLowerID, sourceUV, planeAxis, true);
		result.reserve(result.size() + allPointsUpperID.size() + allPointsLowerID.size() + 2);
		result.push_back(leftPoiID);
		result.insert(result.end(), dist_region1.begin(), dist_region1.end());
		result.push_back(rightPoiID);
		result.insert(result.end(), dist_region3.begin(), dist_region3.end());

		return result;
	}
}
#endif // BIF_TKCMCONVEXHULL_FNC_H
