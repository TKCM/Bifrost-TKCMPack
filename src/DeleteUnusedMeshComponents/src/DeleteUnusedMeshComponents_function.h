#pragma once
#include <vector>
#include <sstream>
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>

namespace TKCM
{
	int Array2DValueNextElement(const Amino::Array<Amino::Array<int>>& arrayVal, int arrayID, int val, Amino::Array<int>& condition){
		if (arrayVal.size() <= arrayID){ return -1; }

		int id;
		for (id = 0; id < arrayVal[arrayID].size(); ++id){
			if (val == arrayVal[arrayID][id]){ break; }
		}
		if (id == arrayVal[arrayID].size()){ return -1; }

		for (int i = 0; i < arrayVal[arrayID].size(); ++i){
			id++;
			if (id == arrayVal[arrayID].size()){ id = 0; }
			int nextVal = arrayVal[arrayID][id];
			if (condition[nextVal] < 0){
				continue;
			} else{
				return nextVal;
			}
		}
		return -1;
	}

	int Array2DValuePreviousElement(const Amino::Array<Amino::Array<int>>& arrayVal, int arrayID, int val, Amino::Array<int>& condition){
		if (arrayVal.size() <= arrayID){ return -1; }

		int id;
		for (id = 0; id < arrayVal[arrayID].size(); ++id){
			if (val == arrayVal[arrayID][id]){ break; }
		}
		if (id == arrayVal[arrayID].size()){ return -1; }

		for (int i = 0; i < arrayVal[arrayID].size(); ++i){
			id--;
			if (id == -1){ id = int(arrayVal[arrayID].size()) - 1; }
			int prevVal = arrayVal[arrayID][id];
			if (condition[prevVal] < 0){
				continue;
			} else{
				return prevVal;
			}
		}
		return -1;
	}
}