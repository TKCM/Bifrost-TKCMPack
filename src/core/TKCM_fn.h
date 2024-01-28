#pragma once
#include <vector>
#include <sstream>
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>

namespace TKCM {
	template<typename T>
	T ArrayIndexNextElement(const Amino::Array<T>& arrayVal, const int index ){
		if(index == arrayVal.size()-1){
			return arrayVal.at(0);
		}else{
			return arrayVal.at(index+1);
		}
	}

	template<typename T>
	T ArrayValueNextElement(const Amino::Array<T>& arrayVal, const T thisElement){
		int id;
		for (id = 0; id < arrayVal.size(); ++id){
			if (thisElement == arrayVal[id]){ break; }
		}

		if (id == arrayVal.size()){
			return -1;
		} else if (id == arrayVal.size() - 1){
			return arrayVal.at(0);
		} else{
			return arrayVal.at(id + 1);
		}
	}
	
	template<typename T>
	T ArrayValuePreviousElement(const Amino::Array<T>& arrayVal, const T thisElement){
		int id;
		for (id = 0; id < arrayVal.size (); ++id) {
			if (thisElement == arrayVal[id]) { break; }
		}

		if (id == arrayVal.size ()) {
			return -1;
		}else if (id == 0) {
			return arrayVal.back();
		}else {
			return arrayVal[id - 1];
		}
	}
		
	template<typename T>
	bool AlmostEqual(const T A, const T B, const T tolerance = T(0.0001)){
		return std::abs(A - B) < tolerance;
	}
	
	template<typename T>
	T Lerp ( const T a, const T b, const float t ){
		return a + ((b - a) * t);
	}

}