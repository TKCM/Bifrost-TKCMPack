#pragma once
#include "../../common/TKCM_core.h"

namespace TKCM{
	namespace Delete {
		namespace Mesh {
			void EraseUnusedComponentTag ( Bifrost::Object& polygonMesh );

			void CreateNormalProp (
				Bifrost::Object& target,
				const Bifrost::Object& source,
				const Amino::Array<int>& poiCondition,
				const Amino::Array<bool>& verCondition,
				const Amino::Array<bool>& faceCondition );

			void CreateUvProp (
				Bifrost::Object& target,
				const Bifrost::Object& source,
				const Amino::Array<bool>& verCondition);

			void CreatePointComponentProp ( 
				Bifrost::Object& target,
				const Bifrost::Object& source,
				const Amino::Array<int>& poiCondition );

			void CreateVertexComponentProp ( 
				Bifrost::Object& target,
				const Bifrost::Object& source,
				const Amino::Array<bool>& verCondition );

			void CreateFaceComponentProp ( 
				Bifrost::Object& target,
				const Bifrost::Object& source,
				const Amino::Array<bool>& faceCondition );

			//////////////////////////////////////////////////////////////////
			template<class T>
			void CreatePointCompPropCore ( 
				Bifrost::Object& target,
				const Bifrost::Object& source,
				const Amino::String& propName, 
				const Amino::Array<int>& poiCondition );

			template<class T>
			void CreateVertexCompPropCore (
				Bifrost::Object& target,
				const Bifrost::Object& source,
				const Amino::String& propName,
				const Amino::Array<bool>& verCondition );

			template<class T>
			void CreateFaceCompPropCore ( 
				Bifrost::Object& target,
				const Bifrost::Object& source,
				const Amino::String& propName, 
				const Amino::Array<bool>& faceCondition );
		}
	}

	namespace Dissolve {
		namespace Mesh {
			void CreateNormalProp (
				Bifrost::Object& target,
				const Bifrost::Object& source,
				const Amino::Array<int>& poiCondition,
				const Amino::Array<UInt32>& verOrder,
				const Amino::Array<bool>& faceCondition );

			void CreateUvProp (
				Bifrost::Object& target,
				const Bifrost::Object& source,
				const Amino::Array<UInt32>& verOrder );

			void CreateVertexComponentProp (
				Bifrost::Object& target,
				const Bifrost::Object& source,
				const Amino::Array<UInt32>& verOrder );

			//////////////////////////////////////////////////////////////////
			template<class T>
			void CreateVertexCompPropCore (
				Bifrost::Object& target,
				const Bifrost::Object& source,
				const Amino::String& propName,
				const Amino::Array<UInt32>& verOrder );
		}
	}
}

