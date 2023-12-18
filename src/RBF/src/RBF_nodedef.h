#ifndef BIF_TKCMRBF_H
#define BIF_TKCMRBF_H
#include <vector>
#include <algorithm>

// Eigen
#include <Eigen/Dense>

// Bifrost Amino
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Any.h>

// TKCM
#include "../../core/TKCM_core.h"

using uInt = unsigned int;

namespace TKCM
{
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	class AMINO_ANNOTATE("Amino::Class") TKCM_DECL RBFSolver{
	public:
		RBFSolver(){};
		~RBFSolver(){};

		Eigen::MatrixXf P, V, W;
		int type;
	};

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	TKCM_DECL
	void setup(
		const Amino::Ptr<Amino::Array<Amino::Ptr<Amino::Array<float>>>>& key,
		const Amino::Ptr<Amino::Array<Amino::Ptr<Amino::Array<float>>>>& value,
		const int type,
		Amino::MutablePtr<TKCM::RBFSolver>& rbf,
		bool& valid
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=RBF::Internal::setup "
		"metadata=[{icon, ../icon/tkcm_internal.png}, {internal, string, true}] "
	);

	TKCM_DECL
	void solve(
		const Amino::Ptr<TKCM::RBFSolver>& rbf,
		const Amino::Ptr<Amino::Array<float>>& driver,
		Amino::MutablePtr<Amino::Array<float>>& result,
		bool& success
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=RBF::solve "
		"metadata=[{icon, ../icon/tkcm.png}, {internal, string, false}] "
	);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	namespace RBFFn
	{
		bool setup(TKCM::RBFSolver& rbf, const Eigen::MatrixXf& P, const Eigen::MatrixXf& F, const int type);
		bool solve(Eigen::VectorXf& result, const TKCM::RBFSolver& rbf, const Eigen::VectorXf& driverPose);

		float Norm(int type, float radius);

		Eigen::MatrixXf array2DToMat(const Amino::Ptr<Amino::Array<Amino::Ptr<Amino::Array<float>>>>& data);
	}
} // namespace TKCM_BIF
#endif // BIF_TKCMRBF_H
