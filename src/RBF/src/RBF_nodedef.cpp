#include "RBF_nodedef.h"

void TKCM::setup(
	const Amino::Ptr<Amino::Array<Amino::Ptr<Amino::Array<float>>>>& key,
	const Amino::Ptr<Amino::Array<Amino::Ptr<Amino::Array<float>>>>& value,
	const int type_,
	Amino::MutablePtr<TKCM::RBFSolver>& rbf,
	bool& valid
){
	rbf = Amino::newMutablePtr<TKCM::RBFSolver>();

	if (key->size() == 0 || value->size() == 0 || key->size() != value->size()){ 
		rbf->type = -1;
		valid = false; 
		return; 
	}

	int keyVCnt = key->at(0)->size();
	int valueVCnt = value->at(0)->size();
	for (int i = 1; i < key->size(); ++i){
		if (keyVCnt != key->at(0)->size()){ return; }
		if (valueVCnt != value->at(0)->size()){ return; }
	}

	Eigen::MatrixXf targetMat = TKCM::RBFFn::array2DToMat(key);
	Eigen::MatrixXf retargetMat = TKCM::RBFFn::array2DToMat(value);

	// RBF setup
	valid = TKCM::RBFFn::setup(*rbf, targetMat, retargetMat, type_);
}

void TKCM::solve(
	const Amino::Ptr<TKCM::RBFSolver>& rbf,
	const Amino::Ptr<Amino::Array<float>>& driver,
	Amino::MutablePtr<Amino::Array<float>>& result,
	bool& success
){
	result = Amino::newMutablePtr<Amino::Array<float>>();
	if (rbf == nullptr){ result->resize(0); return; }

	uInt cnt = uInt(driver->size());
	Eigen::VectorXf driverPos(cnt);
	driverPos = Eigen::VectorXf::Map(&driver->at(0), cnt);

	Eigen::VectorXf result_;
	success = RBFFn::solve(result_, *rbf, driverPos);
	
	if (success){
		for (int i = 0; i < result_.size(); ++i){
			result->push_back(result_[i]);
		}
	}
}
