#ifndef BIF_ARRAYTEST_H
#define BIF_ARRAYTEST_H

#include <vector>
#include <algorithm>

// Bifrost Amino
#include <Amino/Cpp/Annotate.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/Array.h>

#include "../../core/TKCM_export.h"

namespace TKCM
{
	/////////////////////////////////////////////////////////////////////
	TKCM_DECL
	void array_input(
		const	int& get_index,
		const	Amino::Ptr<Amino::Array<int>>& input_array_data,
				int& out_data
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=array_test::input "
	);

	TKCM_DECL
	void array_output(
		const	int& value,
		const	int& out_size,
				Amino::MutablePtr<Amino::Array<int>>& out_data
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=array_test::output "
	);

	TKCM_DECL
	void array_in_out(
				Amino::Ptr<Amino::Array<int>>& io_data AMINO_ANNOTATE("Amino::InOut inName=io_data_in outName=io_data_out"),
		const	int& value,
		const	int& out_size
	) AMINO_ANNOTATE(
		"Amino::Node "
		"name=array_test::io "
	);
}

#endif // BIF_ARRAYTEST_H
