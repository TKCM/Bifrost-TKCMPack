#include "array_test_nodedef.h"

namespace TKCM
{
	// 前提：
	
	// 　const有り引数は入力ポート、const無し引数は出力ポートになります

	// 　配列やジオメトリ、カスタムクラスのデータはAmino::Ptr(Amino::MutablePtr)タイプで宣言する必要があります
	// 　これは「大きなサイズになる可能性のあるデータはポインターを受け渡す」という仕様によるものです


	//////////////////////////////////////////////////////////////////////////////////////
	void array_input(
		const	int& get_index,
		const	Amino::Ptr<Amino::Array<int>>& input_array_data,
				int& out_data
	){
		// 配列データから特定の配列要素を出力する例：
		
		// Amino::PtrはC++標準のスマートポインタと同じ扱い方です
		// 唯一違う点がイミュータブルであるということです
		// 参照は出来ますが編集はできません

		if (!input_array_data || input_array_data->size() <= get_index){
			out_data = 1234;
		} else{
			out_data = input_array_data->at(get_index);
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////
	void array_output(
		const	int& value,
		const	int& out_size,
				Amino::MutablePtr<Amino::Array<int>>& out_data
	){
		//　配列を新規作成して出力する例

		// 配列のメモリを確保して、そのポインタを出力します
		// ミュータブルなので編集が可能です

		out_data = Amino::newMutablePtr<Amino::Array<int>>();
		out_data->resize(out_size, value);
	}

	//////////////////////////////////////////////////////////////////////////////////////
	void array_in_out(
				Amino::Ptr<Amino::Array<int>>& io_data,
		const	int& value,
		const	int& out_size
	){
		// 入力した配列を編集して出力する例
		
		// AMINO_ANNOTATEでInOutを宣言した引数はI/Oポートになります
		// 入力の受け取りはAmino::Ptr（=イミュータブル）です
		// Amino::MutablePtr（=ミュータブル）状態にすることで編集が可能になります

		// ミュータブルなポインタ変数に変換して内容を編集する
		Amino::MutablePtr<Amino::Array<int>> io_data_MPtr = io_data.toMutable();
		io_data_MPtr->resize(out_size, value);

		// ポインタ変数を出力ポートに割り当てる
		io_data = std::move(io_data_MPtr);
	}
}
