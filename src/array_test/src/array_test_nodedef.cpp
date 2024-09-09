#include "array_test_nodedef.h"

namespace TKCM
{
	// �O��F
	
	// �@const�L������͓��̓|�[�g�Aconst���������͏o�̓|�[�g�ɂȂ�܂�

	// �@�z���W�I���g���A�J�X�^���N���X�̃f�[�^��Amino::Ptr(Amino::MutablePtr)�^�C�v�Ő錾����K�v������܂�
	// �@����́u�傫�ȃT�C�Y�ɂȂ�\���̂���f�[�^�̓|�C���^�[���󂯓n���v�Ƃ����d�l�ɂ����̂ł�


	//////////////////////////////////////////////////////////////////////////////////////
	void array_input(
		const	int& get_index,
		const	Amino::Ptr<Amino::Array<int>>& input_array_data,
				int& out_data
	){
		// �z��f�[�^�������̔z��v�f���o�͂����F
		
		// Amino::Ptr��C++�W���̃X�}�[�g�|�C���^�Ɠ����������ł�
		// �B��Ⴄ�_���C�~���[�^�u���ł���Ƃ������Ƃł�
		// �Q�Ƃ͏o���܂����ҏW�͂ł��܂���

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
		//�@�z���V�K�쐬���ďo�͂����

		// �z��̃��������m�ۂ��āA���̃|�C���^���o�͂��܂�
		// �~���[�^�u���Ȃ̂ŕҏW���\�ł�

		out_data = Amino::newMutablePtr<Amino::Array<int>>();
		out_data->resize(out_size, value);
	}

	//////////////////////////////////////////////////////////////////////////////////////
	void array_in_out(
				Amino::Ptr<Amino::Array<int>>& io_data,
		const	int& value,
		const	int& out_size
	){
		// ���͂����z���ҏW���ďo�͂����
		
		// AMINO_ANNOTATE��InOut��錾����������I/O�|�[�g�ɂȂ�܂�
		// ���͂̎󂯎���Amino::Ptr�i=�C�~���[�^�u���j�ł�
		// Amino::MutablePtr�i=�~���[�^�u���j��Ԃɂ��邱�ƂŕҏW���\�ɂȂ�܂�

		// �~���[�^�u���ȃ|�C���^�ϐ��ɕϊ����ē��e��ҏW����
		Amino::MutablePtr<Amino::Array<int>> io_data_MPtr = io_data.toMutable();
		io_data_MPtr->resize(out_size, value);

		// �|�C���^�ϐ����o�̓|�[�g�Ɋ��蓖�Ă�
		io_data = std::move(io_data_MPtr);
	}
}
