#pragma once
#include <vector>
#include <map>
#include <algorithm>
#include <array>

///@file
///@brief Contains the class \ref mackey::ChangeBasis
namespace mackey {

	/// @brief				If vector B is permutation of vector A i.e. B= A * sigma then returns permutation sigma
	/// @tparam pScalar		The value type of the result permutation eg size_t
	/// @tparam T			The value type of the input vectors eg size_t
	/// @param	A			The first input vector
	/// @param	B			The second input vector
	/// @return				Permutation sigma s.t. B[i]=A[sigma[i]]
	template<typename pScalar, typename T>
	std::vector<pScalar> find_permutation(const std::vector<T>& A, const std::vector<T>& B);


	/// @brief				Change of basis matrices needed for the equivariant tensor product
	/// @tparam pScalar		The scalar used for the permutation vectors eg size_t
	template<typename pScalar>
	struct ChangeBasis {
	public:
		std::vector<pScalar> left_to_canon; 	///<The permutation from the left convenient basis to the canonical basis
		std::vector<pScalar> right_to_canon;	///<The permutation from the right convenient basis to the canonical basis
		std::vector<pScalar> canon_to_left;		///<The permutation from the canonical basis to the left convenient basis
		std::vector<pScalar> canon_to_right;	///<The permutation from the canonical basis to the right convenient basis

		/// @brief						Constructs the change of basis matrices given two ranks.
		/// @tparam rank_t				The rank type
		/// @param  A					The first rank
		/// @param  B					The second rank
		/// @param set_conv_to_canon	Computes left_to_canon and right_to_canon if set to 1, and doesn't compute them if set to 0
		/// @param set_canon_to_conv	Computes canon_to_left and canon_to_right if set to 1, and doesn't compute them if set to 0
		template<typename rank_t>
		ChangeBasis(const rank_t& A, const rank_t& B, bool set_conv_to_canon, bool set_canon_to_conv); 
	};
}
#include "impl/ChangeBasis.ipp"
