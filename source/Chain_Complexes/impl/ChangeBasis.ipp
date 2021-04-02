#pragma once
#include "../ChangeBasis.hpp"

///@file
///@brief Contains the construction of the change of basis matrices.


namespace mackey {

	//If vector B is permutation of vector A i.e. B= A *sigma then returns permutation sigma
	template<typename pScalar, typename T>
	std::vector<pScalar> find_permutation(const std::vector<T>& A, const std::vector<T>& B)
	{
		std::conditional_t<std::is_integral_v<T>, std::vector<pScalar>, std::map<T, pScalar>> A_inverse; //The inverse of A
		if constexpr (std::is_integral_v<T>) {
			auto max = *std::max_element(A.begin(), A.end());
			A_inverse.resize(max + 1);
		}
		for (pScalar i = 0; i < A.size(); i++)
			A_inverse[A[i]] = i; //constructing the inverse
		std::vector<pScalar> sigma;
		sigma.reserve(A.size());
		for (const auto& i : B)
			sigma.push_back(A_inverse[i]); //sigma= A_inverse B
		return sigma;
	}

	template<typename pScalar>
	template<typename rank_t>
	ChangeBasis<pScalar>::ChangeBasis(const rank_t& rankC, const rank_t& rankD, bool set_conv_to_canon, bool set_canon_to_conv) {
		//Constructs the 3 bases (left+right convenient and canonical) and returns the two change of basis matrices (as permutations).

		size_t s = rankC.size();
		size_t t = rankD.size();
		if (s == 0 || t == 0) {
			std::cerr << "Can't use empty ranks with this function!";
			abort();
		}

		//Perfect hashing since we have the maximum entries of the array (minimums are all 0)
		Hash<std::array<size_t, 3>>H({ s - 1,t - 1, std::max((size_t)rankC.maxCoeff(), (size_t)rankD.maxCoeff()) - 1 });
		auto totalsize = H.max_hash();

		std::vector<size_t> canonical, leftconv, rightconv;
		canonical.reserve(totalsize);
		leftconv.reserve(totalsize);
		rightconv.reserve(totalsize);

		//first write the canonical in hashed form
		for (size_t d = 0; d < t; d++)
			for (size_t c = 0; c < s; c++)
				for (scalar_t<rank_t> j = 0; j < std::min(rankC[c], rankD[d]); j++)
					for (scalar_t<rank_t> i = 0; i < std::max(rankC[c], rankD[d]); i++)
						canonical.push_back(H(std::array<size_t, 4>({ c, d, (size_t)(i % rankC[c]), (size_t)((i + j) % rankD[d]) })));

		//next write the left convenient in hashed form
		for (size_t d = 0; d < t; d++)
			for (scalar_t<rank_t> j = 0; j < rankD[d];j++)
				for (size_t c = 0; c < s; c++)
					for (scalar_t<rank_t> i = 0; i < rankC[c]; i++)
						leftconv.push_back(H(std::array<size_t, 4>({ c, d, (size_t)i, (size_t)j })));

		//finally write the right convenient in hashed form
		for (size_t c = 0; c < s; c++)
			for (scalar_t<rank_t> i = 0; i < rankC[c]; i++)
				for (size_t d = 0; d < t; d++)
					for (scalar_t<rank_t> j = 0; j < rankD[d];j++)
						rightconv.push_back(H(std::array<size_t, 4>({ c, d , (size_t)i, (size_t)j })));

		//get permutation matrices by comparing the hashed forms
		if (set_conv_to_canon) {
			left_to_canon = find_permutation<pScalar>(canonical, leftconv); // conv[i]=can[conv_to_canon[i]]
			right_to_canon = find_permutation<pScalar>(canonical, rightconv);
		}
		if (set_canon_to_conv) {
			canon_to_left = find_permutation<pScalar>(leftconv, canonical); // can[i]=con[can_to_conv[i]]
			canon_to_right = find_permutation<pScalar>(rightconv, canonical);
		}
	}
}
