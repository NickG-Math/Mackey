#pragma once
#include <vector>

///@file
///@brief Contains Chains and Junction classes for using chain complexes.

namespace Mackey {

	/////////////////////////////////
	/// A Chain complex.

	/// Our chain complexes are always indexed as C[maxindex]--->C[0].
	/////////////////////////////////
	template<typename rank_t, typename diff_t>
	class Chains {
	public:
		std::vector<rank_t> rank; ///< rank[i]=rank(C[i])
		std::vector<diff_t> diff; ///< diff[i] is the differential C[i]->C[i-1] (diff[0] is always empty)
		int maxindex; ///<The maximum index of the Chain complex

		Chains() {}///<Default constructor

		Chains(const std::vector<rank_t>& rank, const std::vector<diff_t>& diff) : rank(rank), diff(diff) 
		{ maxindex = rank.size() - 1; }	///<Constructor given ranks and diffs

		///The dual of a Chain complex reindexed as a Chain complex and up to index k
		Chains dualize(int k) {
			std::vector<rank_t> rank_dual;
			std::vector<diff_t> diff_dual;
			rank_dual.reserve(maxindex + 1);
			diff_dual.reserve(maxindex + 1);
			for (int i = 0; i <= k; i++) {
				auto index = maxindex - i;
				rank_dual.push_back(rank[index]);
				if (i == 0) {
					diff_t empty;
					diff_dual.push_back(empty);
				}
				else
					diff_dual.push_back(diff[index + 1].transpose());
			}
			return Chains(rank_dual, diff_dual);
		}

		///The dual of a Chain complex reindexed as a Chain complex (as opposed to cochains)
		Chains dualize() {
			return dualize(maxindex);
		}

	};

	/////////////////////////////////
	/// The input for Homology, consisting of an entering and an exiting differential.

	/// This is the part C[i+1]-->C[i]-->C[i-1] of a Chain complex C.
	/// C[i]-->C[i-1] is the "Out" part while C[i+1]-->C[i] is the "In" part.
	/////////////////////////////////
	template<typename rank_t, typename diff_t>
	class Junction
	{
	public:

		rank_t rank;///<The rank of the middle group
		rank_t rankOut;///<The rank of the rightmost group
		rank_t rankIn;///<The rank of the leftmost group
		diff_t diffOut;///<The exiting differential
		diff_t diffIn;///<The entering differential


	///Default constructor
		Junction() {};
	/// Extract the Junction C[i+1]->C[i]->C[i-1] from the Chains C.
		Junction(const Chains<rank_t, diff_t>&, int);
	/// Construct Junction by directly setting the elements
		Junction(const rank_t& rank, const rank_t& rankOut, const rank_t& rankIn, const diff_t& diffOut, const diff_t& diffIn) 
		: rank(rank), rankOut(rankOut), rankIn(rankIn), diffOut(diffOut), diffIn(diffIn) {}
	};



	template<typename rank_t, typename diff_t>
	Junction<rank_t, diff_t>::Junction(const Chains<rank_t, diff_t>& C, int i)
	{
		rank = C.rank[i];
		if (0 < i) {
			rankOut = C.rank[i - 1];
			diffOut = C.diff[i];
		}
		if (i < C.maxindex) {
			rankIn = C.rank[i + 1];
			diffIn = C.diff[i + 1];
		}
	}
}
