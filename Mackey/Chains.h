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

		rank_t rankOut, rankIn;
		diff_t diffOut, diffIn;
		
		///Default constructor
		Junction() {};
	/// Extract the Junction C[i+1]->C[i]->C[i-1] from the Chains C.
		Junction(const Chains<rank_t, diff_t>&, int);
	/// Construct Junction by directly setting the elements
		Junction(const rank_t& rank, const rank_t& rankOut, const rank_t& rankIn, const diff_t& diffOut, const diff_t& diffIn) : rank(rank), rankOut(rankOut), rankIn(rankIn), diffOut(diffOut), diffIn(diffIn) {}
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
