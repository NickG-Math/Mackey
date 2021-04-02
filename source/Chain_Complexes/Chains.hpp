#pragma once
#include <vector>
#include <ostream>

///@file
///@brief Contains Chains and Junction classes for using chain complexes.
namespace mackey
{

	///A differential with domain and range
	template <typename _rank, typename _diff>
	struct Arrow
	{
		typedef _rank rank_t; ///<The rank type
		typedef _diff diff_t; ///<The differential type
		rank_t domain; ///<Rank of the domain of the differential
		rank_t range; ///<Rank of the range of the differential
		diff_t diff; ///<The differential
		Arrow(const rank_t &, const rank_t &, const diff_t &);
		Arrow() = default;
		template <typename, typename>
		friend class Tensor;
	};

	/////////////////////////////////
	///	@brief A chain complex.
	///	@details Our chain complexes are always indexed as C[maxindex]--->C[0].
	/////////////////////////////////
	template <typename _rank, typename _diff>
	class Chains
	{
	public:
		typedef _rank rank_t; ///<The rank type
		typedef _diff diff_t; ///<The differential type
		std::vector<rank_t> rank; ///< rank[i]=rank(C[i])
		std::vector<diff_t> diff; ///< diff[i] is the differential C[i]->C[i-1] (diff[0] is always empty)

		Chains() = default; ///<Default constructor

		///The maximum index of the Chain complex
		int maxindex() const;
		Chains(const std::vector<rank_t> &, const std::vector<diff_t> &); ///<Constructor given ranks and diffs

		void reserve(size_t);

		void push_back(const Arrow<rank_t, diff_t> &);

		void push_back(const rank_t &, const diff_t &);

		///The dual cochain complex reindexed as a Chain complex. Optionally stop at index k (if k=-1 then nonstop i.e. k=maxindex)
		Chains dualize(int = -1) const;

		///Checks equality of Chain complexes
		bool operator==(const Chains<rank_t, diff_t> &) const;
	};

	///Prints chain complex
	template <typename rank_t, typename diff_t>
	std::ostream &operator<<(std::ostream &, const Chains<rank_t, diff_t> &);

	///////////////////////////////////////////////////////////////////////////////////////////////////
	///	@brief	Consisting of an entering and an exiting differential.
	/// @details	This is the part C[i+1]-->C[i]-->C[i-1] of a Chain complex C.\n
	/// C[i]-->C[i-1] is the "Out" part while C[i+1]-->C[i] is the "In" part.\n
	///	It's also the input of Homology.
	///////////////////////////////////////////////////////////////////////////////////////////////////
	template <typename _rank, typename _diff>
	class Junction
	{
	public:
		typedef _rank rank_t; ///<The rank type
		typedef _diff diff_t; ///<The differential type

		rank_t rank;	///<The rank of the middle group
		rank_t rankOut; ///<The rank of the rightmost group
		rank_t rankIn;	///<The rank of the leftmost group
		diff_t diffOut; ///<The exiting differential
		diff_t diffIn;	///<The entering differential


		/// Construct Junction by directly setting the elements
		Junction(const rank_t &, const rank_t &, const rank_t &, const diff_t &, const diff_t &);

		///Default constructor
		Junction() = default;

		/// Extract the Junction C[i+1]->C[i]->C[i-1] from the Chains C.
		Junction(const Chains<rank_t, diff_t>&, int);

	private:

		Junction(const Arrow<rank_t, diff_t>&, const Arrow<rank_t, diff_t>&);

		void setIn(const Arrow<rank_t, diff_t>&, bool = 1);

		void setIn(const rank_t&, const diff_t&);

		void setOut(const Arrow<rank_t, diff_t>&, bool = 1);

		void setOut(const rank_t&, const diff_t&);

		template <typename, typename>
		friend class Tensor;
	};
}
#include "impl/Chains.ipp"