#pragma once
#include "Utility/OpenMP_Macros.hpp"
#include "Utility/General.hpp"
#include "Spaces/Point.hpp"

///	@file
///	@brief Contains the class and methods to form the multiplication table.

namespace mackey
{
	/// @brief		The data of the Multiplication Table. 
	/// @details	Separated from the rest so as to be serializable. Has no functions
	template<typename group_t>
	struct MultTableData {
		int level; ///<The Mackey functor level we are working in		///	@brief	The nonzero homology groups in our given range and level and their identification information (if non cyclic). 
		std::map<std::vector<int>, IDGenerators<typename group_t::rank_t>> NonZeroHomology;
		std::vector<std::vector<Green<group_t>>> Greens;			///<Each entry in the multiplication table
		std::map<std::array<int, 3>, Green<group_t> > tripleGreens; ///<Triple products used when all other identification methods fail.
		std::vector<std::vector<int>> degree;		///< Maps each index to the corresponding degree
		std::map<std::vector<int>, int> antidegree; ///< Maps each degree to the corresponding index
		Eigen::Matrix<int, -1, -1> index_product;	///<Maps each index and irreducible to the index of their product
		std::vector<int> minsphere;///<The lower bound on the range of our spheres
		std::vector<int> maxsphere;///<The upper bound on the range of our spheres
		std::vector<std::vector<int>> basicIrreducibles;///<The basic irreducibles we use to produce the factorizations (eg Euler and orientation classes)
	};

	/// The Multiplication Table with both data and methods
	template<typename group_t>
	class MultTable : public MultTableData<group_t>{
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;

	public:

		using MultTableData<group_t>::level;
		using MultTableData<group_t>::NonZeroHomology;
		using MultTableData<group_t>::Greens;
		using MultTableData<group_t>::tripleGreens;
		using MultTableData<group_t>::degree;
		using MultTableData<group_t>::antidegree;
		using MultTableData<group_t>::index_product;
		using MultTableData<group_t>::minsphere;
		using MultTableData<group_t>::maxsphere;
		using MultTableData<group_t>::basicIrreducibles;

		/// @brief		Constructs multiplication table given the fundamental data (say after serialization)
		MultTable(const MultTableData<group_t>& MTD);

		/// @brief		Constructs multiplication table given the fundamental data (say after serialization)
		MultTable(MultTableData<group_t>&& MTD);

		///Constructs the multiplication table given the maximum and minimum spheres and the basic irreducibles.
		MultTable(int, const std::vector<int>&, const std::vector<int>&, const std::vector<std::vector<int>>&);

		/// Retrieve the degree index of the given degree. Returns -1 if no such degree can be found
		int getdegreeindex(const std::vector<int>&) const;

		/// @brief Checks if given degree is within the computed range
		/// @return 1 if it's in range and 0 otherwise
		bool degreewithinrange(const std::vector<int>&) const;

		Green<group_t> triple_product(int i, int j1, int j2) const;

	protected:

		/// Isolate sphere from the degree vector (i.e. get everything but first entry)
		std::vector<int> getsphere(const std::vector<int>& deg) const;


	private:
		Hash<std::vector<int>, std::vector<int>> perfect_hasher;
		std::vector<Chains<rank_t, diff_t>> basicChains;///<The Chains of the basic irreducibles.
		std::map<int, Chains<rank_t, diff_t>> sphereChains;
		void getIrreducibleChains();
		void multiply_all_indices();
		void multiply(int, int);
		void make();
		void generators();
		void compute(const std::vector<std::pair<int, int>>&);
	};
}
#include "impl/Mult_Table.ipp"