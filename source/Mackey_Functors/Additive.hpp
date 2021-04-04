#pragma once
#include "Utility/OpenMP_Macros.hpp"
#include "Spaces/Point.hpp"
#include <iomanip>
///@file
///@brief Contains the class \ref mackey::AdditiveStructure 

namespace mackey {

	///The additive structure of the RO(G) homology of a point
	template<typename group_t>
	class AdditiveStructure {
	public:
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;

		std::vector<int> minsphere;///<The lower bound on the range of our spheres
		std::vector<int> maxsphere;///<The upper bound on the range of our spheres

		///Compute and identify given the range of our spheres
		AdditiveStructure(const std::vector<int>& minsphere, const std::vector<int>& maxsphere);

		///Retrieve the Mackey functors in the homology of the given sphere
		template<typename sphere_t>
		std::vector<MackeyFunctor<rank_t>> getMackey(const sphere_t& sphere) const;

		///Retrieve the Mackey functor for the given degree
		template<typename sphere_t>
		MackeyFunctor<rank_t> getMackey(int degree, const sphere_t& sphere) const;

		///Print the unique Mackey functors that appear
		template<typename T>
		void print_unique(T& stream);

		///Print the unnamed Mackey functors that appear
		template<typename T>
		void print_unknown(T& stream);

		///The identified Mackey functors (one from each 
		const std::vector<MackeyFunctor<rank_t>>& identified();

		///The identified Mackey functors (one from each 
		const std::vector<MackeyFunctor<rank_t>>& unknown();

	private:

		///Identify the Mackey functors
		void identify();

		std::vector<int> true_unknowns;
		std::vector<MackeyFunctor<rank_t>> uniqueMackeys;
		std::vector<MackeyFunctor<rank_t>> unknownMackeys;
		std::vector<std::vector<MackeyFunctor<rank_t>>> allMackeys;
		std::vector<char> found_unknown;
		void compute();
		void identify_unique(MackeyFunctor<rank_t>&, bool);
		void identify_unknown(int);

		void check_unique_and_assign(const MackeyFunctor<rank_t>&);
		void check_unknown_and_assign(const std::vector<MackeyFunctor<rank_t>>&);

		template<typename _group>
		friend std::ostream& operator<<(std::ostream& os, const AdditiveStructure<_group>&);
	};

	/// @brief				Prints all Mackey functors in their universal notation to the output stream
	/// @tparam group_t		The group type of the AdditiveStructure
	template<typename group_t>
	std::ostream& operator<<(std::ostream& os, const AdditiveStructure<group_t>& A);
}

#include "impl/Additive.ipp"
