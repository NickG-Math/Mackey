#pragma once
#include "Mult_Graph.hpp"

///@file
///@brief Contains extra identification methods for the multiplication graph.

namespace mackey {

	/// Provides extra identification methods using triple box products
	template<typename group_t>
	class MultIdentify {
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;

	public:
		/// Uses triple box products to possibly identify the products of the given generators if identification failed before on them
		void pass_product(const std::vector<int>&);

		/// Uses triple box products to possibly identify the products landing in degrees of the given generators if identification failed before on them
		void pass_division(const std::vector<int>&);

		/// Uses triple box products to possibly identify ALL instances where identification failed before. Use with care
		void pass_all_unidentified();

		bool can_do_more; ///<1 if there are more triple products that can be computed for the disconnected generators

		///Uses the Multiplication Graph constructor
		MultIdentify(MultGraph<group_t>& MG);

	private:
		MultGraph<group_t>& MG;
		void pass_identified();
		std::vector<std::array<int, 3>> triples_to_be_done;
		std::vector<std::pair<int, int>> identified;
		void delete_identified();
		/// Uses triple box products to possibly identify given instance where identification failed before.
		bool pass_triple(std::pair<int, int>, bool);
		std::pair<int, std::vector<rank_t>> distinguish(int, const std::vector<rank_t>&);
		rank_t determine_connection_identify(const Green<group_t>&, int, int, int, bool);
	};
}
#include "impl/Mult_Identify.ipp"
