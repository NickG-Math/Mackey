#pragma once
#include "Mult_Graph.hpp"

///@file
///@brief Contains the class \ref mackey::MultIdentify

namespace mackey {

	/// Provides extra identification methods using triple box products
	template<typename group_t>
	class MultIdentify {
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;

	public:
		/// Uses triple box products to possibly identify ALL instances where identification failed before. Use with care
		void pass_all_unidentified();

		///	Stores multiplication graph as a reference
		MultIdentify(MultGraph<group_t>& MG);

	private:
		MultGraph<group_t>& MG;
		std::vector<std::array<int, 3>> triples_to_be_done;
		std::set<std::array<int, 3>> identified;
		std::pair<int, std::vector<rank_t>> distinguish(int, const std::vector<rank_t>&);
		std::vector<rank_t> distinguish(int, const std::vector<rank_t>&, int);
		void identify_connection(int, int);
		void make_connection(int, int, int);
	};
}
#include "impl/Mult_Identify.ipp"
