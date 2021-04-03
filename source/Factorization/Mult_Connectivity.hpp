#pragma once
#include "Mult_Graph.hpp"

///@file
///@brief Contains methods for the finding the connectivity of the multiplication graph should identification fail.

namespace mackey {

	///	@brief	Finds the connectivity of the multiplication graph
	///	@details Provides methods for filling the multiplication graph with "wrong" edges that don't affect the connectivity of the graph
	template<typename group_t>
	class MultConnectivity : private MultGraph<group_t> {
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;

	public:
		template<typename...Args>
		MultConnectivity(Args&&...);

		///The generators disconnected from the sources (element indices)
		std::vector<int> disconnected_indices;

		///The degrees of the disconnected_indices
		std::vector<std::vector<int>> disconnected_degrees;

		///Computes the nodes of the Multiplication Graph disconnected from the given sources
		void compute_with_sources(const std::vector<std::vector<int>>&);

	private:
		MinLength<typename MultGraph<group_t>::graph_t> shortest_paths;
		std::set<int> before_disconnected;
		std::vector<int> sources;
		void set_sources(const std::vector<std::vector<int>>&);
		void doinstages(bool);
		void wrong_connection(int, int, bool);
	};
}
#include "impl/Mult_Connectivity.ipp"