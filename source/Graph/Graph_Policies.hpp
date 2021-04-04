#pragma once
#include "Graph.hpp"

///	@file
///	@brief	Contains the classes \ref mackey::bicolored_edge_with_id, \ref mackey::MinLength, \ref mackey::MinColorsLength

namespace mackey {

	///An edge with two possible colors and a numerical id
	struct bicolored_edge_with_id {
		char color;			
		unsigned char id;
		bicolored_edge_with_id(char color, unsigned char id);
	};

	///Prints edge to output stream in .dot format
	std::ostream& operator << (std::ostream& os, const bicolored_edge_with_id& edge);

	//Prints graph in .dot format with named nodes if provided
	template<typename neighborhood_t>
	std::ostream& operator<<(std::ostream& os, const Graph<neighborhood_t>& G);

	//Prints graph in .dot format with named nodes if provided
	template<typename graph_t, typename policy_t>
	std::ostream& operator<<(std::ostream& os, const ShortestPaths<graph_t, policy_t>& D);

	///ShortestPaths policy that minimizes the length of each path
	template<typename graph_t>
	class MinLength : public ShortestPaths<graph_t, MinLength<graph_t>> {
	public:
		MinLength(const graph_t& G); 	///<Computes shortest path in given graph
	private:
		typedef typename graph_t::node_t node_t;
		using ShortestPaths<graph_t, MinLength<graph_t>>::distance_policy;
		size_t path_length(node_t i) const;
		template<typename iter>
		bool check_and_update_policy(iter it, node_t current);
		void initialize() {};
		void clear() {};
		friend class ShortestPaths<graph_t, MinLength<graph_t>>;
	};

	///Dijkstra policy that minimizes the color alterations and length
	template<typename graph_t>
	class MinColorsLength : public ShortestPaths<graph_t, MinColorsLength<graph_t>> {
	public:
		MinColorsLength(const graph_t& G); ///<Computes shortest path in given graph
	private:
		typedef typename graph_t::node_t node_t;
		using ShortestPaths<graph_t, MinColorsLength<graph_t>>::distance_policy;
		std::vector<size_t> length;
		size_t path_length(node_t i) const;
		template<typename iter>
		bool check_and_update_policy(iter it, node_t current);
		void initialize();
		void clear();
		friend class ShortestPaths<graph_t, MinColorsLength<graph_t>>;
	};
}
#include "impl/Graph_Policies.ipp"
