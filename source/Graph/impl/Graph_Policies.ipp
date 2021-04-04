#pragma once
#include "../Graph_Policies.hpp"

namespace mackey {
	bicolored_edge_with_id::bicolored_edge_with_id(char color, unsigned char id) : color(color), id(id) {}

	std::ostream& operator << (std::ostream& os, const bicolored_edge_with_id& edge) {
		os << "[color=\"";
		if (edge.color == 0)
			os << "red"; //red for multiplication
		else if (edge.color == 1)
			os << "blue"; //blue for division
		else
			os << "black"; //black for multiple
		os << "\"";
		os << ", tooltip=\"" << std::to_string(edge.id) << "\"";
		os << "]";
		return os;
	}

	template<typename graph_t>
	MinLength<graph_t>::MinLength(const graph_t& G) : ShortestPaths<graph_t, MinLength<graph_t>>(G) {}

	template<typename graph_t>
	template<typename iter>
	bool MinLength<graph_t>::check_and_update_policy(iter it, node_t current) {
		auto newdist_policy = distance_policy[current] + 1;
		if (distance_policy[it.node()] == -1 || distance_policy[it.node()] > newdist_policy) {
			distance_policy[it.node()] = newdist_policy;
			return 1;
		}
		return 0;
	}

	template<typename graph_t>
	size_t MinLength<graph_t>::path_length(node_t i) const {
		return distance_policy[i];
	}

	template<typename graph_t>
	MinColorsLength<graph_t>::MinColorsLength(const graph_t& G) : ShortestPaths<graph_t, MinColorsLength<graph_t>>(G) {
		initialize();
	}

	template<typename graph_t>
	template<typename iter>
	bool MinColorsLength<graph_t>::check_and_update_policy(iter it, node_t current) {
		auto newdist_policy = distance_policy[current] + 1 + 1 * (current != this->root && this->previous[current].second.edge().color % 2 != it.edge().color % 2); //add 1 for extra length, 3 for color alteration
		if (distance_policy[it.node()] == -1 || distance_policy[it.node()] > newdist_policy) {
			distance_policy[it.node()] = newdist_policy;
			length[it.node()] = length[current] + 1;
			return 1;
		}
		return 0;
	}

	template<typename graph_t>
	size_t MinColorsLength<graph_t>::path_length(node_t i) const {
		return length[i];
	}

	template<typename graph_t>
	void MinColorsLength<graph_t>::initialize() {
		length = std::vector<size_t>(this->G.number_of_nodes(), 0);
	}

	template<typename graph_t>
	void MinColorsLength<graph_t>::clear() {
		length.clear();
	}

}
