#pragma once
#include "../Graph.hpp"

namespace mackey {
	template<typename node_t, typename edge_t, template<typename ...> typename container_t, typename ...Args>
	template<typename ... Argshere>
	void Neighborhood<node_t, edge_t, container_t, Args...>::insert(node_t end, Argshere&&... args) {
		for (const auto i:edges)
			if (i.first == end) //never insert duplicate edges
				return;
		edges.emplace_back(end, edge_t(args...));
	}

	template<typename node_t, typename edge_t, template<typename ...> typename container_t, typename ...Args>
	void Neighborhood<node_t, edge_t, container_t, Args...>::reserve(size_t n) {
		edges.reserve(n);
	}

	template<typename node_t, typename edge_t, template<typename ...> typename container_t, typename ...Args>
	auto Neighborhood<node_t, edge_t, container_t, Args...>::const_iterator::operator ++() -> const_iterator& {
		++it;
		return *this;
	}
	template<typename node_t, typename edge_t, template<typename ...> typename container_t, typename ...Args>
	auto Neighborhood<node_t, edge_t, container_t, Args...>::const_iterator::node() const -> node_t {
		return it->first;
	}
	template<typename node_t, typename edge_t, template<typename ...> typename container_t, typename ...Args>
	auto Neighborhood<node_t, edge_t, container_t, Args...>::const_iterator::edge() const->const edge_t& {
		return it->second;
	}
	template<typename node_t, typename edge_t, template<typename ...> typename container_t, typename ...Args>
	bool Neighborhood<node_t, edge_t, container_t, Args...>::const_iterator::operator!=(const_iterator other) const {
		return it != other.it;
	}

	template<typename node_t, typename edge_t, template<typename ...> typename container_t, typename ...Args>
	Neighborhood<node_t, edge_t, container_t, Args...>::const_iterator::const_iterator(iter_t it) : it(it) {}

	template<typename node_t, typename edge_t, template<typename ...> typename container_t, typename ...Args>
	auto Neighborhood<node_t, edge_t, container_t, Args...>::begin() const -> const_iterator {
		return const_iterator(edges.begin());
	}

	template<typename node_t, typename edge_t, template<typename ...> typename container_t, typename ...Args>
	auto Neighborhood<node_t, edge_t, container_t, Args...>::end() const -> const_iterator {
		return const_iterator(edges.end());
	}


	template<typename neighborhood_t>
	template<typename ... Args>
	void Graph< neighborhood_t>::draw_edge(node_t start, node_t end, Args&&... args) {
		neighborhoods[start].insert(end, args...);
	}

	template<typename neighborhood_t>
	template<typename T>
	void Graph< neighborhood_t>::reserve_edges_per_node(const T& v) {
		for (size_t i = 0; i < neighborhoods.size(); i++) {
			if constexpr (std::is_integral_v<T>)
				neighborhoods[i].reserve(v);
			else
				neighborhoods[i].reserve(v[i]);
		}
	}

	template<typename neighborhood_t>
	void Graph< neighborhood_t>::reserve_edges_for_node(size_t num, node_t node) {
		neighborhoods[node].reserve(num);
	}

	template<typename neighborhood_t>
	size_t Graph<neighborhood_t>::number_of_nodes() const {
		return neighborhoods.size();
	}
	template<typename neighborhood_t>
	void Graph<neighborhood_t>::reserve_nodes(size_t n) {
		neighborhoods.reserve(n);
	}
	template<typename neighborhood_t>
	void Graph<neighborhood_t>::adjoin_nodes(size_t new_number_of_nodes) {
		neighborhoods.resize(new_number_of_nodes);
	}
	template<typename neighborhood_t>
	auto Graph<neighborhood_t>::neighbors(node_t node) const -> const neighborhood_t& {
		return neighborhoods[node];
	}

	template<typename neighborhood_t>
	void Graph< neighborhood_t>::print_name_labels(std::ostream& os) const {
		if (node_names != nullptr)
		{
			for (size_t i = 0; i < number_of_nodes(); i++)
				os << i << "[label=\"" << node_names(i) << "\"] \n";
		}
	}

	template<typename graph_t>
	bool ShortestPaths<graph_t>::operator()(node_t a, node_t b) const {
		return (distance_policy[a] > distance_policy[b]);
	}

	template<typename graph_t>
	void ShortestPaths<graph_t>::print_path(std::ostream& os, node_t i) const {
		os << "digraph{\ngraph [overlap=false];	\n";
		print_path_minimal(os, i);
		G.print_name_labels(os);
		os << "}";
	}

	template<typename graph_t>
	void ShortestPaths<graph_t>::print_path_minimal(std::ostream& os, node_t i, std::set<node_t>& nodes_printed) const {
		const auto& p = paths[i];
		for (size_t ind = 0; ind + 1 < p.size(); ind++) {
			if (nodes_printed.find(p[ind].node()) == nodes_printed.end()) {
				os << p[ind + 1].node() << " -> " << p[ind].node() << p[ind].edge() << "\n";
				nodes_printed.insert(p[ind].node());
			}
			else
				return;
		}
		if (!p.empty()) {
			if (nodes_printed.find(p.back().node()) == nodes_printed.end()) {
				os << root_per_path[i] << " -> " << p.back().node() << p.back().edge() << "\n";
				nodes_printed.insert(p.back().node());
			}
		}
	}

	template<typename graph_t>
	void ShortestPaths<graph_t>::print_path_minimal(std::ostream& os, node_t i) const {
		const auto& p = paths[i];
		for (size_t ind = 0; ind + 1 < p.size(); ind++)
			os << p[ind + 1].node() << " -> " << p[ind].node() << p[ind].edge() << "\n";
		if (!p.empty())
			os << root_per_path[i] << " -> " << p.back().node() << p.back().edge() << "\n";
	}

	template<typename graph_t>
	ShortestPaths<graph_t>::ShortestPaths(const graph_t& G) : G(G){}

	template<typename graph_t>
	bool ShortestPaths<graph_t>::disconnected(node_t i) const
	{
		return (root_per_path[i] == -1);
	}

	namespace implementation_details {
		template<typename T>
		std::vector<T> reserve(size_t n) {
			std::vector<T> cont;
			cont.reserve(n);
			return cont;
		};
	}

	template<typename graph_t, typename policy_t>
	Dijkstra<graph_t, policy_t>::Dijkstra(const graph_t& G) : ShortestPaths<graph_t>(G), next(*this, implementation_details::reserve<node_t>(G.number_of_nodes())) {}

	template<typename graph_t, typename policy_t>
	void Dijkstra<graph_t, policy_t>::compute_distance() {
		distance_policy = std::vector<size_t>(G.number_of_nodes(), -1);
		previous.resize(G.number_of_nodes());
		static_cast<policy_t*>(this)->initialize();
		current = root;
		has_been_extracted = std::vector<char>(G.number_of_nodes(), 0);
		distance_policy[current] = 0;
		next.push(current);
		while (!next.empty())
		{
			current = next.top();
			next.pop();
			if (has_been_extracted[current])
				continue;
			for (auto it = G.neighbors(current).begin(); it != G.neighbors(current).end(); ++it)
				update(it);
			has_been_extracted[current] = 1;
		}
		has_been_extracted.clear();
	}

	template<typename graph_t, typename policy_t>
	void Dijkstra<graph_t, policy_t>::compute_paths() {
		paths.resize(G.number_of_nodes());
		root_per_path.resize(G.number_of_nodes(),-1);
		for (size_t i = 0; i < G.number_of_nodes(); i++) {
			if (distance_policy[i] == -1)
				continue;
			auto length = static_cast<const policy_t*>(this)->path_length(i);
			if (!paths[i].empty() && paths[i].size() <= length)
				continue;
			root_per_path[i] = root;
			paths[i].clear();
			paths[i].reserve(length);
			auto start = i;
			while (start != root)
			{
				paths[i].push_back(previous[start].second);
				start = previous[start].first;
			}
		}
		previous.clear();
		distance_policy.clear();
		static_cast<policy_t*>(this)->clear();
	}

	template<typename graph_t, typename policy_t>
	void Dijkstra<graph_t, policy_t>::update(edge_iter_t it) {
		if (static_cast<policy_t*>(this)->check_and_update_policy(it)) {
			previous[it.node()] = std::pair(current, it); //<Is this slow?
			next.push(it.node());
		}
	}

	template<typename graph_t, typename policy_t>
	 void Dijkstra<graph_t, policy_t>::compute_with_root(node_t new_root){
		root = new_root;
		compute_distance();
		compute_paths();
	}

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

	template<typename neighborhood_t>
	std::ostream& operator<<(std::ostream& os, const Graph<neighborhood_t>& G)
	{
		os << "digraph{\ngraph [overlap=false];	\n";
		for (size_t i = 0; i < G.number_of_nodes(); i++)
			for (auto it = G.neighbors(i).begin(); it != G.neighbors(i).end(); ++it)
				os << i << "-> " << it.node() << " " << it.edge() << "\n";
		G.print_name_labels(os);
		os << "}";
		return os;
	}

	template<typename graph_t>
	std::ostream& operator<<(std::ostream& os, const ShortestPaths<graph_t>& D)
	{
		os << "digraph{\ngraph [overlap=false];	\n";
		std::set<typename graph_t::node_t> nodes_printed;
		for (size_t i = 0; i < D.G.number_of_nodes(); i++)
			D.print_path_minimal(os, i, nodes_printed);
		D.G.print_name_labels(os);
		os << "}";
		return os;
	}


	template<typename graph_t>
	MinLength<graph_t>::MinLength(const graph_t& G) : Dijkstra<graph_t, MinLength<graph_t>>(G) {}

	template<typename graph_t>
	template<typename iter>
	bool MinLength<graph_t>::check_and_update_policy(iter it) {
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
	MinColorsLength<graph_t>::MinColorsLength(const graph_t& G) : Dijkstra<graph_t, MinColorsLength<graph_t>>(G) {
		initialize();
	}

	template<typename graph_t>
	template<typename iter>
	bool MinColorsLength<graph_t>::check_and_update_policy(iter it) {
		auto newdist_policy = distance_policy[current] + 1 + 1*(current != this->root && this->previous[current].second.edge().color%2 != it.edge().color%2); //add 1 for extra length, 3 for color alteration
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