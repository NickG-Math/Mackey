#pragma once
#include "../Graph.hpp"

namespace mackey {
	template<typename node_t, typename edge_t, template<typename ...> typename container_t, typename ...Args>
	template<typename ... Argshere>
	void Neighborhood<node_t, edge_t, container_t, Args...>::insert(node_t end, Argshere&&... args) {
		for (const auto i : edges)
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

	template<typename graph_t, typename policy_t>
	void ShortestPaths<graph_t, policy_t>::print_path(std::ostream& os, node_t i) const {
		os << "digraph{\ngraph [overlap=false];	\n";
		print_path_minimal(os, i);
		G.print_name_labels(os);
		os << "}";
	}

	template<typename graph_t, typename policy_t>
	void ShortestPaths<graph_t, policy_t>::print_path_minimal(std::ostream& os, node_t i, std::set<node_t>& nodes_printed) const {
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

	template<typename graph_t, typename policy_t>
	void ShortestPaths<graph_t, policy_t>::print_path_minimal(std::ostream& os, node_t i) const {
		const auto& p = paths[i];
		for (size_t ind = 0; ind + 1 < p.size(); ind++)
			os << p[ind + 1].node() << " -> " << p[ind].node() << p[ind].edge() << "\n";
		if (!p.empty())
			os << root_per_path[i] << " -> " << p.back().node() << p.back().edge() << "\n";
	}

	template<typename graph_t, typename policy_t>
	ShortestPaths<graph_t, policy_t>::ShortestPaths(const graph_t& G) : G(G) {}

	template<typename node_t>
	struct Node_Distance_Pair {
		node_t node;
		size_t distance;
		bool operator < (const Node_Distance_Pair& other) const {
			return distance > other.distance;
		}
	};


	template<typename graph_t, typename policy_t>
	bool ShortestPaths<graph_t, policy_t>::disconnected(node_t i) const
	{
		return (root_per_path[i] == -1);
	}

	template<typename graph_t, typename policy_t>
	void ShortestPaths<graph_t, policy_t>::compute_distance() {
		distance_policy = std::vector<size_t>(G.number_of_nodes(), -1); //initialize with infinity
		previous.resize(G.number_of_nodes());
		static_cast<policy_t*>(this)->initialize(); //initialize policy
		distance_policy[root] = 0;

		std::vector<Node_Distance_Pair<node_t>> cont; //lets us reserve in the priority queue
		cont.reserve(G.number_of_nodes());
		std::priority_queue<Node_Distance_Pair<node_t>> next(std::less<Node_Distance_Pair<node_t>>(), cont);

		auto current = root;
		next.push({ current,0 });
		while (!next.empty())
		{
			auto ndp = next.top();
			next.pop();
			current = ndp.node;
			if (ndp.distance <= distance_policy[current]) //since we are not decreasing key, this lets us know that we are not using an outdated node
				for (auto it = G.neighbors(current).begin(); it != G.neighbors(current).end(); ++it) {
					if (static_cast<policy_t*>(this)->check_and_update_policy(it, current)) { //check and update the policy in the inheriting class
						previous[it.node()] = std::pair(current, it);
						next.push({ it.node(),distance_policy[it.node()] }); //if it.node() existed before the new it.node() will have higher priority
					}
				}
		}
	}

	template<typename graph_t, typename policy_t>
	void ShortestPaths<graph_t, policy_t>::compute_paths() {
		paths.resize(G.number_of_nodes());
		root_per_path.resize(G.number_of_nodes(), -1);
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
	void ShortestPaths<graph_t, policy_t>::compute_with_root(node_t new_root) {
		root = new_root;
		compute_distance();
		compute_paths();
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

	template<typename graph_t, typename policy_t>
	std::ostream& operator<<(std::ostream& os, const ShortestPaths<graph_t, policy_t>& D){
		os << "digraph{\ngraph [overlap=false];	\n";
		std::set<typename graph_t::node_t> nodes_printed;
		for (size_t i = 0; i < D.G.number_of_nodes(); i++)
			D.print_path_minimal(os, i, nodes_printed);
		D.G.print_name_labels(os);
		os << "}";
		return os;
	}
}
