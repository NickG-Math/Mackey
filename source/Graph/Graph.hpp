#pragma once
#include <vector>
#include <queue>
#include <set>
#include <iomanip> 

///	@file
///	@brief	Contains the classes \ref mackey::Neighborhood, \ref mackey::Graph, \ref mackey::ShortestPaths

namespace mackey {

	/// 	@brief			The adjecency list of a node in a graph
	///	@tparam _node 		The node type should be an integer type
	///	@tparam	_edge 		Any edge type
	///	@tparam	container	Currently only vector is supported
	///	@tparam Args		Extra arguments ofr the container (eg allocator)
	template<typename _node, typename _edge, template<typename ...> typename container_t, typename ...Args>
	class Neighborhood {
	public:
		typedef _node node_t; ///<The node type
		typedef _edge edge_t; ///<The edge type
		///	Insert node+edge in neighborhood
		template<typename ... Argshere>
		void insert(node_t end, Argshere&&... args);
		/// Reserve number of edges
		void reserve(size_t n);
		/// Constant iterator traversing edges with fixed start
		struct const_iterator {
			auto operator ++()->const_iterator&;
			auto node() const->node_t;	///<The end node of the edge
			auto edge() const->const edge_t&;
			bool operator!=(const_iterator other) const;
			const_iterator() = default;
		private:
			typedef typename container_t<std::pair<node_t, edge_t>, Args...>::const_iterator iter_t;
			const_iterator(iter_t it);
			iter_t it;
			friend class Neighborhood;
		};
		auto begin() const->const_iterator; ///<Starting iterator
		auto end() const->const_iterator;	///<Ending iterator
	private:
		container_t<std::pair<node_t, edge_t>, Args...> edges;
	};

	///	@brief				  A directed graph
	///	@tparam	_neighborhood The Neighborhood type
	template<typename _neighborhood>
	class Graph {
	public:
		typedef _neighborhood neighborhood_t; 			///<The neighborhood type
		typedef typename neighborhood_t::edge_t edge_t;	///<The edge type
		typedef typename neighborhood_t::node_t node_t;	///<The node type
		Graph() = default;
		size_t number_of_nodes() const;

		///Connects start->end. The arguments are forwarded to the edge constructor
		template<typename ... Args>
		void draw_edge(node_t start, node_t end, Args&&... args);

		///Reserves num many edges for given node
		void reserve_edges_for_node(size_t num, node_t node);

		template<typename T>
		///Reserves v[i] many edges for node i
		void reserve_edges_per_node(const T& v);

		///...Reserves n many edges for all nodes
		void reserve_nodes(size_t n);

		///Adds new nodes so that the total is new_number_of_nodes
		void adjoin_nodes(size_t new_number_of_nodes);

		///The adjecency list of a node
		auto neighbors(node_t node) const -> const neighborhood_t&;

		///Function mapping node index to its name
		std::function<std::string(node_t)> node_names;

		///Prints names as .dot style labels to output stream
		void print_name_labels(std::ostream& os) const;
	private:
		std::vector<neighborhood_t> neighborhoods;
	};
	
	/// @brief				The shortest paths from a collection of sources to all points in a graph
	/// @tparam graph_t		The type of directed graph
	/// @tparam policy_t	The policy that determines what we are minimizing against (eg total weight or number of color changes)
	/// @details			Implementation of the Dijkstra algorithm (without decrease key), using the policy to determine if we should update the distance\n
	///						Specify the policy via CRTP: have \c policy_t inherit from this class
	///						Policy examples include \ref MinLength and \ref MinColorsLength				
	template<typename graph_t, typename policy_t>
	class ShortestPaths {
	public:
		typedef typename graph_t::node_t node_t;								///<The node type
		typedef typename graph_t::neighborhood_t::const_iterator edge_iter_t;	///<The iterator type used in the paths

		///Print shortest path to node i to the output stream
		void print_path(std::ostream& os, node_t i) const;

		///	@brief		The shortest paths from the sources to all nodes 
		///	@details	Stored as a sequence of iterators
		///	@warning	Do not invalidate the iterators by changing the graph!
		std::vector<std::vector<edge_iter_t>> paths;

		///The source of each path
		std::vector<node_t> root_per_path;

		///Compute all paths using given root
		void compute_with_root(node_t new_root);

		///	@brief	Checks if node is connected to any one of the sources
		///	@return	Returns \c true if given node is disconnected and \c false if connected
		bool disconnected(node_t) const;

	protected:
		///Constructor computes the shortest paths given the graph
		ShortestPaths(const graph_t& G);

		///Constant reference to the graph
		const graph_t& G;

		///The distance from a source to a node, computed by some policy (eg weight)
		std::vector<size_t> distance_policy;

		///previous[i]=(node,edge) where node is 1 step closer to the source and edge connects node to i
		std::vector<std::pair<node_t, edge_iter_t>> previous;

		///The root (source) of the graph being considered
		node_t root;

	private:
		void print_path_minimal(std::ostream& os, node_t i, std::set<node_t>& nodes_printed) const;
		void print_path_minimal(std::ostream& os, node_t i) const;
		void compute_distance();
		void compute_paths();
		template<typename _graph, typename _pol>
		friend std::ostream& operator<<(std::ostream&, const ShortestPaths<_graph, _pol>&);
	};

	//Prints graph in .dot format with named nodes if provided
	template<typename neighborhood_t>
	std::ostream& operator<<(std::ostream& os, const Graph<neighborhood_t>& G);

	//Prints graph in .dot format with named nodes if provided
	template<typename graph_t, typename policy_t>
	std::ostream& operator<<(std::ostream& os, const ShortestPaths<graph_t, policy_t>& D);

}
#include "impl/Graph.ipp"
