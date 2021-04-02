#pragma once
#include <vector>
#include <queue>
#include <set>
#include <iomanip> 

namespace mackey {

	/// @brief The adjecency list of a node in a graph
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

	///The shortest paths from a collection of sources to all points in a graph
	template<typename graph_t>
	class ShortestPaths {
	public:
		typedef typename graph_t::node_t node_t;								///<The node type
		typedef typename graph_t::neighborhood_t::const_iterator edge_iter_t;	///<The iterator type used in the paths
		///Comparator used for the Dijkstra algorithm
		bool operator()(node_t a, node_t b) const;
		///Print shortest path to node i to the output stream
		void print_path(std::ostream& os, node_t i) const;
		///	@brief	The shortest paths from the sources to all nodes 
		///	@details Stored as a sequence of iterators
		///	@warning	Do not invalidate the iterators by changing the graph!
		std::vector<std::vector<edge_iter_t>> paths;
		///The source of each path
		std::vector<node_t> root_per_path;
		///Returns 1 if given node is disconnected from the sources
		bool disconnected(node_t) const;
	protected:
		///Constructor computes the shortest paths given the graph
		ShortestPaths(const graph_t& G);
		///Constant reference to the graph
		const graph_t& G;
		///The distance from a source to a node, computed by some policy (eg weight)
		std::vector<size_t> distance_policy;
	private:
		void print_path_minimal(std::ostream& os, node_t i, std::set<node_t>& nodes_printed) const;
		void print_path_minimal(std::ostream& os, node_t i) const;
		template<typename _graph>
		friend std::ostream& operator<<(std::ostream&, const ShortestPaths<_graph>&);
	};

	///	@brief 			 The Dijkstra algorithm producing the shortest paths
	///	@tparam graph_t	 The type of the graph 
	///	@tparam policy_t The policy for the distance eg minimizing color changes
	///	@note			 \c policy_t is used via CRTP
	template<typename graph_t, typename policy_t>
	class Dijkstra : public ShortestPaths<graph_t>{
	public:
		typedef typename graph_t::node_t node_t;								///<The node type
		typedef typename ShortestPaths<graph_t>::edge_iter_t edge_iter_t;	///<The iterator type used in the paths
		///Compute all paths using given root
		void compute_with_root(node_t new_root);
		using ShortestPaths<graph_t>::paths;
		using ShortestPaths<graph_t>::root_per_path;
	protected:
		Dijkstra(const graph_t& G); ///<Saves G as const reference
		std::vector<std::pair<node_t, edge_iter_t>> previous; ///<previous[i]= node+edge that is 1 step closer to the source and node->i via edge
		node_t current; ///<The current node in the Dijkstra algorithm
		node_t root;	///<The root i.e. source
		using ShortestPaths<graph_t>::distance_policy;
		using ShortestPaths<graph_t>::G;
	private:
		void compute_distance();
		void compute_paths();
		void update(edge_iter_t it);
		std::priority_queue<node_t, std::vector<node_t>, const ShortestPaths<graph_t>&> next;
		std::vector<char> has_been_extracted;
	};

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
	template<typename graph_t>
	std::ostream& operator<<(std::ostream& os, const ShortestPaths<graph_t>& D);

	///Dijkstra policy that minimizes the length of each path
	template<typename graph_t>
	class MinLength : public Dijkstra<graph_t, MinLength<graph_t>> {
	public:
		MinLength(const graph_t& G); 	///<Computes shortest path in given graph
	private:
		typedef typename graph_t::node_t node_t;
		using Dijkstra<graph_t, MinLength<graph_t>>::current;
		using Dijkstra<graph_t, MinLength<graph_t>>::distance_policy;
		template<typename iter>
		bool check_and_update_policy(iter it);
		void initialize() {};
		void clear() {};
		size_t path_length(node_t i) const;	
		friend class Dijkstra<graph_t, MinLength<graph_t>>;
	};

	///Dijkstra policy that minimizes the color alterations and length
	template<typename graph_t>
	class MinColorsLength : public Dijkstra<graph_t, MinColorsLength<graph_t>> {
	public:
		MinColorsLength(const graph_t& G); ///<Computes shortest path in given graph
	private:
		typedef typename graph_t::node_t node_t;
		using Dijkstra<graph_t, MinColorsLength<graph_t>>::current;
		using ShortestPaths<graph_t>::distance_policy;
		std::vector<size_t> length;
		size_t path_length(node_t i) const;
		template<typename iter>
		bool check_and_update_policy(iter it);
		void initialize();
		void clear();
		friend class Dijkstra<graph_t, MinColorsLength<graph_t>>;
	};
}
#include "impl/Graph.ipp"
