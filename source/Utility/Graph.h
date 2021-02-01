#pragma once
#include<vector>
#include<array>
#include<fstream>
#include<queue>
#include "Types/SFINAE.h"

///@file
///@brief Contains weighted/colored graph classes and Dikjstra algorithms.

namespace Mackey {


	//////////////////////////////////////
	/// A directed Graph

	/// This is a general class that can only be used by inheriting from it.
	/// T is the type of the graph eg weightedGraph or coloredGraph. It's used for the CRTP compile-time polymorphism. By default it's void
	//////////////////////////////////////
	template<typename T = void>
	class Graph {
	public:
		int number_of_nodes; ///<The number of vertices

		std::vector<std::vector<int>> edges; ///<The edges of the graph

		//////////////////////////////////////
		/// The shortest path (in terms of total length of edges) from a selected source to any other point. 

		/// The length of each edge depends on the specilization of Graph

		///path[i] consists of the procession of nodes from i to source.
		//////////////////////////////////////
		std::vector<std::vector<int>> path;

		/// Same as path, but using the procession of edges as opposed to nodes. Not currenltly used/set
		std::vector<std::vector<int>> path_edges;

		/////////////////////////////////////
		/// The distance from a selected source to all points. -1 if disconnected

		/// The length of each edge depends on the specilization of Graph
		//////////////////////////////////////
		std::vector<int> distance;

		/// The disconnected points from each source
		std::vector<int> disconnected;

		///Default constructor
		Graph() {};

		/// Construct given the edges of the graph
		Graph(const std::vector<std::vector<int>>& edges) : number_of_nodes(edges.size()), edges(edges) {}

		/// Construct from a graph with different template
		template<typename S>
		Graph(const Graph<S>& G) : Graph(G.edges) {}


		/// Compute all paths using the given source
		void computeWithSource(int givensource);

		///Compute the disconnected nodes up to a certain node
		void compute_disconnected(int limit);

		///Compute all disconnected nodes
		void compute_disconnected() {
			compute_disconnected(number_of_nodes);
		}

		///Default implementation for graphs uses length=1 for all edges
		void computePath_default();

		///Clears distances and initializes
		void clear();

	protected:
		int source; ///<A source of the graph
	};



	template<typename T>
	void Graph<T>::computeWithSource(int givensource) {
		if (path.size() < number_of_nodes) { //either hasn't been run, or new nodes were added
			path.resize(number_of_nodes);
			distance.assign(number_of_nodes, -1);
		}
		source = givensource;
		if constexpr (SFINAE::has_computePath<T>::value) //if it has an implementation of computePath
			static_cast<T*>(this)->computePath();
		else //otherwise use default implementation with length=1 for all edges
			computePath_default();
	}

	template<typename T>
	void Graph<T>::compute_disconnected(int limit) {
		disconnected.reserve(limit);
		disconnected.clear();
		for (int i = 0; i < limit; i++) {
			if (distance[i] == -1)
				disconnected.push_back(i);
		}
	}

	template<typename T>
	void Graph<T>::clear() {
		path.clear();
		path.resize(number_of_nodes);
		distance.assign(number_of_nodes, -1);
		disconnected.clear();
		if constexpr (SFINAE::has_initialize<T>::value)
			static_cast<T*>(this)->initialize(); //initialize
	}

	/// A directed Graph with weights
	class WeightedGraph : public Graph<WeightedGraph> {
	public:

		///Default constructor
		WeightedGraph() {};

		///Constructor from given graph and weights
		template<typename T>
		WeightedGraph(const Graph<T>&, const std::vector<std::vector<int>>&);

		///Constructor from given edges and weights
		WeightedGraph(const std::vector<std::vector<int>>& edges, const std::vector<std::vector<int>>& weights);

		///Constructor from given graph and zero weights
		template<typename T>
		WeightedGraph(Graph<T>& G);


	private:
		std::vector<std::vector<int>> weights;
		std::vector<int> unweighted_distance, closest;
		std::vector<char> visited, reachable; //not bool for performance
		std::priority_queue<std::pair<int, int>> next;
		std::vector<std::vector<int>> allWeightsAreOne(const std::vector<std::vector<int>>&);
		void initialize();
		void computePath();
		void stepDistance(int);
		void stepPath(int, int);
		void computeDistance();
		friend class Graph<WeightedGraph>; ///<Used to set up the CRTP.
		friend class SFINAE;
	};

	template<typename T>
	WeightedGraph::WeightedGraph(const Graph<T>& G, const std::vector<std::vector<int>>& weights) : Graph<WeightedGraph>(G), weights(weights) { initialize(); }

	WeightedGraph::WeightedGraph(const std::vector<std::vector<int>>& edges, const std::vector<std::vector<int>>& weights)
		: Graph<WeightedGraph>(edges), weights(weights) {
		initialize();
	}

	template<typename T>
	WeightedGraph::WeightedGraph(Graph<T>& G) : WeightedGraph(G, allWeightsAreOne(G.edges)) { initialize(); }


	std::vector<std::vector<int>> WeightedGraph::allWeightsAreOne(const std::vector<std::vector<int>>& e) {
		std::vector<std::vector<int>> w;
		w.resize(e.size());
		for (size_t i = 0; i < e.size(); i++) {
			w[i].assign(e[i].size(),1);
		}
		return w;
	}

	void WeightedGraph::initialize() {
		visited.assign(number_of_nodes, 0);
		unweighted_distance.assign(number_of_nodes, 0);
		closest = unweighted_distance;
	}


	void WeightedGraph::stepDistance(int start) {
		for (size_t j = 0; j < edges[start].size(); j++) {
			auto i = edges[start][j];
			if (!visited[i]) {
				if (distance[i] == -1 || distance[i] > distance[start] + weights[start][j]) {
					distance[i] = distance[start] + weights[start][j];
					unweighted_distance[i] = unweighted_distance[start] + 1;
					closest[i] = start;
					next.push(std::make_pair(-distance[i], i));
				}
				else if (distance[i] == distance[start] + weights[start][j] && unweighted_distance[i] > unweighted_distance[start]+1){
					unweighted_distance[i] = unweighted_distance[start] + 1;
					closest[i] = start;
					next.push(std::make_pair(-distance[i], i));
				}
			}
		}
		visited[start] = 1;
		reachable[start] = 1;
	}

	void WeightedGraph::computeDistance() {
		reachable.assign(number_of_nodes, 0);
		closest[source] = source;
		distance[source] = 0;
		int start = source;
		next.push(std::make_pair(-distance[source], source));
		while (!next.empty()) {
			start = next.top().second;
			next.pop();
			stepDistance(start);
		}
	}

	void WeightedGraph::stepPath(int i, int j) {
		if (closest[j] != source) {
			path[i].push_back(closest[j]);
			stepPath(i, closest[j]);
		}
	}

	void WeightedGraph::computePath() {
		computeDistance();
		for (int i = 0; i < number_of_nodes;i++) {
			if (reachable[i]) {
				path[i].reserve(unweighted_distance[i] + 1);
				path[i].push_back(i);
				stepPath(i, i);
				path[i].push_back(source);
			}
		}
	}



	template<typename T>
	void Graph<T>::computePath_default() {
		WeightedGraph W(*this);
		W.computeWithSource(source);
		path = std::move(W.path);
		distance = std::move(W.distance);
		path_edges = std::move(W.path_edges);
	}


	/// A directed Graph with two colors
	class ColoredGraph : public Graph<ColoredGraph> {
	public:
		std::vector<std::vector<char>> colors;	///<The colors of the graph

		/// Default Constructor
		ColoredGraph() {};

		/// Construct given Graph and colors
		template<typename T>
		ColoredGraph(const Graph<T>& G, std::vector<std::vector<char>>& colors) : Graph<ColoredGraph>(G), colors(colors) {}

	private:
		void computePath();
		void constructDual();
		std::vector<int> adjustpath(std::vector<int>& path);
		WeightedGraph dual;

		friend class Graph<ColoredGraph>; ///<Used to set up the CRTP.
		friend class SFINAE;
	};


	std::vector<int> ColoredGraph::adjustpath(std::vector<int>& path) {
		for (auto& i : path) { i = i / 2; }
		return path;
	}

	void ColoredGraph::constructDual() {
		//we create a new Graph with colored nodes and monochrome edges
		//every node i of the original graph gives two colored nodes i_red, i_blue and 
		//the monochrome edges ending at i_red are only the red edges of the original graph (same for blue)
		//the monochrome edges that alternate node colors are given weight 2

		std::vector<std::vector<int>> dual_edges;
		std::vector<std::vector<int>> weights;
		dual_edges.resize(2 * number_of_nodes);
		weights.resize(2 * number_of_nodes);
		for (int i = 0; i < number_of_nodes; i++) {

			dual_edges[2 * i].reserve(edges[i].size());
			weights[2 * i].reserve(edges[i].size());
			dual_edges[2 * i + 1].reserve(edges[i].size());
			weights[2 * i + 1].reserve(edges[i].size());

			for (std::vector<int>::size_type j = 0; j < edges[i].size(); j++) {
				if (!colors[i][j]) {//red target
					dual_edges[2 * i].push_back(2 * edges[i][j]); //red to red
					weights[2 * i].push_back(1);

					dual_edges[2 * i + 1].push_back(2 * edges[i][j]); //blue to red
					weights[2 * i + 1].push_back(2);
				}
				else { //blue target
					dual_edges[2 * i].push_back(2 * edges[i][j] + 1); //red to blue
					weights[2 * i].push_back(2);

					dual_edges[2 * i + 1].push_back(2 * edges[i][j] + 1); //blue to blue
					weights[2 * i + 1].push_back(1);
				}
			}
		}
		dual = WeightedGraph(dual_edges, weights);
	}

	void ColoredGraph::computePath() {
		constructDual();
		dual.computeWithSource(2 * source);
		auto redpath = std::move(dual.path);
		auto reddistance = std::move(dual.distance);
		dual.clear();
		dual.computeWithSource(2 * source + 1);
		auto bluepath = std::move(dual.path);
		auto bluedistance = std::move(dual.distance);

		for (int i = 0; i < number_of_nodes; i++) {
			int mindistance = -1;
			std::array<int, 4> V = { reddistance[2 * i], reddistance[2 * i + 1], bluedistance[2 * i], bluedistance[2 * i + 1] };
			int finder = -1;
			for (int j = 0; j < 3; j++) {
				if (V[j] != -1 && (mindistance == -1 || V[j] < mindistance)) {
					mindistance = V[j];
					finder = j;
				}
			}
			if (finder == 0 || finder == 1) {//use red
				if (distance[i] == -1 || (distance[i] > reddistance[2 * i + finder] || (distance[i] == reddistance[2 * i + finder] && path[i].size() > redpath[2 * i + finder].size()))) {
					path[i] = adjustpath(redpath[2 * i + finder]);
					distance[i] = reddistance[2 * i + finder];
				}
			}
			else if (finder == 2 || finder == 3) {//use blue
				if (distance[i] == -1 || (distance[i] > bluedistance[2 * i + finder - 2] || (distance[i] == bluedistance[2 * i + finder - 2] && path[i].size() > bluepath[2 * i + finder - 2].size()))) {
					path[i] = adjustpath(bluepath[2 * i + finder - 2]);
					distance[i] = bluedistance[2 * i + finder - 2];
				}
			}
		}
	}

	///////////////////////////
	///Class that allows the printing of graphs in graphviz format with/without provided names for the nodes

	///The template graph_t can be Graph<void> , ColoredGraph etc. 
	///The template node_names_t can be void (if nodes are not named), std::vector<std::string> or more generally any class with an operator std::string [](int) const
	template<typename graph_t, typename node_names_t = void>
	struct GraphPrinter {
		const graph_t& graph; ///<Reference to the graph to be printed
		const node_names_t* const node_names; ///<Pointer to the names of the nodes. This can be either nullptr with type void or not nullptr and type not void
		
		///Constructor with graph and optionally node names
		GraphPrinter(const graph_t& graph, const node_names_t* const node_names = nullptr)
			: graph(graph), node_names(node_names) {}
	};

	///Adjoins node names to graph to get it ready for printing
	template<typename graph_t, typename node_names_t>
	typename std::enable_if_t<SFINAE::is_graph_t<graph_t>::value, GraphPrinter<graph_t, node_names_t>> operator<<(const graph_t& G, const node_names_t& node_names) {
		return GraphPrinter<graph_t, node_names_t>(G, &node_names);
	}

	///Prints graph in .dot format with named nodes if provided
	template<typename graph_t, typename node_names_t >
	std::ostream& operator<<(std::ostream& os, const GraphPrinter<graph_t, node_names_t>& GP) {
		os << "digraph G{ \n";
		//if constexpr (std::is_same<node_names_t, void>::value)
		//	os << "node[shape=point]\n ";
		const std::array<std::string, 2> coloring = { "red","blue" };
		for (int i = 0; i < GP.graph.number_of_nodes; i++) {
			for (int j = 0; j < GP.graph.edges[i].size(); j++) {
					os << i << "-> " << GP.graph.edges[i][j];
				if constexpr (std::is_same<graph_t, ColoredGraph>::value)
					os << " [color=\"" << coloring[GP.graph.colors[i][j]] << "\"]";
				os << "\n";
			}
		}
		if constexpr (!std::is_same<node_names_t, void>::value) {
			for (int i = 0; i < GP.graph.number_of_nodes; i++) 
				os << i << "[label =\"" << (*GP.node_names)[i] << "\"] \n";
		}
		os << "}";
		return os;
	}

	/// Prints the graph in .dot format with unnamed nodes.
	template<typename graph_t>
	typename std::enable_if_t<SFINAE::is_graph_t<graph_t>::value, std::ostream&> operator<<(std::ostream& os, const graph_t& G) {
		os << GraphPrinter<graph_t>(G);
		return os;
	}
}
