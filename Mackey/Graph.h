#pragma once
#include<vector>
#include<array>
#include<fstream>
#include<queue>
#include "SFINAE.h"

///@file
///@brief Contains weighted/colored graph classes and Dikjstra algorithms.

namespace Mackey {


	//////////////////////////////////////
	/// A directed Graph

	/// This is a general class that can only be used by inheriting from it.
	/// T is the type of the graph eg weightedGraph or coloredGraph. It's used for the CRTP compile-time polymorphism. By default it's void
	//////////////////////////////////////
	template<typename T=void>
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

		/// Compute all paths using the given source
		void computeWithSource(int givensource);

		///Compute the disconnected nodes up to a certain node
		void compute_disconnected(int limit);

		///Compute all disconnected nodes
		void compute_disconnected() {
			compute_disconnected(number_of_nodes);
		}
		
		/// Writes a graph.dot file representing the graph. The nodes are unnamed.
		void draw();
		/// Writes a graph.dot file representing the graph using node names.
		void draw(const std::vector<std::string>&);


		///Default implementation for graphs uses length=1 for all edges
		void computePath_default();

		///Clears distances and initializes
		void clear();

	protected:
		int source; ///<A source of the graph
	};



	template<typename T>
	void Graph<T>::computeWithSource(int givensource){
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



	template<typename T>
	void Graph<T>::draw() {
		std::ofstream file;
		file.open("graph.dot");
		file << "digraph G{ \n node[shape=point] \n";
		for (int i = 0; i < number_of_nodes; i++) {
			for (auto j : edges[i]) {
				file << i << "->" << j << "\n";
			}
		}
		file << "}";
		file.close();
	}

	template<typename T>
	void Graph<T>::draw(const std::vector<std::string>& names) {
		std::ofstream file;
		file.open("graph.dot");
		file << "digraph G{ \n";
		for (int i = 0; i < number_of_nodes; i++) {
			for (auto j : edges[i]) {
				file << "\"" << names[i] << "\"" << "->" << "\"" << names[j] << "\" \n";
			}
		}
		file << "}";
		file.close();
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
		WeightedGraph(const std::vector<std::vector<int>> & edges, const std::vector<std::vector<int>>& weights);
		
		///Constructor from given graph and zero weights
		template<typename T>
		WeightedGraph(Graph<T>& G);
		
	
	private:
		std::vector<std::vector<int>> weights;
		std::vector<int> unweighted_distance, closest;
		std::vector<char> visited, reachable; //not bool for performance
		std::priority_queue<std::pair<int, int>> next;
		std::vector<std::vector<int>> zeroWeight(const std::vector<std::vector<int>>&);
		void initialize();
		void computePath();
		void stepDistance(int);
		void stepPath(int, int);
		void computeDistance();
		friend class Graph<WeightedGraph>; ///<Used to set up the CRTP.
		friend class SFINAE;
	};



		template<typename T>
		WeightedGraph::WeightedGraph(const Graph<T>& G, const std::vector<std::vector<int>>& weights): weights(weights){
			number_of_nodes = G.number_of_nodes;
			edges = G.edges;
			path = G.path;
			distance = G.distance;
			path_edges = G.path_edges;
			disconnected = G.disconnected;
			initialize(); 
		}


		WeightedGraph::WeightedGraph(const std::vector<std::vector<int>> & edges, const std::vector<std::vector<int>>& weights)
		: Graph<WeightedGraph>(edges), weights(weights) {
			initialize();
		}
		

		template<typename T>
		WeightedGraph::WeightedGraph(Graph<T>& G) : WeightedGraph(G, zeroWeight(G.edges)) { initialize(); }


	std::vector<std::vector<int>> WeightedGraph::zeroWeight(const std::vector<std::vector<int>>& e) {
		std::vector<std::vector<int>> w;
		w.resize(e.size());
		for (size_t i = 0; i < e.size(); i++) {
			w[i].resize(e[i].size());
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
				if (distance[i] == -1 || distance[i] > distance[start] + weights[start][j] + 1) {
					distance[i] = distance[start] + weights[start][j] + 1;
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
		ColoredGraph(const Graph<T>& G, std::vector<std::vector<char>>& colors): colors(colors) {
			static_cast<Graph<T>&>(*this) = G;}

		/// Writes a graph.dot file representing the colored graph (using red and blue). The nodes are unnamed points.
		void draw();
		/// Writes a graph.dot file representing the colored graph (using red and blue) and named nodes
		void draw(const std::vector<std::string>&);

		/// Writes a graph.dot file representing the colored graph and named nodes and named edges
		void draw(const std::vector<std::string>&, const std::vector<std::vector<int>>&, std::vector<std::string>&);

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
		//the monochrome edges that alternate node colors are given weight one

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
					weights[2 * i].push_back(0);

					dual_edges[2 * i + 1].push_back(2 * edges[i][j]); //blue to red
					weights[2 * i + 1].push_back(1);
				}
				else { //blue target
					dual_edges[2 * i].push_back(2 * edges[i][j] + 1); //red to blue
					weights[2 * i].push_back(1);

					dual_edges[2 * i + 1].push_back(2 * edges[i][j] + 1); //blue to blue
					weights[2 * i + 1].push_back(0);
				}
			}
		}
		dual=WeightedGraph(dual_edges,weights);
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
			if (finder==0 || finder==1) {//use red
				if (distance[i]==-1 || (distance[i] > reddistance[2 * i + finder] || (distance[i]==reddistance[2 * i + finder] && path[i].size()>redpath[2 * i + finder].size()) )) {
					path[i] = adjustpath(redpath[2 * i + finder]);
					distance[i] = reddistance[2 * i + finder];
				}
			}
			else if (finder==2 || finder==3){//use blue
				if (distance[i] == -1 || (distance[i] > bluedistance[2 * i + finder - 2] || (distance[i]==bluedistance[2 * i + finder - 2] && path[i].size()>bluepath[2*i+finder-2].size()))) {
					path[i] = adjustpath(bluepath[2 * i + finder - 2]);
					distance[i] = bluedistance[2 * i + finder - 2];
				}
			}
		}
	}

	void ColoredGraph::draw() {
		std::array<std::string, 2> coloring = { "red","blue" };
		std::ofstream file;
		file.open("graph.dot");
		file << "digraph G{ \n node[shape=point] \n";
		for (int i = 0; i < number_of_nodes; i++) {
			for (std::vector<int>::size_type j = 0; j < edges[i].size(); j++) {
				file << i << "->" << edges[i][j] << "[color=\"" << coloring[colors[i][j]] << "\"]\n";
			}
		}
		file << "}";
		file.close();
	}

	void ColoredGraph::draw(const std::vector<std::string>& names) {
		std::array<std::string, 2> coloring = { "red","blue" };
		std::ofstream file;
		file.open("graph.dot");
		file << "digraph G{ \n";
		for (int i = 0; i < number_of_nodes; i++) {
			for (std::vector<int>::size_type j = 0; j < edges[i].size(); j++) {
				file << "\"" << names[i] << "\"" << "->" << "\"" << names[edges[i][j]] << "\" [color=\"" << coloring[colors[i][j]] << "\"]\n";
			}
		}
		file << "}";
		file.close();
	}

	void ColoredGraph::draw(const std::vector<std::string>& names, const std::vector<std::vector<int>>& edgeid, std::vector<std::string>& edgenames) {
		std::ofstream file;
		file.open("graph.dot");
		file << "digraph G{ \n";
		for (int i = 0; i < number_of_nodes; i++) {
			for (std::vector<int>::size_type j = 0; j < edges[i].size(); j++) {
				if (colors[i][j] == 1 || edgeid[i][j]==-1)
					file << "\"" << names[i] << "\"" << "->" << "\"" << names[edges[i][j]] << "\" [color=\"blue\"]\n";
				else
					file << "\"" << names[i] << "\"" << "->" << "\"" << names[edges[i][j]] << "\" [color=\"red\", label =\"" << edgenames[edgeid[i][j]] << "\"]\n";
			}
		}
		file << "}";
		file.close();
	}

}
