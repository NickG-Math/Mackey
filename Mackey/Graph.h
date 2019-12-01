#pragma once
#include<vector>
#include<array>
#include<fstream>
#include<queue>

///@file
///@brief Contains weighted/colored graph classes and Dikjstra algorithms.

namespace Mackey {

	//////////////////////////////////////
	/// A directed Graph

	/// This is a general class that can only be used by inheriting from it.
	/// T is the type of the graph eg weightedGraph or coloredGraph. It's used for the CRTP compile-time polymorphism.
	//////////////////////////////////////
	template<typename T>
	class Graph {
	public:
		int number_of_nodes; ///<The number of vertices

		//////////////////////////////////////
		/// The "shortest" path from a selected source to any other point. 

		/// The meaning of shortest is determined in the specializations of Graph.
		//////////////////////////////////////
		std::vector<std::vector<int>> path;

		/////////////////////////////////////
		/// The weighted distance from a selected source to all points.

		///The type of weight will depend on the specialization of Graph.
		std::vector<int> weightedDistance;


		/// Construct given the edges of the graph
		Graph(std::vector<std::vector<int>>& edges) : number_of_nodes(edges.size()), edges(edges) { 
			path.resize(number_of_nodes); weightedDistance.assign(number_of_nodes, -1);
		}

		/// A general paradigm to compute the paths using the given source and a Dikjstra algorithm. This does not rewrite any previous path computations.
		void computeWithSource(int givensource) {
			source = givensource;
			static_cast<T*>(this)->computePath();
		};

		/// A general paradigm to compute the paths using the given source and a Dikjstra algorithm, overwritting any previous path computations.
		void computeWithSource_clear(int givensource) {
			path.clear();
			path.resize(number_of_nodes);
			weightedDistance.assign(number_of_nodes, -1);
			static_cast<T*>(this)->initialize();
			computeWithSource(givensource);
		};

		/// Writes a graph.dot file representing the graph. The nodes are unnamed.
		void draw();
		/// Writes a graph.dot file representing the graph using node names.
		void draw(const std::vector<std::string>&);

	protected:
		///Default constructor
		Graph() {};
		std::vector<std::vector<int>> edges; ///<The edges of the graph
		int source; ///<The source of the graph
	private:
		//These are defined in the inherited classes
		void initialize() {};
		void computePath() {};
	};

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

		/// Sets edges and 0 weights
		WeightedGraph(std::vector<std::vector<int>>& edges) : Graph(edges), weights(zeroWeight(edges)) { initialize(); } //we must initalize here and not in a constructor of Graph because we can't static_cast(this) in a constructor

		/// Sets edges and weights
		WeightedGraph(std::vector<std::vector<int>>& edges, std::vector<std::vector<int>>& weights) : Graph(edges), weights(weights) { initialize(); } //see above why this is here

	private:
		std::vector<std::vector<int>> weights;
		std::vector<int> distance, closest;
		std::vector<char> visited, reachable; //not bool for performance
		std::priority_queue<std::pair<int, int>> next;
		std::vector<std::vector<int>> zeroWeight(std::vector<std::vector<int>>);
		void initialize();
		void computePath();
		void stepDistance(int);
		void stepPath(int, int);
		void computeDistance();
		friend class Graph<WeightedGraph>; ///<Used to set up the CRTP.
	};



	std::vector<std::vector<int>> WeightedGraph::zeroWeight(std::vector<std::vector<int>> e) {
		std::vector<std::vector<int>> w;
		w.resize(e.size());
		for (size_t i = 0; i < e.size(); i++) {
			w[i].resize(e[i].size());
		}
		return w;
	}

	void WeightedGraph::initialize() {
		visited.assign(number_of_nodes, 0);
		distance.assign(number_of_nodes, 0);
		closest = distance;
	}


	void WeightedGraph::stepDistance(int start) {
		for (size_t j = 0; j < edges[start].size(); j++) {
			auto i = edges[start][j];
			if (!visited[i]) {
				if (weightedDistance[i] == -1 || weightedDistance[i] > weightedDistance[start] + weights[start][j] + 1) {
					weightedDistance[i] = weightedDistance[start] + weights[start][j] + 1;
					distance[i] = distance[start] + 1;
					closest[i] = start;
					next.push(std::make_pair(-weightedDistance[i], i));
				}
			}
		}
		visited[start] = 1;
		reachable[start] = 1;
	}

	void WeightedGraph::computeDistance() {
		reachable.assign(number_of_nodes, 0);
		closest[source] = source;
		weightedDistance[source] = 0;
		int start = source;
		next.push(std::make_pair(-weightedDistance[source], source));
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
				path[i].reserve(distance[i] + 1);
				path[i].push_back(i);
				stepPath(i, i);
				path[i].push_back(source);
			}
		}

	}


	/// A directed Graph with two colors
	class ColoredGraph : public Graph<ColoredGraph> {
	public:
		/// Construct given the edges and colors
		ColoredGraph(std::vector<std::vector<int>>& edges, std::vector<std::vector<char>>& colors) : Graph(edges), colors(colors) {}

		/// Writes a graph.dot file representing the colored graph (using red and blue). The nodes are unnamed points.
		void draw();
		/// Writes a graph.dot file representing the colored graph (using red and blue) and named nodes
		void draw(const std::vector<std::string>&);

	private:
		std::vector<std::vector<char>> colors;		///<The colors of the graph
		void computePath();
		void constructDual();
		std::vector<int> adjustpath(std::vector<int>& path);

		WeightedGraph dual;
		friend class  Graph<ColoredGraph>; ///<Used to set up the CRTP.
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

		auto number_of_colored_nodes = 2 * number_of_nodes;
		std::vector<std::vector<int>> monochrome_edges, weights;
		monochrome_edges.resize(number_of_colored_nodes);
		weights.resize(number_of_colored_nodes);
		for (int i = 0; i < number_of_nodes; i++) {

			monochrome_edges[2 * i].reserve(edges[i].size());
			weights[2 * i].reserve(edges[i].size());
			monochrome_edges[2 * i + 1].reserve(edges[i].size());
			weights[2 * i + 1].reserve(edges[i].size());

			for (std::vector<int>::size_type j = 0; j < edges[i].size(); j++) {
				if (!colors[i][j]) {//red target
					monochrome_edges[2 * i].push_back(2 * edges[i][j]); //red to red
					weights[2 * i].push_back(0);

					monochrome_edges[2 * i + 1].push_back(2 * edges[i][j]); //blue to red
					weights[2 * i + 1].push_back(1);
				}
				else { //blue target
					monochrome_edges[2 * i].push_back(2 * edges[i][j] + 1); //red to blue
					weights[2 * i].push_back(1);

					monochrome_edges[2 * i + 1].push_back(2 * edges[i][j] + 1); //blue to blue
					weights[2 * i + 1].push_back(0);
				}

			}
		}
		dual=WeightedGraph(monochrome_edges, weights);
	}

	void ColoredGraph::computePath() {
		constructDual();
		dual.computeWithSource(2 * source);
		auto redpath = std::move(dual.path);
		auto reddistance = std::move(dual.weightedDistance);
		dual.computeWithSource_clear(2 * source + 1);
		auto bluepath = std::move(dual.path);
		auto bluedistance = std::move(dual.weightedDistance);

		for (int i = 0; i < number_of_nodes; i++) {
			std::array<int, 3> V = { reddistance[2 * i + 1], bluedistance[2 * i], bluedistance[2 * i + 1] };
			int mindistance = reddistance[2 * i];
			int finder = 0;
			for (int j = 0; j < 2; j++) {
				if (V[j] != -1 && mindistance == -1) {
					mindistance = V[j];
					finder = j + 1;
				}
				else if (V[j] != -1 && V[j] < mindistance) {
					mindistance = V[j];
					finder = j + 1;
				}
			}
			if (finder <= 1) {
				if (weightedDistance[i]==-1 || (weightedDistance[i] > reddistance[2 * i + finder] && reddistance[2 * i + finder]>=0 ) || (weightedDistance[i]==reddistance[2 * i + finder] && path[i].size()>redpath[2 * i + finder].size() && reddistance[2 * i + finder] >= 0)) {
					path[i] = adjustpath(redpath[2 * i + finder]);
					weightedDistance[i] = reddistance[2 * i + finder];
				}
			}
			else {
				if (weightedDistance[i] == -1 || (weightedDistance[i] > bluedistance[2 * i + finder - 2] && bluedistance[2 * i + finder - 2]>=0) || (weightedDistance[i]==bluedistance[2 * i + finder - 2] && path[i].size()>bluepath[2*i+finder-2].size()&& bluedistance[2 * i + finder - 2] >= 0)) {
					path[i] = adjustpath(bluepath[2 * i + finder - 2]);
					weightedDistance[i] = bluedistance[2 * i + finder - 2];
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
}
