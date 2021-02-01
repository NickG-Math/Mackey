#pragma once
#include "Table.h"
#include "Utility/Graph.h"
#include <string>

///@file
///@brief Contains the multiplication graph.

namespace {
	///Comparator to be used for the std::map antielement below. T is a row Eigen matrix
	struct compareElements {
		template<typename T>
		bool operator()(const std::pair<int, T>& a, const std::pair<int, T>& b) const {
			if (a.first < b.first)
				return 1;
			else if (a.first > b.first)
				return 0;
			else {
				for (int i = 0; i < a.second.size(); i++) {
					if (a.second[i] < b.second[i])
						return 1;
					else if (a.second[i] > b.second[i])
						return 0;
				}
			}
			return 0;
		}
	};
}

namespace Mackey {


	/////////////////////////////////////////////////////////
	/// The Multiplication Graph created from the Multiplication Table

	/// A node is an element in a nonzero homology group, a linear combination of the generators. At initialization we only use generators, but as products of generators might not be generators,
	/// after that initial stage we add multiples and linear combinations into the mix.
	/// We have a new index for the elements, and use various vectors and ordered maps to keep track of the element index and the degree index. 
	/////////////////////////////////////////////////////////
	template<typename rank_t, typename diff_t>
	class MultiplicationGraph : public MultiplicationTable<rank_t, diff_t>, public ColoredGraph {
	public:
		int number_of_generators;///<The number of generators in the multiplication graph

		//////////////////////////////////////
		///The identity of each edge is the basic irreducible we multiply / divide by, recorded in this vector
	
		///edgeid[i][s]=j>=0 means that edge[i][s]=k connects i to k by multiplication by j. 
		///
		///edgeid[i][s]=-1 means that we have a multiple.
		///
		///edgeid[i][s]=-2-j means that edge[i][s]=k connects i to k by division by j. 
		/////////////////////////////////////
		std::vector<std::vector<int>> edgeid;


		/// Retrieve the degree of the i-th generator
		inline std::vector<int> getdegree(int i) const {
			return this->degree[this->tracker[i]];
		}

		/// Retrieve the element the i-th generator corresponds to. 
		Eigen::Matrix<int, 1, -1> getelement(int i) const { 
			return this->element[i].template cast<int>(); 
		}

		/// Retrieve the element index of the given degree and element
		int getelementindex(const std::vector<int>& deg, const rank_t& elmnt) const {
			auto deg_i = this->getdegreeindex(deg);
			auto it = antielement.find(std::make_pair(deg_i, elmnt));
			if (it == antielement.end())
				return -1;
			return it->second;
		}

	protected:
		///Constructs the multiplication graph given the maximum and minimum spheres and the basic irreducibles.
		MultiplicationGraph(int, const std::vector<int>&, const std::vector<int>&, const std::vector<std::vector<int>>&, int=0, bool=0);

		std::vector<rank_t> element; ///<An element (linear combination of generators) in each NonZeroHomology group of the table.
		std::vector<int> tracker; ///< Maps element index to degree index

		/// Maps (degree_index, element)->element_index. The comparator first compares the degree_index and if equal then compares the elements first to last entry
		std::map<std::pair<int, rank_t>, int, compareElements> antielement; 

		///Returns k if element[i]=k*element[j]. Returns 0 if element[a] is not a multiple of element[b]
		int isMultiple(int i, int j) const {
			return Mackey::isMultiple(element[i], element[j], this->NonZeroHomology[tracker[i]].group);
		}
		
#ifdef CEREALIZE
		///Constructs the multiplication graph given the multiplication table.
		MultiplicationGraph(MultiplicationTable<rank_t, diff_t>& M) {
			static_cast<MultiplicationTable<rank_t, diff_t>&>(*this) = M;
			make();
		}
#endif

	private:
		void make();
		std::vector<std::pair<int, int>> unidentified; ///All pairs (i,j) for which we couldn't identify the product i*j
		std::vector<std::vector<int>> zeroproduct; ///<zeroproduct[i] consists of all j for which i*j leads to 0.

		void initialize();

		///Adds all edges starting and ending from/to given element
		void add_edges(int);
		///Connect each element to its multiple
		void pass_multiples();

		///Finds the index of the element given degree and presentation and adds it to the list if it doesn't exist there.
		int find_and_add_element(int, const rank_t&);

		///Determines how and if i can be connected to j*i
		rank_t determine_connection(const Green<rank_t, diff_t>&, int, int, int);

		///Connects i to j*i
		void connect(const rank_t&, int, int, int);


		std::vector<int> other_elements(int) const;

		///Checks if the surjective map from the i'th to the k'th generators is an injection
		bool injection(int, int) const;

		///Finds the order of the i'th generator
		int order(int i) const {
			return Mackey::order(element[i], this->NonZeroHomology[tracker[i]].group);
		}

		template<typename s_rank_t, typename s_diff_t>
		friend class MultiplicationGraphIdentify;

		template<typename s_rank_t, typename s_diff_t>
		friend class MultiplicationGraphConnectivity;
	};

	template<typename rank_t, typename diff_t>
	MultiplicationGraph<rank_t, diff_t>::MultiplicationGraph(int level, const std::vector<int>& minsphere, const std::vector<int>& maxsphere, const std::vector<std::vector<int>>& basicIrreducibles, int number_of_teams, bool serialize_each_step)
		: MultiplicationTable<rank_t, diff_t>(level, minsphere, maxsphere, basicIrreducibles, number_of_teams, serialize_each_step) {
		make();
		}
	
	template<typename rank_t, typename diff_t>
	void MultiplicationGraph<rank_t, diff_t>::make(){
		initialize();
		for (int i = 0; i < number_of_generators; i++)
			add_edges(i);
		pass_multiples();
	}

	template<typename rank_t, typename diff_t>
	void MultiplicationGraph<rank_t, diff_t>::initialize() {
		element.reserve(4 * this->NonZeroHomology.size());
		for (int i = 0; i < this->NonZeroHomology.size(); i++) {
			int size = this->NonZeroHomology[i].group.size();
			for (int j = 0; j < size; j++) {
				auto basis = basisElement<rank_t>(size, j);
				element.push_back(basis);
				antielement[std::make_pair(i, basis)] = element.size() - 1;
				tracker.push_back(i);
			}
		}
		number_of_nodes = number_of_generators = element.size();

		edges.reserve((power + 1) * number_of_nodes);
		edges.resize(number_of_nodes);

		edgeid.reserve((power + 1) * number_of_nodes);
		edgeid.resize(number_of_nodes);

		colors.reserve((power + 1) * number_of_nodes);
		colors.resize(number_of_nodes);

		zeroproduct.resize(number_of_generators);

		for (int i = 0; i < number_of_generators; i++) {
			edges[i].reserve(2 * this->number_of_irreducibles);
			edgeid[i].reserve(this->number_of_irreducibles);
			colors[i].reserve(2 * this->number_of_irreducibles);
			zeroproduct[i].reserve(this->number_of_irreducibles);
		}
	}

	///Forms the edges given the desired indices and irreducibles
	template<typename rank_t, typename diff_t>
	void MultiplicationGraph<rank_t, diff_t>::add_edges(int i) {
		auto deg_i = tracker[i];
		for (int j = 0; j < this->number_of_irreducibles; j++) {
			auto deg_prod = this->index_product(deg_i, j);
			if (deg_prod == -1)
				continue;
			auto G = this->Greens[deg_i][j];
			auto v = determine_connection(G, i, j, deg_prod);
			if (v.size() != 0)
				connect(v, i, j, deg_prod);
		}
	}

	///Connect element[i] to element[i] * basicIrreducible[j] if the latter is nonzero and can be identified (if noncyclic group)
	template<typename rank_t, typename diff_t>
	rank_t MultiplicationGraph<rank_t, diff_t>::determine_connection(const Green<rank_t, diff_t>& G, int i, int j, int deg_prod)
	{
		auto basis = G.getNormalBasis(0,element[i]); // j*i
		if (basis.isZero()) {
			if (find(zeroproduct[i], j) == -1)
				zeroproduct[i].push_back(j);
			return rank_t();
		} 
		if (basis.size() == 1)
			return basis;
		auto candidates = id_candidates(basis, G.boxID, this->NonZeroHomology[deg_prod]);
		if (candidates.size() == 0) //not identified
			return rank_t();
		if (candidates.size() == 1)
			return candidates[0];
		if (candidates.size() > 1) {
			auto ij=std::make_pair(i, j);
			if (find(unidentified, ij) == -1)
				unidentified.push_back(ij);
		}
		return rank_t(); 
	}


	///Connect element[i] to element[i]*basicIrreducible[j] and if multiplication is an injection, connect the other way too.
	template<typename rank_t, typename diff_t>
	void MultiplicationGraph<rank_t, diff_t>::connect(const rank_t& basis, int i, int j, int deg_prod)
	{
		auto k = find_and_add_element(deg_prod, basis);
		if (find(edges[i], k) == -1) {// we don't need to find if k=element.size() as it's guaranteed to not be there
			edges[i].push_back(k); // draw edge from i to k
			edgeid[i].push_back(j); // multiplication by the j basic irreducible
			colors[i].push_back(0); //red color
		}
		if (injection(i, k))
		{
			if (find(edges[k], i) == -1) { //draw the inverse edge from k to i with blue color
				edges[k].push_back(i); //draw the inverse edge from k to i
				edgeid[k].push_back(-2 - j); // division by the j basic irreducible, -1 is reserved for multiples
				colors[k].push_back(1); // blue color
			}
		}
	}

	template<typename rank_t, typename diff_t>
	int MultiplicationGraph<rank_t, diff_t>::find_and_add_element(int deg, const rank_t& basis) {
		auto iterator = antielement.find(std::make_pair(deg, basis));
		if (iterator == antielement.end()) {
			auto k = element.size();
			element.push_back(basis);
			tracker.push_back(deg);
			antielement[std::make_pair(deg, basis)] = k;
			number_of_nodes++;

			edges.resize(k + 1);
			edges[k].reserve(2 * this->number_of_irreducibles);

			edgeid.resize(k + 1);
			edgeid[k].reserve(2 * this->number_of_irreducibles + 1);

			colors.resize(k + 1);
			colors[k].reserve(2 * this->number_of_irreducibles);

			zeroproduct.resize(k + 1);
			zeroproduct.reserve(this->number_of_irreducibles);

			return k;
		}
		else
			return iterator->second;
	}

	template<typename rank_t, typename diff_t>
	inline bool MultiplicationGraph<rank_t, diff_t>::injection(int i, int k) const {
		auto a = order(i);
		auto b = order(k);
		return ((a == 0 && b == 0) || (a != 0 && a <= b)); //injection if both have infinite order or a has finite order and b^n=1 --> a^n=1
	}

	template<typename rank_t, typename diff_t>
	std::vector<int> MultiplicationGraph<rank_t, diff_t>::other_elements(int i) const {
		std::vector<int> others;
		others.reserve(tracker.size());
		for (int j = i - 1; j >= 0; j--) {
			if (tracker[j] == tracker[i])
				others.push_back(j);
		}
		for (int j = i + 1; j < tracker.size(); j++) {
			if (tracker[j] == tracker[i])
				others.push_back(j);
		}
		return others;
	}

	template<typename rank_t, typename diff_t>
	void MultiplicationGraph<rank_t, diff_t>::pass_multiples() {
		for (int i = number_of_generators; i < number_of_nodes; i++) {
			auto others = other_elements(i);
			for (const auto& j : others) {
				int m = isMultiple(i, j);
				if (m != 0) {
					edges[j].push_back(i);
					edgeid[j].push_back(-1); //-1 means multiple
					colors[j].push_back(0);
					break;
				}
			}
		}
	}
}
