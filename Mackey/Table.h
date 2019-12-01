#pragma once
#include "Compute.h"
#include "Green.h"
///@file
///@brief Contains the classes and methods to form the multiplication table/graph.

namespace Mackey
{
	///////////////////////////////////////////////////
	/// The input for the multiplication table
	
	/// Since there may be multiple generators in each degree (noncyclic homology) there are two indices to keep track of
	/// The element index for the generators and their multiples, and the degree index
	///////////////////////////////////////////////////
	template<typename rank_t, typename diff_t>
	class TableInput {
	protected:

		const int level; ///<The Mackey functor level we are working in.

		///////////////////////////////////////
		///The nonzero homology groups in our given range and level and their identification information (if non cyclic). 

		///Order is "random" and we keep tabs on it using antidegree
		///////////////////////////////////////
		std::vector<IDGenerators<rank_t, diff_t>> NonZeroHomology;

		std::vector<Chains<rank_t, diff_t>> IndexedChains; ///<The Chains of the nonzero homology groups.

		std::vector<std::vector<int>> degree; ///< The degree of each degree index
		std::map<std::vector<int>, int> antidegree; ///< Maps degrees to their degree index

		std::vector<std::array<int, 2>> element; ///<An element in each NonZeroHomology group. Stores position and multiple
		std::vector<int> tracker; ///< Maps element index to degree index

		std::map<std::array<int, 3>, int> antielement; ///< Maps degree index and the element to the element index

		const std::vector<int> minsphere;///<The lower bound on the range of our spheres
		const std::vector<int> maxsphere;///<The upper bound on the range of our spheres

		/// Constructs the input of the multiplication table.
		TableInput(int level, const std::vector<int>& minsphere, const std::vector<int>& maxsphere);

		inline std::vector<int> getsphere(const std::vector<int>& deg) {
			return tail(deg.data(), deg.size());
		}

		///Hashes the given sphere to produce the index for IndexedChains.
		template<typename sphere_t>
		int hash(const sphere_t&);

		///Finds all element indices with the same degree index as i
		std::vector<int> antitrack(int i) const;
	};


	template<typename rank_t, typename diff_t>
	TableInput<rank_t, diff_t>::TableInput(int level, const std::vector<int>& minsphere, const std::vector<int>& maxsphere) :
		level(level), minsphere(minsphere), maxsphere(maxsphere)
	{
		auto spheres = DegreeConstruction(minsphere, maxsphere);
		auto maxlength = spheres.size() * (dimension(maxsphere - minsphere) + 1);
		NonZeroHomology.reserve(maxlength);
		degree.reserve(maxlength);
		element.reserve(maxlength);
		IndexedChains.resize(spheres.size());

		int degindex = 0;
		for (const auto& i : spheres) {
			auto C = ROChains<rank_t, diff_t>(i);
			for (int k = 0; k <= C.maxindex;k++) {
				Junction<rank_t, diff_t> J(C, k);
				internal::IDGeneratorCompute<rank_t, diff_t>ID(level, J);
				if (!ID.H_level.isZero) {
					auto length = ID.H_level.Groups.size();
					NonZeroHomology.push_back(ID.ID);

					std::vector<int> deg;
					deg.reserve(1 + i.size());
					deg.push_back(invReindex(k, i));
					deg.insert(deg.end(), i.begin(), i.end());
					degree.push_back(deg);
					antidegree[deg] = degindex;

					for (int j = 0; j < length; j++) {
						element.push_back({ j,1 });
						antielement[{degindex, j, 1}] = element.size() - 1;
						tracker.push_back(degindex);
					}
					degindex++;
				}
			}
			IndexedChains[hash(i)] = std::move(C);
		}
	}

	template<typename rank_t, typename diff_t>
	template<typename sphere_t>
	int TableInput<rank_t, diff_t>::hash(const sphere_t& sphere) {
		int hash = sphere[0] - minsphere[0];
		for (size_t i = 1; i < sphere.size(); i++) {
			hash += (sphere[i] - minsphere[i]) * (maxsphere[i-1] - minsphere[i-1] + 1);
		}
		return hash;
	}

	template<typename rank_t, typename diff_t>
	std::vector<int> TableInput<rank_t, diff_t>::antitrack(int i) const {
		std::vector<int> tracked;
		auto index = tracker[i];
		for (int j = i - 1; j >= 0; j--) {
			if (index == tracker[j])
				tracked.push_back(j);
			else
				break;
		}
		for (std::vector<int>::size_type j = i + 1; j < element.size(); j++) {
			if (index == tracker[j])
				tracked.push_back(j);
			else
				break;
		}
		return tracked;
	}
}

namespace Mackey
{

	/// The Multiplication Table
	template<typename rank_t, typename diff_t>
	class MultiplicationTable : public TableInput<rank_t, diff_t> {
	protected:
		int number_of_nodes;///<The total number of nodes in the multiplication graph
		int number_of_generators;///<The number of generators in the multiplication graph
		///////////////////////////////////////////////
		///The edges of the multiplication graph

		///We connect a to ab for basic irreducible b, and also connect ab to a if multiplication by b is an injection
		///////////////////////////////////////////////
		std::vector<std::vector<int>> edges;

		///////////////////////////////////////////////
		///The colors of the edges of the multiplication graph

		///The edge from a to ab is red and the one from ab to a (if it exists) is blue
		///////////////////////////////////////////////
		std::vector<std::vector<char>> colors;

		std::vector<std::vector<int>> basicIrreducibles;///<The basic irreducibles we use to produce the factorizations (eg Euler and orientation classes)
		std::vector<Chains<rank_t, diff_t>> basicChains;///<The Chains of the basic irreducibles.

		///Constructs the multiplication table given the maximum and minimum spheres and the basic irreducibles.
		MultiplicationTable(int, const std::vector<int>&, const std::vector<int>&, const std::vector<std::vector<int>>&);

	private:
		std::vector<std::pair<Green<rank_t, diff_t>, int>> Greens;
		std::vector<int> disconnected;
		void getChains();
		std::pair<Green<rank_t, diff_t>, int> multiply(int, int);
		void connect(int, int , int, int, int);
		std::pair<int, int> determine_connection(const Green<rank_t, diff_t>&, int, int, int);
		void make();
		bool injection(int, int);
		int order(int);
		void multiply_all(const std::vector<int>& element_indices, const std::vector<int>& irreducible_indices);
		void pass(const std::vector<int>&, const std::vector<int>&);
	};

	template<typename rank_t, typename diff_t>
	MultiplicationTable<rank_t, diff_t>::MultiplicationTable(int level, const std::vector<int>& minsphere, const std::vector<int>& maxsphere, const std::vector<std::vector<int>>& basicIrreducibles)
		:TableInput<rank_t, diff_t>(level, minsphere, maxsphere), basicIrreducibles(basicIrreducibles) {
		number_of_nodes = this->element.size();
		number_of_generators = number_of_nodes;
		edges.reserve((prime + 1) * number_of_nodes); //when we need to add extra nodes
		edges.resize(number_of_nodes);
		colors.reserve((prime + 1) * number_of_nodes);
		colors.resize(number_of_nodes);
		getChains();
		make();
	};

	///Set the basicChains given the basicIrreducibles
	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::getChains() {
		basicChains.clear();
		basicChains.reserve(basicIrreducibles.size());
		for (const auto& i : basicIrreducibles) {
			basicChains.push_back(ROChains<rank_t, diff_t>(tail(i.data(), i.size())));
		}
	}

	///Forms the edges of the multiplication table.
	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::make() {

		std::vector<int> element_indices(number_of_nodes);
		std::iota(element_indices.begin(), element_indices.end(), 0);

		std::vector<int> irreducible_indices(basicIrreducibles.size());
		std::iota(irreducible_indices.begin(), irreducible_indices.end(),0);

		multiply_all(element_indices, irreducible_indices);
		pass(element_indices, irreducible_indices);

		if (!disconnected.empty())
			pass(disconnected, irreducible_indices); //try one more time since you have more edges than before.

		number_of_nodes = this->element.size();

		//now time to connect each element to its multiple

		for (int i = number_of_generators; i < number_of_nodes; i++) {
			auto original_i = this->antielement.at({ this->tracker[i], 0, 1 });
			//connect through multiplication
			edges[original_i].push_back(i);
			colors[original_i].push_back(0);
		}

	}

	///Compute the products of all elements and irreducibles provided
	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::multiply_all(const std::vector<int>& element_indices, const std::vector<int>& irreducible_indices) {

		std::vector<int> degree_indices;
		degree_indices.reserve(element_indices.size());
		for (const auto& i : element_indices) {
			if (degree_indices.empty() || degree_indices.back() != i) {
				degree_indices.push_back(i);
			}
		}
		auto maxdegree_plus1 = degree_indices.back() + 1; //it's the size if consecutive
		auto maxirr_plus1 = irreducible_indices.back() + 1;
		Greens.resize(maxdegree_plus1 * maxirr_plus1);

		//#pragma omp parallel for num_threads(12)
		for (decltype(degree_indices.size())  q = 0; q < degree_indices.size(); q++) { //for OpenMP parallelization we don't use a range based loop
			for (const auto& j : irreducible_indices) {
				Greens[j + degree_indices[q] * maxirr_plus1] = multiply(degree_indices[q], j);
			}
		}
	}


	///Forms the edges given the desired indices and irreducibles
	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::pass(const std::vector<int>& element_indices, const std::vector<int>& irreducible_indices) {
		for (const auto& i:element_indices) {
			edges[i].reserve(2 * irreducible_indices.size());
			auto maxirr_plus1 = irreducible_indices.back() + 1;
			auto degree_index = this->tracker[i];
			int counter = 0;
			disconnected.clear();
			for (const auto& j : irreducible_indices) 
			{
				auto G = Greens[j + degree_index * maxirr_plus1].first;
				auto index_product = Greens[j + degree_index * maxirr_plus1].second;
				if (!G.isZero) 
				{
					auto pair = determine_connection(G, i, j, index_product);
					auto position = pair.first;
					auto multiple = pair.second;
					if (position >= 0 && multiple != 0) 
					{
						connect(i, j, index_product, position, multiple);
						counter++;
					}
				}
			}
			if (counter == 0)
				disconnected.push_back(i);
		}
	}


	///Compute the product of generators a_i*b_j
	template<typename rank_t, typename diff_t>
	std::pair<Green<rank_t, diff_t>, int> MultiplicationTable<rank_t, diff_t>::multiply(int index, int j) {

		auto degreeproduct = this->degree[index] + basicIrreducibles[j]; //we are using the degree index i here
		auto iterator = this->antidegree.find(degreeproduct);
		if (iterator == this->antidegree.end()) {
			Green<rank_t, diff_t> G;
			return std::make_pair(G, -1);
		}
		auto index_product = iterator->second;
		auto degreeD = Reindex(this->degree[index]).front();
		std::vector<int> sphere(this->degree[index].begin() + 1, this->degree[index].end());
		auto D = this->IndexedChains[this->hash(sphere)];

		auto degreeC = Reindex(basicIrreducibles[j]).front();
		auto C = basicChains[j];

		Green<rank_t, diff_t> G(C, D, this->level, degreeC, degreeD);
		return std::make_pair(G, index_product);
	}


	///Determine a_i*b_j in terms of position (if noncyclic) and multiple
	template<typename rank_t, typename diff_t>
	std::pair<int,int> MultiplicationTable<rank_t, diff_t>::determine_connection(const Green<rank_t, diff_t>& G, int i, int j, int index)
	{
		auto basis = G.normalBasis[G.select(0, this->element[i][0])];
		int multiple = 0;
		int position = 0;
		rank_t newbasis;
		if (!basis.isZero() && basis.size() == 1)
			multiple = basis[0];
		else if (!basis.isZero()){
			auto candidates = identify(basis, G.boxID, this->NonZeroHomology[index]);
			if (candidates.size() == 0)
				return std::make_pair(-1, -1);
			else if (candidates.size() == 1)
				newbasis = candidates[0];
			else{
				auto other_indices = this->antitrack(index);
				for (const auto& m : other_indices) 
				{
					if (find(edges[i], m) == -1)
						return std::make_pair(-1, -1);
				}
				int p=findBasisElement(candidates);
				if (p == -1)
					return std::make_pair(-1, -1); 
				else
					newbasis = candidates[p];
			}
			position = isBasisElement(newbasis);
			multiple = newbasis[position];
		}
		return std::make_pair(position, multiple);
	}


	///Connect a_i to a_i*b_j and a_i*b_j to a_i if we have an injection
	template<typename rank_t, typename diff_t>
	void MultiplicationTable<rank_t, diff_t>::connect(int i, int j, int index, int position, int multiple)
	{
		auto iterator = this->antielement.find({ index, position, multiple });

		if (iterator == this->antielement.end())
		{ //if it's a new element (multiple)
			auto k = this->element.size();
			this->element.push_back({ position ,multiple });
			this->tracker.push_back(index);
			this->antielement[{ index, position, multiple }] = k;
			edges.resize(k + 1);
			edges[k].reserve(2 * basicIrreducibles.size());
			colors.resize(k + 1);
			colors[k].reserve(2 * basicIrreducibles.size());
			edges[i].push_back(k); // draw edge from i to k with red color
			colors[i].push_back(0);
			if (injection(i, k))
			{
				edges[k].push_back(i); //draw the inverse edge from k to i with blue color
				colors[k].push_back(1);
			}
		}
		else
		{
			auto k = iterator->second;
			if (find(edges[i], k)==-1) {// draw edge from i to k with red color
				edges[i].push_back(k);
				colors[i].push_back(0);
			}
			if (injection(i, k))
			{
				if (find(edges[k], i)==-1){ //draw the inverse edge from k to i with blue color
					edges[k].push_back(i);
					colors[k].push_back(1);
				}
			}
		}
	}

	///Checks if the surjective map from the i'th to the k'th generators is an injection
	template<typename rank_t, typename diff_t>
	inline bool MultiplicationTable<rank_t, diff_t>::injection(int i, int k) {
		auto a = order(i);
		auto b = order(k);
		return ((a == 1 && b == 1) || (a != 1 && a <= b));
	}

	///Finds the order of the i'th generator
	template<typename rank_t, typename diff_t>
	inline int MultiplicationTable<rank_t, diff_t>::order(int i) {
		auto group = this->NonZeroHomology[this->tracker[i]].group;
		auto generator = basisVector(group.size(), this->element[i][0], this->element[i][1]);
		return Mackey::order(generator, group); //defined in Homology.h
	}

}
