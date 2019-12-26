#pragma once
#include "Table.h"
#include "Graph.h"
#include <string>

///@file
///@brief Contains the multiplication graph and the methods for factorizing generators.

namespace Mackey {

	/// The Multiplication Graph
	template<typename rank_t, typename diff_t>
	class MultiplicationGraph : public MultiplicationTable<rank_t,diff_t>, public ColoredGraph {
	protected:
		std::vector<std::vector<int>> edgeid;	///<The identity of each edge is the basic irreducible we multiply/divide by, recorded in this vector 

		int number_of_generators;///<The number of generators in the multiplication graph

		///Constructs the multiplication graph given the maximum and minimum spheres and the basic irreducibles.
		MultiplicationGraph(int, const std::vector<int>&, const std::vector<int>&, const std::vector<std::vector<int>>&);

		std::vector<std::array<int, 2>> element; ///<An element in each NonZeroHomology group. Stores position and multiple
		std::vector<int> tracker; ///< Maps element index to degree index

		std::map<std::array<int, 3>, int> antielement; ///< Maps degree index and the element to the element index

#ifdef CEREALIZE
		///Constructs the multiplication graph given the multiplication table.
		MultiplicationGraph(MultiplicationTable<rank_t,diff_t>& M) {
			static_cast<MultiplicationTable<rank_t, diff_t>&>(*this) = M;
			make(); 
		}
#endif

	private:
		std::vector<int> antitrack(int) const;

		std::vector<int> disconnected;
		void connect(int, int, int, int, int);
		std::pair<int, int> determine_connection(const Green<rank_t, diff_t>&, int, int, int);
		bool injection(int, int);
		int order(int);
		void make();
		void pass();
		void pass(const std::vector<int>&, const std::vector<int>&);
		void initialize();
		void pass_disconnected();
		void pass_multiples();
	};


	template<typename rank_t, typename diff_t>
	std::vector<int> MultiplicationGraph<rank_t, diff_t>::antitrack(int index) const {
		std::vector<int> tracked;
		tracked.reserve(tracker.size());
		for (int i = 0; i < tracker.size(); i++) {
			if (tracker[i] == index)
				tracked.push_back(i);
		}
		return tracked;
	}




	template<typename rank_t, typename diff_t>
	MultiplicationGraph<rank_t, diff_t>::MultiplicationGraph(int level, const std::vector<int>& minsphere, const std::vector<int>& maxsphere, const std::vector<std::vector<int>>& basicIrreducibles)
	: MultiplicationTable<rank_t,diff_t>(level, minsphere, maxsphere, basicIrreducibles) {
		make();
	}


	template<typename rank_t, typename diff_t>
	void MultiplicationGraph<rank_t, diff_t>::make() {
		initialize();
		pass();
		pass_disconnected();
		number_of_nodes = element.size();
		pass_multiples();
	}

	template<typename rank_t, typename diff_t>
	void MultiplicationGraph<rank_t, diff_t>::initialize() {

		element.reserve(2 * this->NonZeroHomology.size());
		for (int i = 0; i < this->NonZeroHomology.size(); i++) {
			for (int j = 0; j < this->NonZeroHomology[i].group.size(); j++) {
				element.push_back({ j,1 });
				antielement[{i, j, 1}] = element.size() - 1;
				tracker.push_back(i);
			}
		}
		number_of_nodes=number_of_generators = element.size();


		edges.reserve((power + 1) * number_of_nodes);
		edges.resize(number_of_nodes);

		edgeid.reserve((power + 1) * number_of_nodes);
		edgeid.resize(number_of_nodes);

		colors.reserve((power + 1) * number_of_nodes);
		colors.resize(number_of_nodes);
	}

	///Repeat pass since we have more edges now and we might be able to identify better.
	template<typename rank_t, typename diff_t>
	void MultiplicationGraph<rank_t, diff_t>::pass_disconnected() {
		if (!disconnected.empty()) {
			//all irreducibles
			std::vector<int> irreducible_indices(this->number_of_irreducibles);
			std::iota(irreducible_indices.begin(), irreducible_indices.end(), 0);
			pass(disconnected, irreducible_indices);
		}
	}

	///Connect each element to its multiple
	template<typename rank_t, typename diff_t>
	void MultiplicationGraph<rank_t, diff_t>::pass_multiples() {
		for (int i = number_of_generators; i < number_of_nodes; i++) {
			auto original_i = antielement.at({ tracker[i], 0, 1 });
			//connect through multiplication
			edges[original_i].push_back(i);
			edgeid[original_i].push_back(-1); //-1 means multiple
			colors[original_i].push_back(0);
		}
	}

	///Forms the edges given for all generators
	template<typename rank_t, typename diff_t>
	void MultiplicationGraph<rank_t, diff_t>::pass() {

		//all element indices
		std::vector<int> element_indices(number_of_nodes);
		std::iota(element_indices.begin(), element_indices.end(), 0);

		//all irreducibles
		std::vector<int> irreducible_indices(this->number_of_irreducibles);
		std::iota(irreducible_indices.begin(), irreducible_indices.end(), 0);

		pass(element_indices, irreducible_indices);
	}

	///Forms the edges given the desired indices and irreducibles
	template<typename rank_t, typename diff_t>
	void MultiplicationGraph<rank_t, diff_t>::pass(const std::vector<int>& element_indices, const std::vector<int>& irreducible_indices) {
		disconnected.clear();
		for (const auto& i : element_indices) {
			edges[i].reserve(2 * this->number_of_irreducibles);
			edgeid[i].reserve(2 * this->number_of_irreducibles);
			colors[i].reserve(2 * this->number_of_irreducibles);

			auto maxirr_plus1 = irreducible_indices.back() + 1;
			auto degree_index = tracker[i];
			int counter = 0;
			for (const auto& j : irreducible_indices)
			{
				auto G = this->Greens[j + degree_index * maxirr_plus1].first;
				auto index_product = this->Greens[j + degree_index * maxirr_plus1].second;
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



	///Determine a_i*b_j in terms of position (if noncyclic) and multiple
	template<typename rank_t, typename diff_t>
	std::pair<int, int> MultiplicationGraph<rank_t, diff_t>::determine_connection(const Green<rank_t, diff_t>& G, int i, int j, int index)
	{
		auto basis = G.normalBasis[G.select(0, element[i][0])];
		int multiple = 0;
		int position = 0;
		rank_t newbasis;
		if (!basis.isZero() && basis.size() == 1)
			multiple = basis[0];
		else if (!basis.isZero()) {
			auto candidates = id_candidates(basis, G.boxID, this->NonZeroHomology[index]);
			if (candidates.size() == 0) //not identified
				return std::make_pair(-1, -1);
			else {
				newbasis = candidates[0];
				position = isBasisElement(newbasis);
				if (position >= 0 && candidates.size() == 1)
					multiple = newbasis[position];
				else {//not basis element or multiple possibilities, but as long as all but one basis elements are hit we might be good
					auto elements_product = antitrack(index);
					auto elements_i = antitrack(tracker[i]);

					std::vector<int> nothit;
					nothit.reserve(elements_product.size());
					int k = 0;
					for (const auto& m : elements_product) {
						bool hit = 0;
						for (const auto& t : elements_i) {
							hit += (find(edges[t], m) >= 0);
						}
						if (!hit)
							nothit.push_back(k);
						k++;
					}
					if (nothit.size() == 0) {
						throw("what?");
					}
					else if (nothit.size() == 1) {
						multiple = newbasis[nothit[0]];
						if (multiple == 0)
							return std::make_pair(-1, -1);
						else {
							for (const auto& i : candidates) {
								if (i[nothit[0]] != multiple)
									return std::make_pair(-1, -1);
							}
						}
					}
					else
						return std::make_pair(-1, -1);
				}
			}
		}
		return std::make_pair(position, multiple);
	}


	///Connect a_i to a_i*b_j and a_i*b_j to a_i if we have an injection
	template<typename rank_t, typename diff_t>
	void MultiplicationGraph<rank_t, diff_t>::connect(int i, int j, int index, int position, int multiple)
	{
		auto iterator = antielement.find({ index, position, multiple });

		if (iterator == antielement.end())
		{ //if it's a new element (multiple)
			auto k = element.size();
			element.push_back({ position ,multiple });
			tracker.push_back(index);
			antielement[{ index, position, multiple }] = k;

			edges.resize(k + 1);
			edges[k].reserve(2 * this->number_of_irreducibles);

			edgeid.resize(k + 1);
			edgeid[k].reserve(2 * this->number_of_irreducibles);

			colors.resize(k + 1);
			colors[k].reserve(2 * this->number_of_irreducibles);

			edges[i].push_back(k); // draw edge from i to k
			edgeid[i].push_back(j); // multiplication by the j basic irreducible
			colors[i].push_back(0); //red color

			if (injection(i, k))
			{
				edges[k].push_back(i); //draw the inverse edge from k to i
				edgeid[k].push_back(j); // division by the j basic irreducible
				colors[k].push_back(1); // blue color
			}
		}
		else
		{
			auto k = iterator->second;
			if (find(edges[i], k) == -1) {// draw edge from i to k with red color
				edges[i].push_back(k);
				edgeid[i].push_back(j);
				colors[i].push_back(0);
			}
			if (injection(i, k))
			{
				if (find(edges[k], i) == -1) { //draw the inverse edge from k to i with blue color
					edges[k].push_back(i);
					edgeid[k].push_back(j);
					colors[k].push_back(1);
				}
			}
		}
	}

	///Checks if the surjective map from the i'th to the k'th generators is an injection
	template<typename rank_t, typename diff_t>
	inline bool MultiplicationGraph<rank_t, diff_t>::injection(int i, int k) {
		auto a = order(i);
		auto b = order(k);
		return ((a == 1 && b == 1) || (a != 1 && a <= b));
	}

	///Finds the order of the i'th generator
	template<typename rank_t, typename diff_t>
	inline int MultiplicationGraph<rank_t, diff_t>::order(int i) {
		auto group = this->NonZeroHomology[tracker[i]].group;
		auto generator = basisVector(group.size(), element[i][0], element[i][1]);
		return Mackey::order(generator, group); //defined in Homology.h
	}

	/// Factorizes generators into the given basic irreducibles and sources
	template<typename rank_t, typename diff_t>
	class Factorization : public MultiplicationGraph<rank_t, diff_t> {
	public:
		/// The total number of generators and their multiples we can factorize
		int size;
		/// The total number of generators we can factorize
		int realsize;
		/// The names of the basic irreducibles.
		std::vector<std::string> basicIrr_names;

		/// Retrieve the factorization of the i-th generator
		std::string getname(int);

		/// Retrieve the degree of the i-th generator
		std::vector<int> getdegree(int) const;

		/// Retrieve the position of the i-th generator, always 0 if cyclic
		int getposition(int i) const;

		/// Find the element with given degree, position and multiple
		int getelement(int, int, int) const;

		/// Retrieve the degree of the i-th generator
		int getdegreeindex(const std::vector<int>&) const;

		/// Form the multiplication table and graph given the max and min spheres and the basic irreducibles
		Factorization(int, const std::vector<int>&, const std::vector<int>&, const std::vector<std::vector<int>>&, const std::vector<std::string>&);

		/// Compute the factorizations using the given sources for the multiplication graph and their given names.
		void compute_with_sources(const std::vector<std::vector<int>>&, const std::vector<std::string>&);


#ifdef CEREALIZE
		///Constructor given the multiplication table.
		Factorization(MultiplicationTable<rank_t, diff_t>& M, const std::vector<std::string>& basicIrr_names) 
			:MultiplicationGraph<rank_t,diff_t>(M), basicIrr_names(basicIrr_names) { initialize(); }
#endif

	private:
		void initialize();
		void factorize(int i);
		int multiple;
		std::vector<std::vector<int>> factorization;
		std::vector<int> orderOfoperations, sources;
		std::map<int, std::string> source_names;
		void set_sources(const std::vector<std::vector<int>>&, const std::vector<std::string>&);

	};


	template<typename rank_t, typename diff_t>
	Factorization<rank_t, diff_t>::Factorization(int level, const std::vector<int>& minsphere, const std::vector<int>& maxsphere, const std::vector<std::vector<int>>& basicIrreducibles, const std::vector<std::string>& basicIrr_names)
		: MultiplicationGraph<rank_t, diff_t>(level, minsphere, maxsphere, basicIrreducibles), basicIrr_names(basicIrr_names)	{initialize();}

	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::initialize() {
		multiple = 1;
		size = this->number_of_nodes;
		realsize = this->number_of_generators;
		factorization.reserve(3);
		orderOfoperations.reserve(3);
	}

	template<typename rank_t, typename diff_t>
	inline std::vector<int> Factorization<rank_t, diff_t>::getdegree(int i) const  {
		return this->degree[this->tracker[i]];
	}

	template<typename rank_t, typename diff_t>
	inline int Factorization<rank_t, diff_t>::getposition(int i)  const  {
		return this->element[i][this->element[i].size()-2];
	}

	template<typename rank_t, typename diff_t>
	inline int Factorization<rank_t, diff_t>::getelement(int index, int pos, int multiple) const {
		auto iterator = this->antielement.find({ index,pos,multiple });
		if (iterator == this->antielement.end())
			return -1;
		else
			return iterator->second;
	}

	template<typename rank_t, typename diff_t>
	inline int Factorization<rank_t, diff_t>::getdegreeindex(const std::vector<int>& degree) const {
		auto iterator = this->antidegree.find(degree);
		if (iterator == this->antidegree.end())
			return -1;
		else
			return iterator->second;
	}

	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::set_sources(const std::vector<std::vector<int>>& given_sources, const std::vector<std::string>& names) {
		sources.reserve(given_sources.size());
		for (const auto& i : given_sources) {
			auto deg = this->antidegree[i];
			sources.push_back(this->antielement[{deg, 0, 1}]);
		}
		for (std::vector<int>::size_type i = 0; i < names.size(); i++)
			source_names[sources[i]] = names[i];
	}

	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::compute_with_sources(const std::vector<std::vector<int>>& given_sources, const std::vector<std::string>& names) {
		set_sources(given_sources,names);
		for (const auto& i : sources)
			this->computeWithSource(i);
	}

	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::factorize(int i) {
		if (this->path[i].empty())
			return;

		factorization.clear();
		orderOfoperations.clear();
		multiple = 1;

		bool flagpos, flagneg;
		flagpos = flagneg = 0;
		std::vector<int> poscounter(this->number_of_irreducibles);
		std::vector<int> negcounter(this->number_of_irreducibles);

		for (int j = this->path[i].size() - 2; j >= 0; j--) {
			int indexstart = this->tracker[this->path[i][j]];
			int indexend= this->tracker[this->path[i][j+1]];

			auto difference = this->degree[indexstart] - this->degree[indexend];
			auto difference_last = this->element[this->path[i][j]].back() - this->element[this->path[i][j+1]].back();
			if (difference_last != 0 && isZero(difference) ) {
				multiple *= this->element[this->path[i][j]].back();
			}


			for (decltype(this->number_of_irreducibles) k = 0; k < this->number_of_irreducibles; k++) {
				if (difference == this->basicIrreducibles[k]) {
					if (flagneg) {
						orderOfoperations.push_back(-1);
						factorization.push_back(negcounter); //flush
						negcounter.assign(negcounter.size(), 0);
						flagneg = 0;
					}
					poscounter[k]++;
					flagpos = 1;
					break;
				}
				else if (difference == - this->basicIrreducibles[k]) {
					if (flagpos) {
						orderOfoperations.push_back(1);
						factorization.push_back(poscounter); //flush
						poscounter.assign(poscounter.size(), 0);
						flagpos = 0;
					}
					negcounter[k]++;
					flagneg = 1;
					break;
				}
			}
		}
		if (flagpos) {
			orderOfoperations.push_back(1);
			factorization.push_back(poscounter);
		}
		else if (flagneg) {
			orderOfoperations.push_back(-1);
			factorization.push_back(negcounter);
		}
	}


	template<typename rank_t, typename diff_t>
	std::string Factorization<rank_t, diff_t>::getname(int i)
	{
		std::string name;
		if (this->path[i].empty())
			return name;
		factorize(i);
		if (multiple != 1) {
			name.append(std::to_string(multiple));
			name.append("*");
		}
		name.append(source_names[this->path[i].back()]);
		for (std::vector<int>::size_type j = 0; j < factorization.size(); j++) {
			if (orderOfoperations[j] == 1)
				name.append("*(");
			else 
				name.append("/(");
			for (std::vector<int>::size_type k = 0; k < factorization[j].size(); k++) {
				if (factorization[j][k] > 0) {
					name.append(basicIrr_names[k]);
					if (factorization[j][k] > 1) {
						name.append("^");
						name.append(std::to_string(factorization[j][k]));
					}
					name.append("*");
				}
			}
			if (name.back() == '*')
				name.pop_back();
			name.append(")");
		}
		return name;
	}
}
