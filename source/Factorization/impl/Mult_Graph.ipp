#pragma once
#include "../Mult_Graph.hpp"

///@file
///@brief Contains the multiplication graph.

namespace mackey {

	namespace implementation_details {
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

	template<typename group_t>
	std::vector<int> MultGraph<group_t>::getdegree(int i) const {
		return this->degree[this->tracker[i]];
	}

	template<typename group_t>
	Eigen::Matrix<int, 1, -1> MultGraph<group_t>::getelement(int i) const {
		return this->element[i].template cast<int>();
	}

	template<typename group_t>
	int MultGraph<group_t>::getelementindex(const std::vector<int>& deg, const rank_t& elmnt) const {
		auto deg_i = this->getdegreeindex(deg);
		if (deg_i == -1)
			return -1;
		auto it = antielement.find(std::make_pair(deg_i, elmnt));
		if (it == antielement.end())
			return -1;
		return it->second;
	}

	template<typename group_t>
	int MultGraph<group_t>::isMultiple(int i, int j) const {
		return this->NonZeroHomology.at(getdegree(i)).group.isMultiple(element[i], element[j]);
	}

	template<typename group_t>
	int MultGraph<group_t>::order(int i) const {
		return this->NonZeroHomology.at(getdegree(i)).group.order(element[i]);
	}

	template<typename group_t>
	template<typename ...Args>
	MultGraph<group_t>::MultGraph(Args&&... args): MultTable<group_t>(std::forward<Args>(args)...) {
		make();
	}

	template<typename group_t>
	void MultGraph<group_t>::make() {
		initialize();
		for (int i = 0; i < number_of_generators; i++)
			add_edges(i);
		pass_multiples();
	}

	template<typename group_t>
	void MultGraph<group_t>::initialize() {
		element.reserve(4 * this->NonZeroHomology.size());
		for (size_t i = 0; i < this->NonZeroHomology.size(); i++) {
			auto size = this->NonZeroHomology.at(this->degree[i]).group.number_of_summands();
			for (size_t j = 0; j < size; j++) {
				auto basis = basisElement<rank_t>(size, j);
				element.push_back(basis);
				antielement[std::make_pair(i, basis)] = element.size() - 1;
				tracker.push_back(i);
			}
		}
		number_of_generators = element.size();
		graph.reserve_nodes((group_t::power + 1) * number_of_generators);
		graph.adjoin_nodes(number_of_generators);
		graph.reserve_edges_per_node(2 * this->basicIrreducibles.size());
		zeroproduct.resize(number_of_generators);

		for (int i = 0; i < number_of_generators; i++) 
			zeroproduct[i].reserve(this->basicIrreducibles.size());
	}

	///Forms the edges given the desired indices and irreducibles
	template<typename group_t>
	void MultGraph<group_t>::add_edges(int i) {
		auto deg_i = tracker[i];
		for (size_t j = 0; j < this->basicIrreducibles.size(); j++) {
			auto prod = product_element_irreducible(i, j);
			if (prod.deg == -1)
				continue;
			auto v = determine_connection(i, j, prod);
			if (v.size() != 0)
				connect(v, i, j, prod.deg);
		}
	}

	template<typename group_t>
	auto MultGraph<group_t>::product_element_irreducible(int i, int j) -> Product_Element_Irreducible {
		Product_Element_Irreducible prod;
		int deg_i = this->tracker[i];
		prod.deg= this->index_product(deg_i, j);
		if (prod.deg == -1)
			return prod;
		prod.G= this->Greens[deg_i][j];
		prod.basis= prod.G.getNormalBasis(0, element[i]);
		return prod;
	}

	template<typename group_t>
	auto MultGraph<group_t>::product_element_irreducible(int i, int j1, int j2) -> Product_Element_Irreducible {
		Product_Element_Irreducible prod_all;
		auto prod_ij1 = product_element_irreducible(i, j1);
		prod_all.deg = this->index_product(prod_ij1.deg, j2);
		if (prod_all.deg == -1)
			return prod_all;
		prod_all.G = this->tripleGreens.at({ this->tracker[i], j1, j2 });
		prod_all.basis = prod_all.G.getNormalBasis(0, prod_ij1.basis);
		return prod_all;
	}

	template<typename group_t>
	auto MultGraph<group_t>::product_candidates(const Product_Element_Irreducible& prod) -> std::vector<rank_t> {
		return id_candidates(prod.basis, prod.G.boxID, this->NonZeroHomology.at(this->degree[prod.deg]));
	}




	///Connect element[i] to element[i] * basicIrreducible[j] if the latter is nonzero and can be identified (if noncyclic group)
	template<typename group_t>
	auto MultGraph<group_t>::determine_connection(int i, int j, const Product_Element_Irreducible& prod)
	{
		if (prod.basis.isZero()) {
			if (std::find(zeroproduct[i].begin(), zeroproduct[i].end(), j) == zeroproduct[i].end())
				zeroproduct[i].push_back(j);
			return rank_t();
		}
		if (prod.basis.size() == 1)
			return prod.basis;
		auto candidates = product_candidates(prod);
		if (candidates.size() == 0) //not identified
			return rank_t();
		if (candidates.size() == 1)
			return candidates[0];
		if (candidates.size() > 1) 
			unidentified.insert(std::pair(i, j));
		return rank_t();
	}


	///Connect element[i] to element[i]*basicIrreducible[j] and if multiplication is an injection, connect the other way too.
	template<typename group_t>
	void MultGraph<group_t>::connect(const rank_t& basis, int i, int j, int deg_prod)
	{
		auto k = find_and_add_element(deg_prod, basis);
		graph.draw_edge(i, k, 0, j);
		if (injection(i, k)) 
			graph.draw_edge(k, i, 1, j);
	}


	template<typename group_t>
	int MultGraph<group_t>::find_and_add_element(int deg, const typename group_t::rank_t& basis) {
		auto iterator = antielement.find(std::make_pair(deg, basis));
		if (iterator == antielement.end()) {
			auto k = element.size();
			element.push_back(basis);
			tracker.push_back(deg);
			antielement[std::make_pair(deg, basis)] = k;
			graph.adjoin_nodes(k + 1);
			graph.reserve_edges_for_node(2 * this->basicIrreducibles.size(), k);
			zeroproduct.resize(k + 1);
			zeroproduct.reserve(this->basicIrreducibles.size());
			return k;
		}
		else
			return iterator->second;
	}

	template<typename group_t>
	inline bool MultGraph<group_t>::injection(int i, int k) const {
		auto a = order(i);
		auto b = order(k);
		return ((a == 0 && b == 0) || (a != 0 && a <= b)); //injection if both have infinite order or a has finite order and b^n=1 --> a^n=1
	}

	template<typename group_t>
	std::vector<int> MultGraph<group_t>::other_elements(int i) const {
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

	template<typename group_t>
	void MultGraph<group_t>::pass_multiples() {
		for (auto i = number_of_generators; i < graph.number_of_nodes(); i++) {
			auto others = other_elements(i);
			for (const auto& j : others) {
				auto m = isMultiple(i, j);
				if (m != 0) {
					graph.draw_edge(j, i, 2, m);
					break;
				}
			}
		}
	}
}
