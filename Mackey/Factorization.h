#pragma once
#include "Table.h"
#include "Graph.h"
#include <string>

///@file
///@brief Contains the Factorization class for factorizing generators.

namespace Mackey {

	/// Factorizes generators into the given basic irreducibles and sources
	template<typename rank_t, typename diff_t>
	class Factorization : public MultiplicationTable<rank_t, diff_t>, public ColoredGraph {
	public:
		/// The total number of generators and their multiples we can factorize
		int size;
		/// The total number of generators we can factorize
		int realsize;
		/// The names of the basic irreducibles.
		const std::vector<std::string> basicIrr_names;

		/// Retrieve the factorization of the i-th generator
		std::string getname(int i);
		/// Retrieve the degree of the i-th generator
		std::vector<int> getdegree(int i);
		/// Form the multiplication table and graph given the max and min spheres and the basic irreducibles
		Factorization(const std::vector<int>&, const std::vector<int>&, const std::vector<std::vector<int>>&, const std::vector<std::string>&);
		/// Compute the factorizations using the given sources for the multiplication graph and their given names.
		void compute_with_sources(const std::vector<std::vector<int>>&, const std::vector<std::string>&);

	private:
		void factorize(int i);
		int multiple;
		std::vector<std::vector<int>> factorization;
		std::vector<int> orderOfoperations, sources;
		std::map<int, std::string> source_names;
		void set_sources(const std::vector<std::vector<int>>&);
		void set_source_names(const std::vector<std::string>&);

	};

	template<typename rank_t, typename diff_t>
	std::vector<int> Factorization<rank_t, diff_t>::getdegree(int i) {
		std::vector<int> degree((this->antimap[i]).begin(), (this->antimap[i]).end() - 2);
		return degree;
	}


	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::set_sources(const std::vector<std::vector<int>>& given_sources) {
		sources.reserve(given_sources.size());
		for (auto i : given_sources) {
			i.push_back(0);
			i.push_back(1);
			sources.push_back(this->map[i]);
		}
	}

	template<typename rank_t, typename diff_t>
	Factorization<rank_t, diff_t>::Factorization(const std::vector<int>& minsphere, const std::vector<int>& maxsphere, const std::vector<std::vector<int>>& basicIrreducibles, const std::vector<std::string>& basicIrr_names)
		: MultiplicationTable<rank_t, diff_t>(minsphere, maxsphere, basicIrreducibles), basicIrr_names(basicIrr_names), ColoredGraph(MultiplicationTable<rank_t, diff_t>::edges, MultiplicationTable<rank_t, diff_t>::colors)
	{
		size = ColoredGraph::number_of_nodes;
		realsize = this->second_pass_nodes;
		factorization.reserve(3);
		orderOfoperations.reserve(3);
	}


	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::set_source_names(const std::vector<std::string>& names) {
		for (int i = 0; i < names.size(); i++)
		{
			source_names[sources[i]] = names[i];
		}
	}

	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::compute_with_sources(const std::vector<std::vector<int>>& given_sources, const std::vector<std::string>& names) {
		set_sources(given_sources);
		set_source_names(names);
		for (auto i : sources) {
			computeWithSource(i);
		}
	}
	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::factorize(int i) {
		if (path[i].empty()) {
			return;
		}
		std::vector<int> difference;
		std::vector<int> zerovector(prime + 1);

		factorization.clear();
		orderOfoperations.clear();
		multiple = 1;

		bool flagpos, flagneg;
		flagpos = flagneg = 0;
		std::vector<int> poscounter(this->basicIrreducibles.size());
		std::vector<int> negcounter(this->basicIrreducibles.size());

		for (int j = path[i].size() - 2; j >= 0; j--) {
			difference = this->antimap[path[i][j]] - this->antimap[path[i][j + 1]];
			auto difference_last = difference.back();
			difference.pop_back();
			difference.pop_back();

			if (difference_last != 0 && difference == zerovector) {
				multiple *= this->antimap[path[i][j]].back();
			}


			for (int k = 0; k < this->basicIrreducibles.size(); k++) {
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
				else if (difference == zerovector - this->basicIrreducibles[k]) {
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
	std::string Factorization<rank_t, diff_t>::getname(int i) {
		std::string name;
		if (this->path[i].empty()) {
			return name;
		}
		factorize(i);
		name.append(std::to_string(multiple));
		name.append("*");
		name.append(source_names[path[i].back()]);
		for (int j = 0; j < factorization.size(); j++) {
			if (orderOfoperations[j] == 1) {
				name.append("*(");
			}
			else {
				name.append("/(");
			}
			for (int k = 0; k < factorization[j].size(); k++) {
				if (factorization[j][k] > 0) {
					name.append(basicIrr_names[k]);
					if (factorization[j][k] > 1) {
						name.append("^");
						name.append(std::to_string(factorization[j][k]));
					}
					name.append("*");
				}
			}
			name.append(")");
		}
		return name;
	}
}
