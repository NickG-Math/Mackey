#pragma once
#include "MultiplicationGraph_Identify.h"
#include "MultiplicationGraph_Connectivity.h"
#include <string>

///@file
///@brief Contains the multiplication graph and the methods for factorizing generators.

namespace Mackey {

	/// Factorizes generators into the given basic irreducibles and sources
	template<typename rank_t, typename diff_t>
	class Factorization : public MultiplicationGraphIdentify<rank_t, diff_t> {
	public:
		/// The names of the basic irreducibles.
		std::vector<std::string> basicIrr_names;

		/// Retrieve the factorization of the i-th generator
		std::string getname(int);

		/// Form the multiplication table and graph given the max and min spheres and the basic irreducibles
		Factorization(int, const std::vector<int>&, const std::vector<int>&, const std::vector<std::vector<int>>&, const std::vector<std::string>&, int = 0, bool = 0);

		/// Compute the factorizations using the given sources for the multiplication graph and their given names.
		void compute_with_sources(const std::vector<std::vector<int>>&, const std::vector<std::string>&);

		/// Uses triple box products to possibly identify ALL instances where identification failed before. Use with caution!
		void pass_unidentified(bool serialize_each_step);

		/// Uses triple box products to possibly identify the given generators.
		void pass(const std::vector<int>&, bool serialize_each_step);


		/// Uses triple box products to possibly identify the generators not connected to the sources.
		void pass_disconnected(bool serialize_each_step);




#ifdef CEREALIZE
		///Constructor given the multiplication table.
		Factorization(MultiplicationTable<rank_t, diff_t>& M, const std::vector<std::string>& basicIrr_names)
			:MultiplicationGraphIdentify<rank_t, diff_t>(M), basicIrr_names(basicIrr_names) { initialize(); }
#endif

	private:
		int multiple;
		void initialize();
		void factorize(int i);
		std::vector<std::vector<int>> factorization;
		std::vector<int> orderOfoperations, sources;
		std::map<int, std::string> source_names;
		void set_sources(const std::vector<std::vector<int>>&, const std::vector<std::string>&);
		void find_disconnected_generators();

	};


	template<typename rank_t, typename diff_t>
	Factorization<rank_t, diff_t>::Factorization(int level, const std::vector<int>& minsphere, const std::vector<int>& maxsphere, const std::vector<std::vector<int>>& basicIrreducibles, const std::vector<std::string>& basicIrr_names, int number_of_teams, bool serialize_each_step)
		: MultiplicationGraphIdentify<rank_t, diff_t>(level, minsphere, maxsphere, basicIrreducibles, number_of_teams, serialize_each_step), basicIrr_names(basicIrr_names) { initialize(); }

	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::initialize() {
		factorization.reserve(3);
		orderOfoperations.reserve(3);
	}

	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::set_sources(const std::vector<std::vector<int>>& given_sources, const std::vector<std::string>& names) {
		sources.reserve(given_sources.size());
		for (const auto& i : given_sources) {
			auto deg = this->antidegree[i];
			auto basis = basisElement<rank_t>(1, 0);
			sources.push_back(this->antielement.at(std::make_pair(deg, basis)));
		}
		for (std::vector<int>::size_type i = 0; i < names.size(); i++)
			source_names[sources[i]] = names[i];
	}

	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::compute_with_sources(const std::vector<std::vector<int>>& given_sources, const std::vector<std::string>& names) {
		set_sources(given_sources, names);
		for (const auto& i : sources)
			this->computeWithSource(i);
	}


	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::pass_unidentified(bool serialize_each_step) {
		this->pass_all_unidentified(serialize_each_step);
		for (const auto& i : sources)
			this->computeWithSource(i);
	}


	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::pass(const std::vector<int>& desired_generators, bool serialize_each_step) {
		do
			this->pass_division(desired_generators, serialize_each_step);
		while (this->can_do_more);
		do
			this->pass_product(desired_generators, serialize_each_step);
		while (this->can_do_more);
	}




	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::pass_disconnected(bool serialize_each_step) {
		find_disconnected_generators(); //find the disconnected generators
		//first use division to identify
		if (this->disconnected.size() > 0) {
			do {
				this->pass_division(this->disconnected, serialize_each_step);
				find_disconnected_generators();
			} while (this->disconnected.size() > 0 && this->can_do_more);
		}
		//next use multiplication to identify
		if (this->disconnected.size() > 0) {
			do {
				this->pass_product(this->disconnected, serialize_each_step);
				find_disconnected_generators();
			} while (this->can_do_more);
		}
	}

	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::find_disconnected_generators() {
		for (const auto& i : sources)
			this->computeWithSource(i);
		this->compute_disconnected(this->number_of_generators);
	}

	template<typename rank_t, typename diff_t>
	void Factorization<rank_t, diff_t>::factorize(int i) {
		if (this->path[i].empty())
			return;

		factorization.clear();
		orderOfoperations.clear();

		bool flagpos, flagneg;
		flagpos = flagneg = 0;
		std::vector<int> poscounter(this->number_of_irreducibles);
		std::vector<int> negcounter(this->number_of_irreducibles);

		multiple = 1;
		for (int j = this->path[i].size() - 2; j >= 0; j--) {
			int indexstart = this->tracker[this->path[i][j]];
			int indexend = this->tracker[this->path[i][j + 1]];

			auto difference = this->degree[indexstart] - this->degree[indexend];
			if (isZero(difference)) {
				int m = this->isMultiple(this->path[i][j], this->path[i][j + 1]);
				if (m != 0)
					multiple *= m;
			}
			for (int k = 0; k < this->number_of_irreducibles; k++) {
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
				else if (difference == -this->basicIrreducibles[k]) {
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