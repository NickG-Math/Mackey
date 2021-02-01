#pragma once
#include "Factorization/MultiplicationGraph_Identify.h"
#include "Factorization/MultiplicationGraph_Connectivity.h"

///@file
///@brief Contains the multiplication graph and the methods for factorizing generators.

namespace {
	struct factorization_triple {
		int multiple;
		std::vector<std::vector<int>> factorization;
		std::vector<int> orderOfoperations;
	};
}

namespace Mackey {

	//forward declaration
	template<typename rank_t, typename diff_t>
	class FactorizationPrinter;

	/// Factorizes generators into the given basic irreducibles and sources
	template<typename rank_t, typename diff_t>
	class Factorization : public MultiplicationGraphIdentify<rank_t, diff_t> {
	public:
		/// The names of the basic irreducibles.
		const std::vector<std::string> basicIrr_names;

		/// Retrieve the factorization of all elements in a given degree
		std::vector<std::string> getname(const std::vector<int>&) const;

		/// Retrieve the factorization of the i-th element
		std::string getname(int) const;

		/// Get reference to underlying colored graph
		const ColoredGraph& getgraph() const {
			return static_cast<const ColoredGraph&>(*this);
		}

		/// Returns FactorizationPrinter object that's ready to be read to an output stream
		FactorizationPrinter<rank_t, diff_t> generator_names() const {
			return FactorizationPrinter<rank_t, diff_t>(*this);
		}

		/// Form the multiplication table and graph given the max and min spheres and the basic irreducibles
		Factorization(int, const std::vector<int>&, const std::vector<int>&, const std::vector<std::vector<int>>&, const std::vector<std::string>&, int = 0, bool = 0);

		/// Compute the factorizations using the given sources for the multiplication graph and their given names.
		void compute_with_sources(const std::vector<std::vector<int>>&, const std::vector<std::string>&);

		/// Uses triple box products to try and identify ALL instances where identification failed before. Use with caution!
		void pass_unidentified(bool serialize_each_step=0);

		/// Uses triple box products to try and identify the given generators.
		void pass(const std::vector<int>&, bool serialize_each_step=0);


		/// Uses triple box products to try and identify the generators not connected to the sources.
		void pass_disconnected(bool serialize_each_step=0);

#ifdef CEREALIZE
		///Constructor given the multiplication table.
		Factorization(MultiplicationTable<rank_t, diff_t>& M, const std::vector<std::string>& basicIrr_names)
			:MultiplicationGraphIdentify<rank_t, diff_t>(M), basicIrr_names(basicIrr_names) {}
#endif

	private:
		factorization_triple factorize(int i) const;
		std::vector<int> sources;
		std::map<int, std::string> source_names;
		void set_sources(const std::vector<std::vector<int>>&, const std::vector<std::string>&);
		void find_disconnected_generators();
	};


	template<typename rank_t, typename diff_t>
	Factorization<rank_t, diff_t>::Factorization(int level, const std::vector<int>& minsphere, const std::vector<int>& maxsphere, const std::vector<std::vector<int>>& basicIrreducibles, const std::vector<std::string>& basicIrr_names, int number_of_teams, bool serialize_each_step)
		: MultiplicationGraphIdentify<rank_t, diff_t>(level, minsphere, maxsphere, basicIrreducibles, number_of_teams, serialize_each_step), basicIrr_names(basicIrr_names) {}

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
	factorization_triple Factorization<rank_t, diff_t>::factorize(int i) const {
		factorization_triple f;
		if (this->path[i].empty())
			return f;
		f.factorization.reserve(8); //usually good
		f.orderOfoperations.reserve(8); //usually good
		bool flagpos, flagneg;
		flagpos = flagneg = 0;
		std::vector<int> poscounter(this->number_of_irreducibles);
		std::vector<int> negcounter(this->number_of_irreducibles);
		f.multiple = 1;
		for (int j = this->path[i].size() - 2; j >= 0; j--) {
			int indexstart = this->tracker[this->path[i][j]];
			int indexend = this->tracker[this->path[i][j + 1]];

			auto difference = this->degree[indexstart] - this->degree[indexend];
			if (isZero(difference)) {
				int m = this->isMultiple(this->path[i][j], this->path[i][j + 1]);
				if (m != 0)
					f.multiple *= m;
			}
			for (int k = 0; k < this->number_of_irreducibles; k++) {
				if (difference == this->basicIrreducibles[k]) {
					if (flagneg) {
						f.orderOfoperations.push_back(-1);
						f.factorization.push_back(negcounter); //flush
						negcounter.assign(negcounter.size(), 0);
						flagneg = 0;
					}
					poscounter[k]++;
					flagpos = 1;
					break;
				}
				else if (difference == -this->basicIrreducibles[k]) {
					if (flagpos) {
						f.orderOfoperations.push_back(1);
						f.factorization.push_back(poscounter); //flush
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
			f.orderOfoperations.push_back(1);
			f.factorization.push_back(poscounter);
		}
		else if (flagneg) {
			f.orderOfoperations.push_back(-1);
			f.factorization.push_back(negcounter);
		}
		return f;
	}


	template<typename rank_t, typename diff_t>
	std::string Factorization<rank_t, diff_t>::getname(int i) const
	{
		std::string name;
		if (this->path[i].empty())
			return name;
		auto f=factorize(i);
		if (f.multiple != 1) {
			name.append(std::to_string(f.multiple));
			name.append("*");
		}
		name.append(source_names.at(this->path[i].back()));
		for (std::vector<int>::size_type j = 0; j < f.factorization.size(); j++) {
			if (f.orderOfoperations[j] == 1)
				name.append("*(");
			else
				name.append("/(");
			for (std::vector<int>::size_type k = 0; k < f.factorization[j].size(); k++) {
				if (f.factorization[j][k] > 0) {
					name.append(basicIrr_names[k]);
					if (f.factorization[j][k] > 1) {
						name.append("^");
						name.append(std::to_string(f.factorization[j][k]));
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


	template<typename rank_t, typename diff_t>
	std::vector<std::string> Factorization<rank_t, diff_t>::getname(const std::vector<int>& degree) const
	{
		auto deg_index=this->getdegreeindex(degree);
		if (deg_index == -1)
			return { "0" };
		std::vector<int> element_indices;
		for (int i = 0; i < this->tracker.size(); i++)
			if (this->tracker[i] == deg_index)
				element_indices.push_back(i);
		std::vector<std::string> names;
		names.reserve(element_indices.size());
		for (const auto& i : element_indices) {
			auto s = getname(i);
			if (s.empty())
				names.push_back("?");
			else
				names.push_back(s);
		}
		return names;
	}

	///Class that facilitates printing of all factorizations produced by Factorization
	template<typename rank_t, typename diff_t>
	class FactorizationPrinter {
		const Factorization<rank_t, diff_t>& F; 
	public:
		///Constructor taking factorization object
		FactorizationPrinter(const Factorization<rank_t, diff_t>& F) :F(F) {}
		///Returns factorization of i-th generator
		std::string operator[](int i) const {
			auto name = F.getname(i);
			if (name.empty())
				return "?";
			return name;
		}
		///Returns the number of all generators
		int size() const {
			return F.number_of_generators;
		}
		///Returns the degree of the i-th generator
		std::vector<int> getdegree(int i) const {
			return F.getdegree(i);
		}
	};

	///Prints Factorization
	template<typename rank_t, typename diff_t>
	std::ostream& operator<<(std::ostream& os, const FactorizationPrinter<rank_t, diff_t>& F) {
		for (int i = 0; i < F.size(); i++) {
			os << "At degree " << F.getdegree(i) << " we have: " << F[i] << "\n";
		}
		return os;
	}
}
