#pragma once
#include "../Factorization.hpp"

///@file
///@brief Contains the multiplication graph and the methods for factorizing generators.

namespace mackey {
	template<typename group_t>
	Factorization<group_t>::Factorization(int level, const std::vector<int>& minsphere, const std::vector<int>& maxsphere, const std::vector<std::vector<int>>& basicIrreducibles, const std::vector<std::string>& basicIrr_names)
		: MultGraph<group_t>(level, minsphere, maxsphere, basicIrreducibles), basicIrr_names(basicIrr_names), shortest_paths(this->graph) {
		this->graph.node_names = [&](size_t i)->std::string {return getname(i);};
	}

	template<typename group_t>
	Factorization<group_t>::Factorization(const MultTableData<group_t>& MTD, const std::vector<std::string>& basicIrr_names)
		: MultGraph<group_t>(MTD), basicIrr_names(basicIrr_names), shortest_paths(this->graph) {
		this->graph.node_names = [&](size_t i)->std::string {return getname(i);};
	}

	template<typename group_t>
	Factorization<group_t>::Factorization(MultTableData<group_t>&& MTD, const std::vector<std::string>& basicIrr_names)
		: MultGraph<group_t>(std::move(MTD)), basicIrr_names(basicIrr_names), shortest_paths(this->graph) {
		this->graph.node_names = [&](size_t i)->std::string {return getname(i);};
	}

	template<typename group_t>
	void Factorization<group_t>::set_sources(const std::vector<std::vector<int>>& given_sources, const std::vector<std::string>& names) {
		sources.clear();
		sources.reserve(given_sources.size());
		for (const auto& i : given_sources) {
			auto deg = this->antidegree[i];
			auto basis = basisElement<rank_t>(1, 0);
			sources.push_back(this->antielement.at(std::make_pair(deg, basis)));
		}
		source_names.clear();
		for (size_t i = 0; i < names.size(); i++)
			source_names[sources[i]] = names[i];
	}

	template<typename group_t>
	void Factorization<group_t>::compute_with_sources(const std::vector<std::vector<int>>& given_sources, const std::vector<std::string>& names) {
		set_sources(given_sources, names);
		for (auto i : sources)
			shortest_paths.compute_with_root(i);
	}

	template<typename group_t>
	std::string Factorization<group_t>::print_power(const std::vector<size_t>& power) const {
		std::string name;
		for (size_t i = 0; i < power.size(); i++) {
			if (power[i] != 0) {
				name += basicIrr_names[i];
				if (power[i] != 1)
					name += "^" + std::to_string(power[i]);
				name += "*";
			}
		}
		if (!name.empty() && name.back() == '*')
			name.pop_back();
		return name;
	}

	namespace implementation_details {
		void append_power(std::string& name, char operation, int multiple, const std::string& power) {
			if (multiple!=1)
				name += "*" + std::to_string(multiple);
			if (!power.empty()) {
				if (operation == 0)
					name.append("*(");
				else
					name.append("/(");
				name += power + ")";
			}
		}
	}

	template<typename group_t>
	std::string Factorization<group_t>::getname(size_t i) const
	{
		std::string name;
		if (shortest_paths.root_per_path[i] == -1)
			return name;
		name.append(source_names.at(shortest_paths.root_per_path[i]));
		std::vector<size_t> powers(this->basicIrreducibles.size());
		char operation = 0;
		int multiple = 1;
		for (auto it = shortest_paths.paths[i].rbegin(); it != shortest_paths.paths[i].rend(); ++it) {
			if (it->edge().color != operation) {
				implementation_details::append_power(name, operation, multiple, print_power(powers));
				operation = it->edge().color;
				multiple = 1;
				powers = std::vector<size_t>(this->basicIrreducibles.size());
			}
			if (it->edge().color == 2)
				multiple *= it->edge().id;
			else
				powers[it->edge().id]++;
		}
		implementation_details::append_power(name, operation, multiple, print_power(powers));
		return name;
	}


	template<typename group_t>
	std::vector<std::string> Factorization<group_t>::getname(const std::vector<int>& degree) const
	{
		auto deg_index = this->getdegreeindex(degree);
		if (deg_index == -1) {
			if (!this->degreewithinrange(degree)) {
				std::cerr << "You must ensure the degree is within the computed range before passing to getname";
				abort();
			}
			else
				return { "0" };
		}
		std::vector<int> element_indices;
		for (int i = 0; i < this->number_of_generators; i++)
			if (this->tracker[i] == deg_index)
				element_indices.push_back(i);
		std::vector<std::string> names;
		names.reserve(element_indices.size());
		for (const auto& i : element_indices)
			names.push_back(getname(i));
		return names;
	}


	template<typename group_t>
	std::vector<std::vector<int>> Factorization<group_t>::disconnected_degrees() const
	{
		auto disc_indices = find_disconnected();
		std::vector<std::vector<int>> disc_degrees;
		disc_degrees.reserve(disc_indices.size());
		for (auto i : disc_indices)
			disc_degrees.push_back(this->getdegree(i));
		return disc_degrees;
	}

	template<typename group_t>
	std::vector<int> Factorization<group_t>::find_disconnected() const {
		std::vector<int> disconnected;
		for (int i = 0; i < this->number_of_generators; i++) {
			if (shortest_paths.paths[i].empty())
				disconnected.push_back(i);
		}
		return disconnected;
	}

	///Prints Factorization
	template<typename group_t>
	std::ostream& operator<<(std::ostream& os, const Factorization<group_t>& F) {
		for (size_t i = 0; i < F.number_of_generators; i++) {
			std::string name = F.getname(i);
			if (name.empty())
				name = "???";
			os << "At degree " << F.getdegree(i) << " we have: " << name << "\n";
		}
		return os;
	}
}
