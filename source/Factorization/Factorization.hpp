#pragma once
#include "Mult_Graph.hpp"

///@file
///@brief Contains the multiplication graph and the methods for factorizing generators.

namespace mackey {

	/// Factorizes generators into the given basic irreducibles and sources
	template<typename group_t>
	class Factorization : public MultGraph<group_t> {

		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;
		using typename MultGraph<group_t>::graph_t;
	public:
		/// The names of the basic irreducibles.
		const std::vector<std::string> basicIrr_names;

		/// Retrieve the factorization of all elements in a given degree
		std::vector<std::string> getname(const std::vector<int>&) const;

		/// Retrieve the factorization of the i-th element
		std::string getname(size_t) const;

		/// Form the multiplication table and graph given the max and min spheres and the basic irreducibles
		Factorization(int, const std::vector<int>&, const std::vector<int>&, const std::vector<std::vector<int>>&, const std::vector<std::string>&);


		/// @brief		Constructs multiplication table given the fundamental data (say after serialization)
		Factorization(const MultTableData<group_t>& MTD, const std::vector<std::string>&);

		/// @brief		Constructs multiplication table given the fundamental data (say after serialization)
		Factorization(MultTableData<group_t>&& MTD, const std::vector<std::string>&);


		/// Compute the factorizations using the given sources for the multiplication graph and their given names.
		void compute_with_sources(const std::vector<std::vector<int>>&, const std::vector<std::string>&);

		std::vector<std::vector<int>> disconnected_degrees() const;

		///	@brief	The paths from sources to all nodes that minimize multiplication/division alteration and overall length
		MinColorsLength<graph_t> shortest_paths;
	private:
		std::string print_power(const std::vector<size_t>& power) const;
		std::vector<int> sources;
		std::map<int, std::string> source_names;
		void set_sources(const std::vector<std::vector<int>>&, const std::vector<std::string>&);
		std::vector<int> find_disconnected() const;
	};

}
#include "impl/Factorization.ipp"
