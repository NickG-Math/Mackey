#pragma once
#include "Point.hpp"

///	@file
///	@brief Contains the class \ref mackey::Space

namespace mackey {

	///	@brief			An equivariant space with an equivariant cellular decomposition. 
	///	@tparam	space_t	The specialization class type eg \c BC4S2
	///	@tparam	group_t	The group type eg \c C4 
	///	@details		User must inherit from this class to use it (CRTP). 
	///					They must also construct \ref cells, and specialize the functions
	///					\ref boundary and \ref estimate_nonzero_entries
	template<typename space_t, typename group_t>
	class Space {
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;
	public:

		///Returns const reference to the equivariant cellular chain complex of the space
		const auto& getChains();

		///Returns const reference to the equivariant cellular cochain complex of the space
		const auto& getCoChains();
		///Returns \f$H_{*+V}^{C_{2^n}}(X)\f$ where \f$l\f$=level and \f$V\f$=homologysphere are provided
		auto ROHomology(int level, const std::vector<int>& homologysphere);

		///Returns \f$H_{*+V}(X)\f$ where \f$V\f$=homologysphere is provided
		auto ROHomology(const std::vector<int>& homologysphere);

		///Returns \f$H^{*+V}_{C_{2^n}}(X)\f$ where \f$n\f$=level and \f$V\f$=cohomologysphere are provided
		auto ROCohomology(int level, const std::vector<int>& cohomologysphere);

		///Returns \f$H^{*+V}(X)\f$ where \f$V\f$=cohomologysphere is provided
		auto ROCohomology(const std::vector<int>& cohomologysphere);

	protected:

		///	@brief The equivariant cells in each dimension. Needs to be constructed in any child class!
		///	@details Eg for \f$G=C_4\f$, cell[1]= [2,4,1] means that the 1-dimensional cells form: \f$(C_4/C_2 \coprod C_4/e \coprod C_4/C_4)_+\wedge S^1\f$
		std::vector<rank_t> cells;

		///	@brief A pair of a cell  and a coefficient.
		///	@details The cell is given by its index
		class cell_coefficient_pair {
			int cell;
			scalar_t<diff_t> coefficient; 
		public:
			///Constructor setting cell and coefficient directly
			cell_coefficient_pair(int cell, scalar_t<diff_t> coefficient) :cell(cell), coefficient(coefficient) {}
		};

		///	@brief	 					The boundary map of a cell. Must be overriden in any child class!
		///	@details 					The boundary of the j-th i-dimensional cell is the linear combination of i-1-dimensional cells.
		///	@param	dimension_domain	The dimension of the domain cell
		///	@param	cell				The index of the domain cell in the equivariant decomposition. 
		///	@return						Up to the specialization:
		///								It must be a class with a const iterator returning a \ref Space::cell_coefficient_pair
		///								(eg a vector of \ref Space::cell_coefficient_pair)
		auto boundary(int dimension_domain, int cell);

		///	@brief An estimate of the nonzero entries of the differential matrix at the given dimension. Must be overriden in any child class!
		///	@details If no estimate is desired, just return -1 in any child class
		long estimate_nonzero_entries(int dimension_domain);

	private:

		auto getROChains(const std::vector<int>& sphere);
		auto getROCoChains(const std::vector<int>& sphere);
		Chains<rank_t,diff_t> complex, cocomplex;

		auto getdiff(int domain);

		void setcomplexes();

	};

}
#include "impl/Space.ipp"
