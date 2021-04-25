#pragma once
#include "Point.hpp"

///	@file
///	@brief Contains the class \ref mackey::Space

namespace mackey {

	///	@brief			An equivariant space with an equivariant cellular decomposition. 
	///	@tparam	space_t	The specialization class type eg \ref BC4S2
	///	@tparam	group_t	The group type eg \ref C4 
	///	@details		User must inherit from this class to use it (CRTP). 
	///					They must also construct \ref cells, and specialize the functions
	///					\ref boundary and \ref estimate_nonzero_entries
	///	@par			Here's how cells are encoded:
	///					Assume we have in dimension 2 the cells \f[(C_4/C_2 \coprod C_4/e \coprod C_4/C_4)_+\wedge S^2\f] 
	///					This is encoded as ```cell[2]=[2,4,1]``` \n
	///					To refer to the individual cells there are two types of indices, equivariant and nonequivariant.
	///					- The equivariant index 1 corresponds to \f$C_4/C_e\f$.
	///					- The nonequivariant index 1 corresponds to \f$gC_4/C_2\f$. \n
	///					Nonequivariant indices may be needed in the output of the boundary map eg if \f$d(x)=gy\f$ and not \f$y\f$.
	///					The boundary map takes equivariant indices and returns nonequivariant ones
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

		///	@brief		The equivariant index of a cell
		///	@details	Eg given \c [2,4,1] then index \c 1 means \f$C_4/e\f$
		struct equivariant_index {
			int index; ///<index number
			equivariant_index() = default; ///<Default constructor
			/// @brief Direct initialization
			equivariant_index(int index);
		};

		///	@brief		The nonequivariant index of a cell
		///	@details	Eg given ```[2,4,1]``` then index \c 1 means \f$gC_4/C_2\f$
		struct nonequivariant_index {
			int index;	///<index number
			nonequivariant_index() = default; ///<Default constructor
			/// @brief Direct initialization
			nonequivariant_index(int index);
		};

		/// @brief		Converts nonequivariant index to equivariant one given dimension.
		///	@warning	Only works if \c n was an equivariant index to begin with. Otherwise aborts with error
		equivariant_index convert(nonequivariant_index n, int dimension) const;

		/// @brief	Converts equivariant index to nonequivariant one given cell
		nonequivariant_index  convert(equivariant_index e, int dimension) const;

		///	@brief		A pair of a cell and a coefficient.
		/// @details	The cell has \c nonequivariant_index (use \ref convert if you have \c equivariant_index as input)
		struct cell_coefficient_pair {
			nonequivariant_index cell;
			scalar_t<diff_t> coefficient;
		};

		///	@brief	 					The boundary map of each cell.
		///	@attention					Must be overriden in any child class!
		///	@details 					The boundary of an i-dimensional cell is the linear combination of i-1-dimensional cells.
		///	@param	dimension_domain	The dimension of the domain cell
		///	@param	cell				The nonequivariant index of the "starting" cell. Eg if we have [2,4,1] then the starting cells have indices 0,2,6
		///	@return						Up to the specialization:
		///								It must be a class with a const iterator returning a \ref Space::cell_coefficient_pair
		///								(eg a vector of \ref Space::cell_coefficient_pair)
		auto boundary(int dimension_domain, equivariant_index cell);

		///	@brief An estimate of the nonzero entries of the differential matrix at the given dimension. Must be overriden in any child class!
		///	@details If no estimate is desired, just return -1 in any child class
		int64_t estimate_nonzero_entries(int dimension_domain);

		/// @brief	Produces g*n given cell n
		nonequivariant_index action(nonequivariant_index n, int dimension) const;

	private:

		auto getROChains(const std::vector<int>& sphere);
		auto getROCoChains(const std::vector<int>& sphere);
		diff_t getdiff(int domain);
		void setcomplexes();
		Chains<rank_t, diff_t> complex, cocomplex;


	};

}
#include "impl/Space.ipp"
