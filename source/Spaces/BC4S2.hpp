#pragma once
#include "Space.hpp"
#include "Utility/Generators.hpp"

///	@file
///	@brief	Contains the class \ref mackey::BC4S2

namespace mackey {

	///	@brief		The space \f$B_{C_4}\Sigma_2\f$ up to given dimension
	///	@attention	Currently the differentials are only supported for \f$\mathbf Z/2\f$ coefficients
	template<typename group_t>
	class BC4S2 : public Space<BC4S2<group_t>, group_t> {

	public:
		///	@brief		Constructs BC4S2 as a cellular equivariant space up to given dimension.
		///	@attention	The dimension must not be 2 mod 4
		BC4S2(int dimension);

	private:
		typedef typename Space<BC4S2<group_t>, group_t>::cell_coefficient_pair ccp;
		typedef typename std::vector<ccp> boundary_t;
		typedef typename Space<BC4S2<group_t>, group_t>::equivariant_index eq_t;
		typedef typename Space<BC4S2<group_t>, group_t>::nonequivariant_index neq_t;
		typedef typename group_t::rank_t rank_t;
		typedef typename group_t::diff_t diff_t;

		struct cell_expanded {
			std::vector<char> cell; ///<Vector of 0,1 entries representing the homogeneous coordinates
			cell_expanded action() const;
			/// @brief	Equality mod \f$Sigma_2\f$, ignoring the placeholder coordinates 2
			bool operator==(const cell_expanded& c) const;
			bool isfaceof(const cell_expanded& c) const;
		};
		struct cell_compressed {
			std::vector<uint32_t> cell; ///< Vector whose entries represent the nonzero homogeneous coordinates
			bool operator==(const cell_compressed& c) const;
		};

		const int dimension;
		std::vector<std::vector<cell_compressed>> cells_equivariant;
		std::vector<std::vector<std::vector<int>>> cuts;
		std::vector<std::vector<cell_expanded>> cells_nonequivariant;

		int getrank(const cell_compressed& v) const;
		rank_t getrank(const std::vector<cell_compressed>& v) const;

		void setcells();
		void setexpandedcells();

		std::vector<cell_expanded> getfaces(int dim, eq_t k);

		cell_expanded getexpandedcell(int dim, eq_t i) const;

		int64_t estimate_nonzero_entries(int i);
		auto boundary(int dim, eq_t k);

		template<typename T>
		std::pair<cell_compressed, int> erase_by_criterion(const cell_compressed& v, const T& criterion) const;

		std::array<std::pair<cell_compressed, int>, 2> quickcut(const cell_compressed& v) const;
		std::pair<std::vector<cell_compressed>, std::vector<std::vector<int>>> quickcut(const std::vector<cell_compressed>& v) const;

		///Befriending parent class
		friend class Space<BC4S2<group_t>, group_t>;
	};

}
#include "impl/BC4S2.ipp"
