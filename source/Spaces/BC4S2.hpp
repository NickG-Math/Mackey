#pragma once
#include "Space.hpp"

///	@file
///	@brief	Contains the class \ref mackey::BC4S2

namespace mackey{

///	@brief	The spaces \f$B_{C_4}\Sigma_2(j)\f$
template<typename group_t>
class BC4S2 : public mackey::Space<BC4S2<group_t>, group_t> {

public:
	///Constructs BC4S2(j) as a cellular equivariant space.
	BC4S2(int j);

private:
	typedef typename group_t::rank_t rank_t;
	typedef typename group_t::diff_t diff_t;
	typedef std::vector<char> cell_t; // a cell is 1,-1,0,-1,... etc.
	typedef typename std::vector<typename Space<BC4S2<group_t>, group_t>::cell_coefficient_pair> boundary_t;

	std::vector<std::vector<cell_t>> all_cells_in_dimension; //All cells in each dimension
	std::vector<std::map<cell_t, int>> cell_to_location; //A map that takes each cell to its index

	///Implementing the boundary map of Space
	auto boundary(int i, int k);

	///Implementing the estimate map of Space
	long estimate_nonzero_entries(int i);

	///Befriending parent class
	friend class Space<BC4S2<group_t>, group_t>;
};
}
#include "impl/BC4S2.ipp"
