#pragma once
#include "Products/Massey.hpp"
#include "Chain_Complexes/Morse.hpp"

///	@file
///	@brief		Contains the methods computing the RO(G) homology of a point
/// @details	That includes the additive (Mackey functor) and multiplicative (Green functor) structures and Massey products. 

namespace mackey {


	///Compute the dimension of a representation sphere
	template<typename T, typename deg_t = std::vector<int>>
	int dimension(const deg_t& sphere);

	///Given degree (k,sphere) returns degree (k',sphere) where k' is a reindexing of k s.t. 0<=k'<=dimension(sphere)
	template<typename T, typename deg_t = std::vector<int>>
	deg_t Reindex(deg_t degree);
	///Inverse of \ref Reindex "Reindex"
	template<typename T, typename deg_t = std::vector<int>>
	deg_t invReindex(deg_t degree);

	///Returns the standard Chains of the given sphere. Optionally stop at index i (i=-1 means nonstop i.e. i=dimension(sphere) )
	template<typename T, typename deg_t = std::vector<int>>
	auto getChains(const deg_t& sphere, int i = -1);

	///Returns the Chains of given sphere. Optionally stop at index i (i=-1 means nonstop i.e. i=dimension(sphere) )
	template<typename group_t, typename deg_t>
	auto ROChains(const deg_t& sphere, int i = -1);

	///The Mackey functor Homology of a Junction of Mackey functors
	template<typename group_t>
	auto HomologyLevels(const Levels<Junction<typename group_t::rank_t, typename group_t::diff_t>, group_t >& J);

	///Computes the homology for the given level and degree
	template<typename group_t, typename deg_t = std::vector<int> >
	auto ROHomology(int level, const deg_t& degree);

	///Computes the Mackey functor homology of the given sphere
	template<typename group_t, typename deg_t = std::vector<int>>
	auto ROHomology(const deg_t& sphere);

	///Computes the product of two generators in the RO(G) homology given their level, degrees and selections (if noncyclic). Selections are 0,0 by default.
	template<typename group_t, typename deg_t = std::vector<int>>
	auto ROGreen(int level, const deg_t& first, const deg_t& second, int selectFirst = 0, int selectSecond = 0);

	///Computes the Massey product of three generators in the RO(G) homology given their level, degrees and selections (if noncyclic). Selections are 0,0,0 by default.
	template<typename group_t, typename deg_t = std::vector<int>>
	Massey<group_t> ROMassey(int level, const deg_t& first, const deg_t& second, const deg_t& third, int selectFirst = 0, int selectSecond = 0, int selectThird = 0);
}
#include "impl/Point.ipp"
