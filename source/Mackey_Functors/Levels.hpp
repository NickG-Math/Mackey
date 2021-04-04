#pragma once
#include "Homology/Homology.hpp" 

///@file
///@brief Contains the methods for transfering, restricting and Weyl-group-acting.

namespace mackey {

	///Transfer the rank to given level.
	template<typename group_t>
	auto transfer(const typename group_t::rank_t&, int);

	///Transfer the differential to given level. We need the ranks both at the original level and the given level, for both domain and range.
	template<typename group_t>
	auto transfer(const typename group_t::diff_t& diff, const typename group_t::rank_t& domain, typename group_t::rank_t& domain_top, const typename group_t::rank_t& range, typename group_t::rank_t& range_top, int level);

	/////////////////////////////////////////////////
	/// Storage of various levels of Junction/Chains

	/// T is any type with a transfer function. Example: Junction, Chains.
	//////////////////////////////////////
	template<typename T, typename group_t>
	struct Levels {
		std::vector<T> level; ///<The various levels.
		Levels(T& bottom); 	///Transfer the bottom and get all levels.
	};

	template<typename group_t>
	chains_t<group_t> transfer(const chains_t<group_t>& C, int level);

	///Transfer Junction to the desired level.
	template<typename group_t>
	junction_t<group_t> transfer(const junction_t<group_t>& J, int level);

	///Transfer generator to level given the ranks at the original level (domain) and the target level (range).
	template<typename group_t, typename Derived>
	Derived transfer(const Eigen::MatrixBase<Derived>& generator, const typename group_t::rank_t& domain, const typename group_t::rank_t& range);

	///Restrict generator to level given the ranks at the original level (domain) and the target level (range).
	template<typename rank_t, typename Derived>
	Derived restriction(const Eigen::MatrixBase<Derived>& generator, const rank_t& domain, const rank_t& range);

	///Compute the act group action on a generator given the rank of the group it lives in.
	template<typename rank_t, typename Derived>
	Derived action(const Eigen::MatrixBase<Derived>& generator, const rank_t& rank);

	///The inverse of the restriction function on a generator given the ranks at the original level (domain) and the target level (range).
	///
	/// Recall that free Mackey functors have injective restrictions, so we only need the generator to be in the image.
	template<typename rank_t, typename T>
	T invRes(const T& generator, const rank_t& domain, const rank_t& range);

	///Writing the transfer of each generator in terms of the generators in the image.
	template<typename group_t>
	dense_t<typename group_t::rank_t> transfer(const Homology<typename group_t::rank_t, typename group_t::diff_t>& low, const Homology<typename group_t::rank_t, typename group_t::diff_t>& high, const typename group_t::rank_t& rank_low, const typename group_t::rank_t& rank_high);

	///Writing the restriction of each generator in terms of the generators in the image.
	template<typename group_t>
	dense_t<typename group_t::rank_t> restriction(const Homology<typename group_t::rank_t, typename group_t::diff_t>& high, const Homology<typename group_t::rank_t, typename group_t::diff_t>& low, const typename group_t::rank_t& rank_high, const typename group_t::rank_t& rank_low);

	///Writing the act group action on each generator in terms of the other generators.
	template<typename group_t>
	dense_t<typename group_t::rank_t> action(const Homology<typename group_t::rank_t, typename group_t::diff_t>& H, const typename group_t::rank_t& rank);
}
#include "impl/Levels.ipp"
