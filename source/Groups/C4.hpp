#pragma once
#include "Chain_Complexes/Chains.hpp"
#include "Utility/General.hpp"

///	@file
///	@brief		Defines group \f$C_4\f$
///	@details	Implementation is simpler compared to all \f$C_{2^n}\f$

//First set the variables

namespace mackey
{

	///Quick and easy implementation of group \f$C_4\f$
	template <typename _rank, typename _diff>
	struct C_4
	{
		typedef _rank rank_t;
		typedef _diff diff_t;
		static const int prime = 2;
		static const int power = 2;
		constexpr static std::array<int, 2> sphere_dimensions = {1, 2};
		template <typename deg_t>
		static Chains<rank_t, diff_t> PositiveChains(const deg_t &, int);
	};
}
#include "impl/C4.ipp"
