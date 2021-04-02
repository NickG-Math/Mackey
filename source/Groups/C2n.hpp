#pragma once
#include "Chain_Complexes/Chains.hpp"
#include "Utility/General.hpp"

///	@file
///	@brief	Defines the groups \f$C_2^n\f$

namespace mackey
{

	///Template for groups \f$C_{2^n}\f$
	template<int N, typename...>
	struct C2Power;

	///Groups \f$C_{2^n}\f$ given rank and differential types
	template<int N, typename _rank, typename _diff>
	struct C2Power<N,_rank,_diff> {
		typedef _rank rank_t;
		typedef _diff diff_t;
		static constexpr int prime = 2;
		static constexpr int power = N;
		static const std::array<int, N> sphere_dimensions;
		template <typename deg_t>
		static Chains<rank_t, diff_t> PositiveChains(const deg_t&, int);
	};

	///Groups \f$C_{2^n}\f$ given coefficient type
	template<int N, typename coeff>
	struct C2Power<N, coeff> : public C2Power<N, Eigen::Matrix<short,1,-1>, Eigen::SparseMatrix<coeff>> {
		typedef Eigen::Matrix<short, 1, -1> rank_t;
		typedef Eigen::SparseMatrix<coeff> diff_t;
	};

	///Group \f$C_2\f$
	template <typename ... Args>
	using C2 = C2Power<1, Args...>;

	///Group \f$C_4\f$
	template <typename ... Args>
	using C4 = C2Power<2, Args...>;

	///Group \f$C_8\f$
	template <typename ... Args>
	using C8= C2Power<3, Args...>;

	///Group \f$C_16\f$
	template <typename ... Args>
	using C16 = C2Power<4, Args...>;

}

#include "impl/C2n.ipp"
