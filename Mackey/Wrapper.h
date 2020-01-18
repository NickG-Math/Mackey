#pragma once
#include "Chains.h"
#include "MackeyFunctor.h"
#include "General.h"

#if !defined __INTEL_COMPILER && defined _MSC_VER
#include <intrin.h>
#endif 

///@file
///@brief Wraps group-specific implementations into a general interface.

///Contains everything that's group specific and needs to be set manually by the user.
namespace GroupSpecific {
	///The global variables
	class Variables {
	public:
		const static int prime; ///< If G=C_p^n then prime=p
		const static int power; ///< If G=C_p^n then power=n
		const static int reps; ///< The number of representations of the group
		const static std::vector<int> sphere_dimensions; ///< The dimensions of the representations of the group (the order has to be fixed beforehand)
	};

	///The function for the Chains up to a given index corresponding to non-virtual representations.
	template<typename rank_t, typename diff_t, typename deg_t>
	class Function {
	public:
		typedef Mackey::Chains<rank_t, diff_t>(*functype)(int, const deg_t&);	///<The type of a function pointer to PositiveChains
		const static functype PositiveChains;		///<The function pointer to PositiveChains

	};
}

namespace Mackey {

	//Alias for shorter names
	const auto prime = GroupSpecific::Variables::prime; ///< If G=C_p^n then prime=p
	const auto power = GroupSpecific::Variables::power; ///< If G=C_p^n then power=n
	const auto reps = GroupSpecific::Variables::reps; ///< The number of representations of the group

	///Raise power to given exponent
	template<typename T>
	inline T intexp(const T exponent) {
		if (prime == 2)
			return 1 << exponent;
		else
			return static_cast<T>(std::pow(prime, exponent)); //this should be improved...
	}

	///Compute log(exponentiated) with base power
	template<typename T>
	inline T intlog(const T exponentiated) {
		if (prime == 2) { //count the trailing bits
#ifdef __INTEL_COMPILER
			return _bit_scan_forward(exponentiated);
#elif _MSC_VER
			unsigned long index;
			_BitScanReverse64(&index,exponentiated);
			return index;
#else //__GNUC__ || _clang_
			return __builtin_ctz(exponentiated);
#endif
		}
		else
			return static_cast<T>(std::log(exponentiated) / std::log(prime)); //this should be improved...
	}

	///Compute the dimension of a representation sphere
	template<typename deg_t>
	inline int dimension(const deg_t& sphere) {
		int sum = 0;
		for (int i = 0; i < reps; i++)
			sum += abs(sphere[i]) * GroupSpecific::Variables::sphere_dimensions[i];
		return sum;
	}

	///Reindexes the homological degree k so that it is always 0<=k<=dimension(sphere)
	template <typename T>
	inline int Reindex(int k, const T& sphere) {
		for (int i = 0; i < reps; i++) {
			if (sphere[i] < 0)
				k -= sphere[i] * GroupSpecific::Variables::sphere_dimensions[i];
		}
		return k;
	}

	///Inverse of \ref Reindex "Reindex"
	template <typename T>
	inline int invReindex(int k, const T& sphere) {
		for (int i = 0; i < reps; i++) {
			if (sphere[i] < 0)
				k += sphere[i] * GroupSpecific::Variables::sphere_dimensions[i];
		}
		return k;
	}

	///Returns the standard Chains of the given sphere up to index i.
	template<typename rank_t, typename diff_t, typename deg_t>
	Chains<rank_t, diff_t> StandardChains(int i, const deg_t& sphere) {
		bool cohomology = 0;
		for (const auto& i : sphere) {
			if (i < 0) {
				cohomology = 1;
				break;
			}
		}
		if (!cohomology)
			return GroupSpecific::Function<rank_t, diff_t, deg_t>::PositiveChains(i, sphere);
		else {
			auto C = GroupSpecific::Function<rank_t, diff_t, deg_t>::PositiveChains(dimension(-sphere), -sphere);
			return C.dualize(i);
		}
	}


	int computeweight(const std::vector<int>& deg) {
		int boxnumber = 0;
		for (int i = 1; i < deg.size() - 1; i++)
			if (deg[i] * deg[i + 1] < 0)
				boxnumber++;
		int weight=0;
		for (int i = 1; i < deg.size(); i++)
			weight += std::pow(abs(deg[i]), GroupSpecific::Variables::sphere_dimensions[i - 1]);
		weight *= (abs(deg[0]) + 1);
		return std::pow(weight, boxnumber+1);
	}
}