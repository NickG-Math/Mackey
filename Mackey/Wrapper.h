#pragma once
#include "Chains.h"
#include "MackeyFunctor.h"

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

	///The function for the standard differentials.
	template<typename rank_t, typename diff_t, typename deg_t>
	class Function {
	public:
		typedef std::pair<rank_t, diff_t>(*functype)(int, const deg_t&);	///<The type of a function pointer to StandardDiff
		const static functype StandardDiff;		///<The function pointer to StandardDiff

	};
}

namespace Mackey {

	//Alias for shorter names
	const auto prime = GroupSpecific::Variables::prime; ///< If G=C_p^n then prime=p
	const auto power = GroupSpecific::Variables::power; ///< If G=C_p^n then power=n
	const auto reps = GroupSpecific::Variables::reps; ///< The number of representations of the group

	///Computes the i-th differential of the standard chains given GroupSpecific::Function::StandardDiff.
	template<typename rank_t, typename diff_t, typename deg_t>
	std::pair<rank_t, diff_t> StandardDiff(int i, const deg_t& sphere) {
		return GroupSpecific::Function<rank_t, diff_t, deg_t>::StandardDiff(i, sphere);
	}

	///Raise power to given exponent
	template<typename T>
	inline T intexp(const T exponent) {
		if (prime == 2) {
			return 1 << exponent;
		}
		else {
			return static_cast<T>(std::pow(prime, exponent)); //this should be improved...
		}
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
		else {
			return static_cast<T>(std::log(exponentiated) / std::log(prime)); //this should be improved...
		}
	}

	///Compute the dimension of a representation sphere
	template<typename deg_t>
	inline int dimension(const deg_t& sphere) {
		int sum = 0;
		for (int i = 0; i < reps; i++) {
			sum += abs(sphere[i]) * GroupSpecific::Variables::sphere_dimensions[i];
		}
		return sum;
	}

	///Reindexes the homological degree k so that it is always 0<=k<=dimension(sphere)
	template <typename T>
	inline int Reindex(int k, const T& sphere) {
		for (int i = 0; i < reps; i++) {
			if (sphere[i] < 0) {
				k -= sphere[i] * GroupSpecific::Variables::sphere_dimensions[i];
			}
		}
		return k;
	}

	///Inverse of \ref Reindex "Reindex"
	template <typename T>
	inline int invReindex(int k, const T& sphere) {
		for (int i = 0; i < reps; i++) {
			if (sphere[i] < 0) {
				k += sphere[i] * GroupSpecific::Variables::sphere_dimensions[i];
			}
		}
		return k;
	}

	///Returns the standard Chains of the given sphere up to index i.
	template<typename rank_t, typename diff_t, typename deg_t>
	Chains<rank_t, diff_t> StandardChains(int i, const deg_t& sphere) {
		std::vector<rank_t> rank;
		std::vector<diff_t> diff;
		rank.reserve(i + 1);
		diff.reserve(i + 1);
		for (int j = 0; j <= i; j++) {
			std::pair<rank_t, diff_t> B = StandardDiff<rank_t, diff_t>(j, sphere);
			rank.push_back(B.first);
			diff.push_back(B.second);
		}
		return Chains<rank_t, diff_t>(rank, diff);
	}

}


#ifdef MACKEY_NAMES
///Contains optional group specific variables for printing Mackey functors.
namespace GroupSpecificOptional {

	using namespace Mackey;
	//////////////////////////////
/// A list of Mackey functors and their names.

/// Only used if macro MACKEY_NAMES is defined (list is optional). The diff_t template is used to get the coefficient type
//////////////////////////////
	template<typename rank_t, typename diff_t>
	class MackeyList {
	public:
		static const std::vector<MackeyFunctor<rank_t>> Mackeys;///<The list of Mackey functors
		static const std::vector<std::string> names;///<The list of names
	};
}

namespace Mackey{

	///Identifies the Mackey Functor with its Latex notation using the MackeyList. Only used if macro MACKEY_NAMES is defined.
	template<typename rank_t, typename diff_t>
	std::string identify(MackeyFunctor<rank_t>& M)
	{
		M.normalize();
		GroupSpecificOptional::MackeyList<rank_t, diff_t> List;
		for (decltype(List.Mackeys.size()) i = 0; i < List.Mackeys.size(); i++) {
			if (M.compare(List.Mackeys[i])) {
				return List.names[i];
			}
		}
		return std::string();
	}
}
#endif


