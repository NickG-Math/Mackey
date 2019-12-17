#include <Mackey/Compute.h>
#include <Mackey/Polynomial.h>

//Demo for a group specific implementation for C_2^n. Currently for n=2, you only have to set the power variable below to the n you want

//First set the variables

const int GroupSpecific::Variables::prime = 2;
const int GroupSpecific::Variables::power = 2;
const int GroupSpecific::Variables::reps = GroupSpecific::Variables::power;

std::vector<int> power2_dimensions(int power) {
	std::vector<int> dims(power, 2);
	dims[0] = 1;
	return dims;
}

const std::vector<int> GroupSpecific::Variables::sphere_dimensions = power2_dimensions(GroupSpecific::Variables::power);


using namespace Mackey;

//The d and x sequences are defined in HHR17 pg 392.

template<typename deg_t>
std::vector<int> dsequence(int n, const deg_t& sphere) {
	std::vector<int> d(sphere.size() + 1);
	for (int i = 0; i < sphere.size(); i++) {
		d[i + 1] = d[i] + sphere[i] * GroupSpecific::Variables::sphere_dimensions[i];
	}
	return d;
}

template< typename Scalar, typename deg_t>
std::vector<Polynomial<Scalar>> xsequence(int n, const deg_t& sphere, const std::vector<int>& d) {
	std::vector<Polynomial<Scalar>> x(n +1);
	if (n <2)
		return x;
	x[2] = Polynomial<Scalar>(1);
	for (int i = 3; i <= n; i++) {
		int j = 1;
		while (2 + d[j - 1] >= i || i > 2 + d[j])
			j++;
		x[i] = geometric<Scalar>((1 << j) - 1) / x[i - 1];
	}
	return x;
}

//Then write the standard chains. The function signature and name need not be changed, only the definition
template<typename rank_t, typename diff_t, typename deg_t>
Chains<rank_t, diff_t> myfunction(int n, const deg_t& sphere) {

	auto d = dsequence(n,sphere);
	auto x = xsequence<typename diff_t::Scalar>(n,sphere, d);
	std::vector<rank_t> rank;
	std::vector<diff_t> diff;

	rank.reserve(n + 1);
	rank.push_back(altmatrix<rank_t>(1, 1, { 1 }));
	for (int i = 1; i <= n; i++) {
		int j = 1;
		while (d[j - 1] >= i || i > d[j])
			j++;
		rank.push_back(altmatrix<rank_t>(1, 1, { static_cast<typename rank_t::Scalar>(1 << j) }));
	}
	diff.resize(n + 1);
	if (rank.size() == 1) {
		return Chains<rank_t, diff_t>(rank, diff);
	}
	else {
		diff[1]=altmatrix<diff_t>(rank[0][0], rank[1][0], { 1 });
		std::vector<typename diff_t::Scalar> v = { 1,-1 };
		Polynomial<typename diff_t::Scalar> r(v); //1-g polynomial
		for (int i = 2; i <= n; i++) {
			auto relation = Polynomial<typename diff_t::Scalar>(1) - Polynomial<typename diff_t::Scalar>(1, rank[i - 1][0]);
			if (i % 2) {//odd
				auto s = x[i] % relation; 
				s.pad(rank[i - 1][0]);
				diff[i]=altmatrix<diff_t>(rank[i - 1][0], rank[i][0], s.p);
			}
			else { //even
				auto rxd = (r * x[i]) % relation;
				rxd.pad(rank[i - 1][0]);
				diff[i]=altmatrix<diff_t>(rank[i - 1][0], rank[i][0], rxd.p);
			}
		}
		return Chains<rank_t, diff_t>(rank, diff);
	}
}

//Finally set the standard chains. No need to change the syntax in what follows:

template<typename rank_t, typename diff_t, typename deg_t>
const typename GroupSpecific::Function<rank_t, diff_t, deg_t>::functype GroupSpecific::Function<rank_t, diff_t, deg_t>::PositiveChains = &myfunction;
