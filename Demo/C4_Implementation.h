#include <Mackey/Compute.h>
//Demo for a group specific implementation for C4

//First set the variables

const int GroupSpecific::Variables::prime = 2;
const int GroupSpecific::Variables::power = 2;
const int GroupSpecific::Variables::reps = 2;
const std::vector<int> GroupSpecific::Variables::sphere_dimensions = { 1,2 };

using namespace Mackey;


//Then write the standard chains. The function signature and name need not be changed, only the definition


template<typename rank_t, typename diff_t, typename deg_t>
Chains<rank_t, diff_t> myfunction(int k, const deg_t& sphere) {


	//A list of all the matrices appearing for convenience/small performance gain. Alternatively just define them directly in what follows

	const static std::array<diff_t, 10> diffList = \
	{ altmatrix<diff_t>(2, 2, { 1,-1 }), altmatrix<diff_t>(2, 4, { 1,-1 }), altmatrix<diff_t>(4, 4, { 1,0,0,-1 }), altmatrix<diff_t>(4, 4, { 1,-1 }), \
		altmatrix<diff_t>(1, 2, { 1 }), altmatrix<diff_t>(1, 4, { 1 }), altmatrix<diff_t>(2, 2, { 1 }), altmatrix<diff_t>(2, 4, { 1 }), altmatrix<diff_t>(4, 4, { 1 }), altmatrix<diff_t>(4, 4, { 1,0,0,1 }) };

	std::vector<rank_t> rank(k + 1);
	std::vector<diff_t> diff(k + 1);
	const auto n = sphere[0];
	const auto m = sphere[1];

	for (int i = 0; i <= k; i++) {
		if (!(i % 2)) //even
		{
			if (i == 0) {
				rank[0] = Eigen::MatrixBase<rank_t>::Constant(1,1,1); //constant 1
				continue;
			}
			else if (i <= n) {
				diff[i] = diffList[0];
			}
			else if (i == n + 1 && m) {
				diff[i] = diffList[1];
			}
			else if (i <= n + 2 * m && !(n % 2)) {
				diff[i] = diffList[2];
			}
			else if (i <= n + 2 * m) {
				diff[i] = diffList[3];
			}
		}
		else if (i == 1)
		{
			if (n != 0) {
				diff[i] = diffList[4];
			}
			else if (m != 0) {
				diff[i] = diffList[5];
			}
		}
		else //odd >=3
		{
			if (i <= n) {
				diff[i] = diffList[6];
			}
			else if (i == n + 1) {
				diff[i] = diffList[7];
			}
			else if (i <= n + 2 * m && !(n % 2)) {
				diff[i] = diffList[8];
			}
			else if (i <= n + 2 * m) {
				diff[i] = diffList[9];
			}
		}
		rank[i] = rank_t::Constant(1, 1, diff[i].cols()); //the number of columns of the differential
	}
	return Chains<rank_t,diff_t>(rank, diff);
}

//Finally set the standard chains. No need to change the syntax in what follows:

template<typename rank_t, typename diff_t, typename deg_t>
const typename GroupSpecific::Function<rank_t, diff_t, deg_t>::functype GroupSpecific::Function<rank_t, diff_t, deg_t>::PositiveChains = &myfunction;
