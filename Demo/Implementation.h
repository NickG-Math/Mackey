#include <Mackey/Compute.h>
//Demo for a group specific implementation for C4

//First set the variables

const int GroupSpecific::Variables::prime = 2;
const int GroupSpecific::Variables::power = 2;
const int GroupSpecific::Variables::reps = 2;
const std::vector<int> GroupSpecific::Variables::sphere_dimensions = { 1,2 };

//Then write the standard chains. The function signature and name need not be changed, only the definition

template<typename rank_t, typename diff_t, typename deg_t>
std::pair<rank_t, diff_t> myfunction(int i, const deg_t& sphere) {

	using namespace Mackey;

	//A list of all the matrices appearing for convenience/small performance gain. Alternatively just define them directly in what follows
	const static std::array<diff_t, 10> diffList = \
	{ altmatrix<diff_t>(2, 2, { 1,-1 }), altmatrix<diff_t>(2, 4, { 1,-1 }), altmatrix<diff_t>(4, 4, { 1,0,0,-1 }), altmatrix<diff_t>(4, 4, { 1,-1 }), \
		altmatrix<diff_t>(1, 2, { 1 }), altmatrix<diff_t>(1, 4, { 1 }), altmatrix<diff_t>(2, 2, { 1 }), altmatrix<diff_t>(2, 4, { 1 }), altmatrix<diff_t>(4, 4, { 1 }), altmatrix<diff_t>(4, 4, { 1,0,0,1 }) };

	rank_t rank(1);
	diff_t diff;
	const auto n = sphere[0];
	const auto m = sphere[1];
	if (n > 0 || m > 0)
	{
		if (!(i % 2)) //even
		{
			if (i == 0) {
				rank << 1;
				return std::make_pair(rank, diff);
			}
			else if (i <= n) {
				diff = diffList[0];
			}
			else if (i == n + 1 && m) {
				diff = diffList[1];
			}
			else if (i <= n + 2 * m && !(n % 2)) {
				diff = diffList[2];
			}
			else if (i <= n + 2 * m) {
				diff = diffList[3];
			}
		}
		else if (i == 1)
		{
			if (n != 0) {
				diff = diffList[4];
			}
			else if (m != 0) {
				diff = diffList[5];
			}
		}
		else //odd >=3
		{
			if (i <= n) {
				diff = diffList[6];
			}
			else if (i == n + 1) {
				diff = diffList[7];
			}
			else if (i <= n + 2 * m && !(n % 2)) {
				diff = diffList[8];
			}
			else if (i <= n + 2 * m) {
				diff = diffList[9];
			}
		}
	}
	else if (n == 0 && m == 0)
	{
		rank << 1;
		return std::make_pair(rank, diff);
	}
	else
	{
		if (i == 0) {
			if (m != 0) {
				rank << 4;
				return std::make_pair(rank, diff);
			}
			else {
				rank << 2;
				return std::make_pair(rank, diff);
			}
		}
		else {
			std::vector<int> dualsphere = { -n,-m };
			std::pair<rank_t, diff_t> B = StandardDiff<rank_t, diff_t>(-n - 2 * m - i + 1, dualsphere);
			diff = B.second.transpose();
		}
	}
	rank << diff.cols();
	return std::make_pair(rank, diff);
}

//Finally set the standard chains. No need to change the syntax in what follows:

template<typename rank_t, typename diff_t, typename deg_t>
const typename GroupSpecific::Function<rank_t, diff_t, deg_t>::functype GroupSpecific::Function<rank_t, diff_t, deg_t>::StandardDiff=&myfunction;
