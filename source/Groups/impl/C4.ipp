#include "../C4.hpp"

namespace mackey
{
	//PositiveChains for C_4
	template <typename rank_t, typename diff_t>
	template <typename deg_t>
	Chains<rank_t, diff_t> C_4<rank_t, diff_t>::PositiveChains(const deg_t &sphere, int k)
	{

		//A list of all the matrices appearing for convenience/small performance gain.

		const static std::array<diff_t, 10> diffList =
			{altmatrix<diff_t>(2, 2, {1, -1}), altmatrix<diff_t>(2, 4, {1, -1}), altmatrix<diff_t>(4, 4, {1, 0, 0, -1}), altmatrix<diff_t>(4, 4, {1, -1}),
			 altmatrix<diff_t>(1, 2, {1}), altmatrix<diff_t>(1, 4, {1}), altmatrix<diff_t>(2, 2, {1}), altmatrix<diff_t>(2, 4, {1}), altmatrix<diff_t>(4, 4, {1}), altmatrix<diff_t>(4, 4, {1, 0, 0, 1})};

		std::vector<rank_t> rank(k + 1);
		std::vector<diff_t> diff(k + 1);
		const auto n = sphere[0];
		const auto m = sphere[1];

		for (int i = 0; i <= k; i++)
		{
			if (!(i % 2)) //even
			{
				if (i == 0)
				{
					rank[0] = rank_t::Constant(1, 1, 1); //constant 1
					continue;
				}
				else if (i <= n)
					diff[i] = diffList[0];
				else if (i == n + 1 && m)
					diff[i] = diffList[1];
				else if (i <= n + 2 * m && !(n % 2))
					diff[i] = diffList[2];
				else if (i <= n + 2 * m)
				{
					diff[i] = diffList[3];
				}
			}
			else if (i == 1)
			{
				if (n != 0)
					diff[i] = diffList[4];
				else if (m != 0)
					diff[i] = diffList[5];
			}
			else //odd >=3
			{
				if (i <= n)
					diff[i] = diffList[6];
				else if (i == n + 1)
					diff[i] = diffList[7];
				else if (i <= n + 2 * m && !(n % 2))
					diff[i] = diffList[8];
				else if (i <= n + 2 * m)
					diff[i] = diffList[9];
			}
			rank[i] = rank_t::Constant(1, 1, diff[i].cols()); //the number of columns of the differential
		}
		return Chains<rank_t, diff_t>(rank, diff);
	}
}