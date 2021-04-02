#include "../BC4S2.hpp"

namespace mackey
{

	namespace implementation_details
	{

		///	@param	v	Vector v
		///	@return 	v if its first nonzero index is >0 and -v if it's <0
		std::vector<char> first_pos(const std::vector<char> &v)
		{
			bool flag = 0;
			for (auto i : v)
				if (i != 0)
				{
					flag = (i < 0);
					break;
				}
			if (flag)
				return mackey::operator-(v);
			else
				return v;
		}

		///	@brief	BC4S2 helper functions
		///	@return	The binomial $\binom{n,m}$
		long binomial(int top, int bot)
		{
			if (bot > top - bot)
				bot = top - bot;
			if (bot == 0)
				return 1;
			int binom = top;
			for (int i = 1; i <= bot - 1; i++)
			{ // $n/1 * (n-1)/2 * \cdots (n-k+1)/k$
				binom *= (top - i);
				binom /= i + 1;
			}
			return binom;
		}

		///	@brief			Generates all combinations of -1, 0, 1
		///	@param total	The total number of -1,0,1
		///	@param nonzeros	The number of -1,1
		std::vector<std::vector<char>> all_plus_0_minus_combinations(int total, int nonzeros)
		{
			std::vector<std::vector<char>> combs;
			combs.reserve(2 * binomial(total, nonzeros));
			std::vector<char> min(total, -1);
			std::vector<char> max(total, 1);
			InterpolatingVectorGenerator<std::vector<char>> V(min, max);
			for (const auto &i : V)
			{
				int count = 0;
				for (const auto &e : i)
					if (e != 0)
					{
						count++;
						if (count > nonzeros)
							break;
					}
				if (count == nonzeros)
					combs.push_back(i);
			}
			return combs;
		}

		//C4 action on a BC4S2 coordinate
		template <typename T>
		T C4action(const T &x)
		{
			auto gx = x;
			for (size_t i = 0; i < x.size(); i++)
			{
				if (i % 4 == 1)
					gx[i] = -x[i]; //antipode
				else if (i % 4 == 2)
				{ //rotation
					gx[i] = -x[i + 1];
					gx[i + 1] = x[i];
				}
			}
			return gx;
		}

		//C4 robits of a BC4S2 coordinate
		template <typename T>
		std::vector<T> C4orbits(const T &x)
		{
			std::vector<T> orbits;
			orbits.push_back(x);
			for (int i = 1; i < 4; i++)
			{
				auto gx = C4action(orbits[i - 1]);
				if (x != gx && x != mackey::operator-(gx))
					orbits.push_back(gx);
				else
					break;
			}
			return orbits;
		}
	}

	template <typename group_t>
	auto BC4S2<group_t>::boundary(int i, int k)
	{
		boundary_t bd;
		bd.reserve(all_cells_in_dimension[i][k].size());
		for (int r = 0; r < all_cells_in_dimension[i][k].size(); r++)
		{
			if (all_cells_in_dimension[i][k][r] != 0)
			{
				auto newelement = all_cells_in_dimension[i][k];
				newelement[r] = 0;
				auto loc = cell_to_location[i - 1][implementation_details::first_pos(newelement)];
				bd.emplace_back(loc, 1);
			}
		}
		return bd;
	}

	//Implementing the estimate map!!
	template <typename group_t>
	long BC4S2<group_t>::estimate_nonzero_entries(int i)
	{
		return i * summation(this->cells[i]);
	}

	template <typename group_t>
	BC4S2<group_t>::BC4S2(int j)
	{
		this->cells.reserve(4 * j);
		all_cells_in_dimension.resize(4 * j);
		cell_to_location.resize(4 * j);
		for (int i = 1; i <= 4 * j; i++)
		{
			std::vector<scalar_t<rank_t>> cells_in_dimension_i;
			auto comb = implementation_details::all_plus_0_minus_combinations(4 * j, i);
			all_cells_in_dimension[i - 1].reserve(comb.size());
			cells_in_dimension_i.reserve(comb.size());
			int count = 0;
			for (const auto &q : comb)
			{
				if (cell_to_location[i - 1].find(implementation_details::first_pos(q)) == cell_to_location[i - 1].end())
				{
					auto oq = implementation_details::C4orbits(q);
					cells_in_dimension_i.push_back(oq.size());
					for (const auto &o : oq)
					{
						all_cells_in_dimension[i - 1].push_back(o);
						cell_to_location[i - 1][implementation_details::first_pos(o)] = count;
						count++;
					}
				}
			}
			this->cells.push_back(Eigen::Map<rank_t>(cells_in_dimension_i.data(), cells_in_dimension_i.size()));
		}
	}
}
