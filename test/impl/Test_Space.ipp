#pragma once
#include <iostream>
#include "Spaces/BC4S2.hpp"

namespace mackey
{
	namespace test
	{
		template <typename group_t>
		void test_BC4S2()
		{
			int j = 3;
			auto B = BC4S2<group_t>(j);
			auto r = B.ROCohomology({2, 0});
			std::cout << r;
		}
	}
}