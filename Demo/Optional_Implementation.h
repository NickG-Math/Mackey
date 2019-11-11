#ifdef MACKEY_NAMES //optional

#include "Mackey/Compute.h"
#include <array>

//Demo for a group specific implementation for C4. 


//first set the names
template <typename rank_t>
const std::vector<std::string>  GroupSpecificOptional::MackeyList<rank_t>::names = { "0", "Z/2", "Z/4", "overline Z/2", "Z" , "Z_-", "L",  "L_-", "Q",  "Q^sharp", "Z/2+overline Z/2", "Z/2+Z_-", "L^sharp", "p^*L", "p^*L_-", "Z_-^flat", "Z/2+L", "Z/2+L" };

//now set the Mackey functors.
using namespace Mackey;
template <typename rank_t>
const std::vector<MackeyFunctor<rank_t>> listMackeys() {
	std::vector<MackeyFunctor<rank_t>> Mackeys;
	Mackeys.resize(18);
	rank_t a(9);

	a.setZero();
	Mackeys[0] = MackeyFunctor<rank_t>(a, 3);

	a << 0, 0, 2, 0, 0, 0, 0, 0, 0;
	Mackeys[1] = MackeyFunctor<rank_t>(a, 3);

	a << 0, 2, 4, 0, 2, 0, 1, 0, 1;
	Mackeys[2] = MackeyFunctor<rank_t>(a, 3);

	a << 0, 2, 0, 0, 0, 0, 0, 0, 1;
	Mackeys[3] = MackeyFunctor<rank_t>(a, 3);

	a << 1, 1, 1, 2, 2, 1, 1, 1, 1;
	Mackeys[4] = MackeyFunctor<rank_t>(a, 3);

	a << 1, 1, 0, 2, 0, 1, 0, -1, -1;
	Mackeys[5] = MackeyFunctor<rank_t>(a, 3);

	a << 1, 1, 1, 1, 1, 2, 2, 1, 1;
	Mackeys[6] = MackeyFunctor<rank_t>(a, 3);

	a << 1, 1, 2, 1, 1, 2, 0, -1, -1;
	Mackeys[7] = MackeyFunctor<rank_t>(a, 3);

	a << 0, 2, 2, 0, 1, 0, 0, 0, 1;
	Mackeys[8] = MackeyFunctor<rank_t>(a, 3);

	a << 0, 2, 2, 0, 0, 0, 1, 0, 1;
	Mackeys[9] = MackeyFunctor<rank_t>(a, 3);

	a << 0, 2, 2, 0, 0, 0, 0, 0, 1;
	Mackeys[10] = MackeyFunctor<rank_t>(a, 3);

	a << 1, 1, 2, 2, 0, 1, 0, -1, -1;
	Mackeys[11] = MackeyFunctor<rank_t>(a, 3);

	a << 1, 1, 1, 1, 2, 2, 1, 1, 1;
	Mackeys[12] = MackeyFunctor<rank_t>(a, 3);

	a << 1, 1, 1, 2, 1, 1, 2, 1, 1;
	Mackeys[13] = MackeyFunctor<rank_t>(a, 3);

	a << 1, 1, 2, 2, 1, 1, 0, -1, -1;
	Mackeys[14] = MackeyFunctor<rank_t>(a, 3);

	a << 1, 1, 0, 1, 0, 2, 0, -1, -1;
	Mackeys[15] = MackeyFunctor<rank_t>(a, 3);

	a << 1, 1, 0, 1, 0, 2, 0, 1, 1;
	auto A = MackeyFunctor<rank_t>(a, 3);
	A.Groups[2].resize(2);
	A.Groups[2] << 2, 1;
	A.Tr[1].resize(1);
	A.Tr[1][0].resize(2);
	A.Tr[1][0] << 0, 1;
	A.Res[1].resize(2);
	A.Res[1][0].resize(1);
	A.Res[1][0] << 0;
	A.Res[1][1].resize(1);
	A.Res[1][1] << 2;
	Mackeys[16] = A;

	A.Tr[1][0](0) = 1; //the other version of Z/2+L
	Mackeys[17] = A;
	return Mackeys;
}


///Set the list of Mackey functors
template <typename rank_t>
const std::vector<Mackey::MackeyFunctor<rank_t>>  GroupSpecificOptional::MackeyList<rank_t>::Mackeys = listMackeys<rank_t>();
#endif
