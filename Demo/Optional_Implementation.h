#ifdef MACKEY_NAMES //optional
#include <Mackey/Compute.h>

//Demo for a group specific implementation for C4. 
//
//first set the names
template <typename rank_t, typename diff_t>
const std::vector<std::string> listnames() {
	if constexpr (Mackey::is_finite_cyclic<typename diff_t::Scalar>()) {
		if constexpr (diff_t::Scalar::order == 2) {
			return { "0", "Z/2", "Q",  "Q^sharp", "overline Z/2", "R" , "R_-", "L",  "Z/2+overline Z/2",  "L^sharp", "p^*L", "R_-^flat", "Z/2+R_-", "Z/2+Q", "Z/2+Q^{sharp}", "Z/2+L" , "Z/2+R", "Z/2+L" };
		}
	}
	else if (std::is_integral_v<typename diff_t::Scalar>)
		return { "0", "Z/2", "Z/4", "overline Z/2", "Z" , "Z_-", "L",  "L_-", "Q",  "Q^sharp",  "Z/2+overline Z/2", "Z/2+Z_-", "L^sharp", "p^*L", "p^*L_-", "Z_-^flat", "Z/2+L", "Z/2+L" };
}


//now set the Mackey functors.
using namespace Mackey;
template <typename rank_t, typename diff_t>
const std::vector<MackeyFunctor<rank_t>> listMackeys() {
	std::vector<MackeyFunctor<rank_t>> Mackeys;

	if constexpr (Mackey::is_finite_cyclic<typename diff_t::Scalar>()) {
		if constexpr (diff_t::Scalar::order == 2) {
			Mackeys.resize(18);
			rank_t a(9);

			a.setZero();
			Mackeys[0] = MackeyFunctor<rank_t>(a, 3);

			a << 0, 0, 1, 0, 0, 0, 0, 0, 0;// Z/2
			Mackeys[1] = MackeyFunctor<rank_t>(a, 3);

			a << 0, 1, 1, 0, 1, 0, 0, 0, 1; //Q
			Mackeys[2] = MackeyFunctor<rank_t>(a, 3);

			a << 0, 1, 1, 0, 0, 0, 1, 0, 1; //Q^sharp
			Mackeys[3] = MackeyFunctor<rank_t>(a, 3);

			a << 0, 1, 0, 0, 0, 0, 0, 0, 1; //overline Z/2
			Mackeys[4] = MackeyFunctor<rank_t>(a, 3);

			a << 1, 1, 1, 0, 0, 1, 1, 1, 1; //R
			Mackeys[5] = MackeyFunctor<rank_t>(a, 3);

			a << 1, 1, 0, 0, 0, 1, 0, 1, 1; //R_-
			Mackeys[6] = MackeyFunctor<rank_t>(a, 3);

			a << 1, 1, 1, 1, 1, 0, 0, 1, 1; //L
			Mackeys[7] = MackeyFunctor<rank_t>(a, 3);


			a << 0, 1, 1, 0, 0, 0, 0, 0, 1; //Z/2+overline Z/2
			Mackeys[8] = MackeyFunctor<rank_t>(a, 3);

			a << 1, 1, 1, 1, 0, 0, 1, 1, 1; //L^sharp
			Mackeys[9] = MackeyFunctor<rank_t>(a, 3);

			a << 1, 1, 1, 0, 1, 1, 0, 1, 1; //p^*L
			Mackeys[10] = MackeyFunctor<rank_t>(a, 3);

			a << 1, 1, 0, 1, 0, 0, 0, 1, 1; //R_-^flat
			Mackeys[11] = MackeyFunctor<rank_t>(a, 3);

			a << 1, 1, 1, 0, 0, 1, 0, 1, 1; //Z/2+R_{-}
			Mackeys[12] = MackeyFunctor<rank_t>(a, 3);

			a << 0, 1, 0, 0, 0, 0, 0, 0, 1; // Z/2+Q
			auto A = MackeyFunctor<rank_t>(a, 3);
			A.Groups[2].resize(2);
			A.Groups[2] << 1, 1;
			A.Tr[1].resize(1);
			A.Tr[1][0].resize(2);
			A.Tr[1][0] << 0, 1;
			A.Res[1].resize(2);
			A.Res[1][0].resize(1);
			A.Res[1][0] << 0;
			A.Res[1][1].resize(1);
			A.Res[1][1] << 0;
			Mackeys[13] = A;

			A = MackeyFunctor<rank_t>(a, 3); // Z/2+Q^{sharp}
			A.Groups[2].resize(2);
			A.Groups[2] << 1, 1;
			A.Tr[1].resize(1);
			A.Tr[1][0].resize(2);
			A.Tr[1][0] << 0, 0;
			A.Res[1].resize(2);
			A.Res[1][0].resize(1);
			A.Res[1][0] << 0;
			A.Res[1][1].resize(1);
			A.Res[1][1] << 1;
			Mackeys[14] = A;


			a << 1, 1, 0, 1, 0, 0, 0, 1, 1; // Z/2+L
			A = MackeyFunctor<rank_t>(a, 3); //
			A.Groups[2].resize(2);
			A.Groups[2] << 1, 1;
			A.Tr[1].resize(1);
			A.Tr[1][0].resize(2);
			A.Tr[1][0] << 1, 0;
			A.Res[1].resize(2);
			A.Res[1][0].resize(1);
			A.Res[1][0] << 0;
			A.Res[1][1].resize(1);
			A.Res[1][1] << 0;
			Mackeys[15] = A;

			a << 1, 1, 0, 0, 0, 1, 0, 1, 1; // Z/2+R
			A = MackeyFunctor<rank_t>(a, 3); //
			A.Groups[2].resize(2);
			A.Groups[2] << 1, 1;
			A.Tr[1].resize(1);
			A.Tr[1][0].resize(2);
			A.Tr[1][0] << 0, 0;
			A.Res[1].resize(2);
			A.Res[1][0].resize(1);
			A.Res[1][0] << 0;
			A.Res[1][1].resize(1);
			A.Res[1][1] << 1;
			Mackeys[16] = A;

			A = Mackeys[15];  //the other version of Z/2+L
			A.Tr[1][0](0) = 0;
			A.Tr[1][0](1) = 1;
			Mackeys[17] = A;
			return Mackeys;
		}
	}
	else if (std::is_integral_v<typename diff_t::Scalar>) {

		Mackeys.resize(18);
		rank_t a(9);

		a.setZero();
		Mackeys[0] = MackeyFunctor<rank_t>(a, 3);

		a << 0, 0, 2, 0, 0, 0, 0, 0, 0; // Z/2
		Mackeys[1] = MackeyFunctor<rank_t>(a, 3);

		a << 0, 2, 4, 0, 2, 0, 1, 0, 1; // Z/4
		Mackeys[2] = MackeyFunctor<rank_t>(a, 3);

		a << 0, 2, 0, 0, 0, 0, 0, 0, 1; // overline Z/2
		Mackeys[3] = MackeyFunctor<rank_t>(a, 3);

		a << 1, 1, 1, 2, 2, 1, 1, 1, 1; // Z
		Mackeys[4] = MackeyFunctor<rank_t>(a, 3);

		a << 1, 1, 0, 2, 0, 1, 0, -1, -1; // Z_-
		Mackeys[5] = MackeyFunctor<rank_t>(a, 3);

		a << 1, 1, 1, 1, 1, 2, 2, 1, 1; // L
		Mackeys[6] = MackeyFunctor<rank_t>(a, 3);

		a << 1, 1, 2, 1, 1, 2, 0, -1, -1; // L_-
		Mackeys[7] = MackeyFunctor<rank_t>(a, 3);

		a << 0, 2, 2, 0, 1, 0, 0, 0, 1; //Q
		Mackeys[8] = MackeyFunctor<rank_t>(a, 3);

		a << 0, 2, 2, 0, 0, 0, 1, 0, 1; //Q^{sharp}
		Mackeys[9] = MackeyFunctor<rank_t>(a, 3);

		a << 0, 2, 2, 0, 0, 0, 0, 0, 1; //Z/2+overline Z/2
		Mackeys[10] = MackeyFunctor<rank_t>(a, 3);

		a << 1, 1, 2, 2, 0, 1, 0, -1, -1; //Z/2+Z_-
		Mackeys[11] = MackeyFunctor<rank_t>(a, 3);

		a << 1, 1, 1, 1, 2, 2, 1, 1, 1; //L^sharp
		Mackeys[12] = MackeyFunctor<rank_t>(a, 3);

		a << 1, 1, 1, 2, 1, 1, 2, 1, 1; //p^*L
		Mackeys[13] = MackeyFunctor<rank_t>(a, 3);

		a << 1, 1, 2, 2, 1, 1, 0, -1, -1; // p^*L_-
		Mackeys[14] = MackeyFunctor<rank_t>(a, 3);

		a << 1, 1, 0, 1, 0, 2, 0, -1, -1; //Z_-^flat 
		Mackeys[15] = MackeyFunctor<rank_t>(a, 3);

		a << 1, 1, 0, 1, 0, 2, 0, 1, 1; //Z/2+L
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
	}
	return Mackeys;
}



template <typename rank_t, typename diff_t>
const std::vector<std::string>  GroupSpecificOptional::MackeyList<rank_t, diff_t>::names = listnames<rank_t, diff_t>();

template <typename rank_t, typename diff_t>
const std::vector<Mackey::MackeyFunctor<rank_t>>  GroupSpecificOptional::MackeyList<rank_t,diff_t>::Mackeys = listMackeys<rank_t,diff_t>();
#endif
