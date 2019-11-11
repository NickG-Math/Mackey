#pragma once
#include <string>
#include <vector>
#include <Eigen/Dense>

///@file
///@brief Contains the class describing Mackey Functors.




namespace {

	template<typename Derived>
	inline bool genericcompare(const Eigen::MatrixBase<Derived>& A, const Eigen::MatrixBase<Derived>& B) {
		if (A.size() != B.size() || A != B) {
			return 0;
		}
		return 1;
	}

	template<typename T>
	inline bool genericcompare(const std::vector<T>& A, const std::vector<T>& B) {
		if (A.size() != B.size()) {
			return 0;
		}
		else {
			for (size_t i = 0; i < A.size(); i++) {
				if (!genericcompare(A[i], B[i])) {
					return 0;
				}
			}
		}
		return 1;
	}
}

namespace Mackey
{

	///A Mackey Functor
	template<typename rank_t>
	class MackeyFunctor {
	public:
		std::vector<rank_t> Groups;///<The Groups starting at level 0 (see Homology for the labeling)

		///The transfers starting at level 0. It's encoded as a vectors of vectors of arrays because the answer may be noncyclic
		///Example (C4): If Tr[0]={[2]} and Tr[1]={[1,1],[2,3]} then Tr_0^2(gen)=2*gen and Tr_2^4(gen0)=gen0+gen1 and Tr_2^4(gen1)=2*gen0+3*gen1
		std::vector<std::vector<rank_t>> Tr;

		std::vector<std::vector<rank_t>> Res;///<Same as transfers, but now using restrictions
		std::vector<std::vector<rank_t>> Weyl;///<Same as transfers, but now using the Weyl group action

		void resize(int levels) {
			Groups.resize(levels); Tr.resize(levels - 1); Res.resize(levels - 1); Weyl.resize(levels - 1);
		}
		
		MackeyFunctor() {};///<Default constructor		
		MackeyFunctor(const rank_t&, int);///<When there is only one generator at each level we can use this convenient constructor

		
		bool compare(const MackeyFunctor<rank_t>&);///<Compares two Mackey Functors (after normalizing) and returns 1 if equal.
		void normalize();///<Normalizes a Mackey Functor.
	};



	template<typename rank_t>
	inline bool MackeyFunctor<rank_t>::compare(const MackeyFunctor<rank_t>& M)
	{
		if (!genericcompare(Groups, M.Groups) || !genericcompare(Tr, M.Tr) || !genericcompare(Res, M.Res) || !genericcompare(Weyl, M.Weyl)) {
			return 0;
		}
		return 1;
	}


	template<typename rank_t>
	void MackeyFunctor<rank_t>::normalize()
	{
		for (size_t i = 0; i < Groups.size() - 1; i++) {
			if (Groups[i].cols() == 1 && Groups[i + 1].cols() == 1) {//both cyclic
				if (Tr[i][0](0) < 0 && Res[i][0](0) < 0) {
					Tr[i][0](0) = -Tr[i][0](0);
					Res[i][0](0) = -Res[i][0](0);
				}
			}
		}
	}



	template<typename rank_t>
	MackeyFunctor<rank_t>::MackeyFunctor(const rank_t& bigarray, int length)
	{
		resize(length);
		for (int i = 0; i < length; i++) {
			if (bigarray(i) != 0) {
				Groups[i] = bigarray.segment(i, 1);
			}
		}
		for (int i = 0; i < length - 1; i++) {
			if (bigarray(i) != 0 && bigarray(i + 1) != 0) {
				Tr[i].push_back(bigarray.segment(i + length, 1));
				Res[i].push_back(bigarray.segment(i + 2 * length - 1, 1));
			}
			if (bigarray(i) != 0) {
				Weyl[i].push_back(bigarray.segment(i + 3 * length - 2, 1));
			}
		}
	}

}
