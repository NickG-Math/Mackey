#pragma once
#include <vector>
#include <cmath>
#include "Aliases.h"
#include "SFINAE.h"

///@file
///@brief Contains the base class for all Smith implementations. 

namespace {

	template<typename Scalar>
	inline Scalar floor_division(Scalar x, Scalar y) {
		if constexpr (std::is_floating_point_v<Scalar>)
			return floor(x / y);
		else if (std::is_integral_v<Scalar>) { //C++ rounds integer division as opposed to floor it!
			auto q = x / y;
			auto r = x % y;
			if (r != 0)
				return q - ((x < 0) ^ (y < 0));
			return q;
		}
		else
			return x;
	}

	///Produces the transpositions that sort a given vector of pairs
	template<typename T>
	std::vector<int> transpositions(std::vector<std::pair<T, int>> a) {
		std::vector<int> transpositions;
		if (a.size() <= 1)
			return transpositions;
		transpositions.reserve((a.size() - 1) * (a.size() - 1));
		for (int j = 0; j < a.size() - 1; j++) {
			for (int i = 0; i < a.size() - j - 1; i++) {
				if (a[i].first < a[i + 1].first) {
					transpositions.push_back(i);
					auto temp = a[i];
					a[i] = a[i + 1];
					a[i + 1] = temp;
				}
			}
		}
		return transpositions;
	}


	template<typename T>
	void swapCol(Eigen::SparseMatrix<T, 0>& A, int i, int j) {
		Eigen::SparseMatrix<T, 0> Ai = A.col(i);
		Eigen::SparseMatrix<T, 0> Aj = A.col(j);
		A.col(i) = Aj;
		A.col(j) = Ai;
	}

	template<typename T>
	void swapRow(Eigen::SparseMatrix<T, 1>& A, int i, int j) {
		Eigen::SparseMatrix<T, 1> Ai = A.row(i);
		Eigen::SparseMatrix<T, 1> Aj = A.row(j);
		A.row(i) = Aj;
		A.row(j) = Ai;
	}

	template<typename T>
	void swapCol(Eigen::MatrixBase<T>& A, int i, int j) {
		A.col(i).swap(A.col(j));
	}

	template<typename T>
	void swapRow(Eigen::MatrixBase<T>& A, int i, int j) {
		A.row(i).swap(A.row(j));
	}

}

namespace Mackey {

		template <typename sS_t, typename sR_t, typename sC_t>
		class SmithMP;
		template <typename sS_t, typename sR_t, typename sC_t>
		class SmithFP;
		template <typename sS_t, typename sR_t, typename sC_t>
		class SmithSP;
		template <typename sS_t, typename sR_t, typename sC_t>
		class SmithSparse;



	///////////////////////////////////////////////////////////////
	/// The Smith Normal Form and its coefficient matrices

	/// The SNF of A is a diagonal matrix S with invertible coefficient matrices P,Q,Pi,Qi s.t. P*A*Q=S and Pi*S*Qi=A
	///
	/// P,Pi are inverses and so are Q,Qi and they are all computed through row-column elimination
	///
	///S_t is the type of the original matrix.
	///R_t is S_t but row major for best performance (required for sprase)
	///C_t is S_t but column major for best performance (required for sprase)
	/////////////////////////////////////////////////////////////
	template <typename S_t, typename R_t, typename C_t>
	class Smith {
	public:
		row_t<S_t> diagonal;			///< The diagonal of the Smith normal form
		R_t P;					///< One of the coefficient matrices (S=P*A*Q)
		R_t Qi;					///< One of the coefficient matrices (A=Pi*S*Qi)
		C_t Q;					///< One of the coefficient matrices (S=P*A*Q)
		C_t Pi;					///< One of the coefficient matrices (A=Pi*S*Qi)

	protected:
		///Initializes the member variables. Actual computation is done by the implementation classes
		Smith(const S_t&, bool, bool, bool);

		S_t S;					///< The Smith Normal Form as a diagonal matrix
		const int M; 				///< Row dimension of original matrix
		const int N;				///< Column dimension of original matrix
		const bool wantP;			///< Whether we want to compute the P and Pi coefficient matrices
		const bool wantQ;			///< Whether we want to compute the Q and Qi coefficient matrices
		void sorter();				///< Sorts the Smith Normal Form to have decreasing entries in absolute value (except +1,-1)

	};

	template<typename S_t, typename R_t, typename C_t>
	Smith<S_t, R_t, C_t> ::Smith(const S_t& A, bool wantP, bool wantQ, bool sort)
		: S(A), M(S.rows()), N(S.cols()), wantP(wantP), wantQ(wantQ) {
		if (wantP) {
			P.resize(M, M);
			Pi.resize(M, M);
			P.setIdentity();
			Pi.setIdentity();
		}
		if (wantQ) {
			Q.resize(N, N);
			Qi.resize(N, N);
			Q.setIdentity();
			Qi.setIdentity();
		}
		diagonal.resize(std::min(M, N));
	}


	template<typename T>
	void swap(T& a, T& b) {
		T temp = a;
		a = b;
		b = temp;
	}

	template<typename S_t, typename R_t, typename C_t>
	void Smith<S_t, R_t, C_t> ::sorter() {
		std::vector<std::pair<Scalar_t<S_t>, int>> toBeSorted;
		toBeSorted.reserve(diagonal.size());
		for (int i = 0; i < diagonal.size(); i++) {
			if (abs(diagonal[i]) != 1)
				toBeSorted.push_back(std::make_pair(abs(diagonal[i]), i));
		}
		auto transp = transpositions(toBeSorted);
		for (const auto& i : transp) {
			auto first = toBeSorted[i].second;
			auto second = toBeSorted[i + 1].second;
			swap(diagonal[first], diagonal[second]);
			if (wantP) {
				swapRow(P,first,second);
				swapCol(Pi, first, second);
			}
			if (wantQ) {
				swapCol(Q, first, second);
				swapRow(Qi, first, second);
			}
		}
	}

	///Returns the Smith Normal Form of a matrix, with the coefficient matrices and sorted if desired
	template <typename S_t, typename R_t, typename C_t>
	Smith<S_t, R_t, C_t> diagonalize(const S_t& A, bool wantP, bool wantQ, bool sort) {
		if constexpr (is_Sparse<S_t>::value)
			return SmithSparse<S_t, R_t, C_t>(A, wantP, wantQ, sort);
		else {
			if constexpr (is_Finite_Cyclic<Scalar_t<S_t>>::value)
				return SmithFP<S_t, R_t, C_t>(A, wantP, wantQ, sort);
			else {
				if (A.rows() >0 && A.cols() >0)
					return SmithSP<S_t, R_t, C_t>(A, wantP, wantQ, sort);
				else
					return SmithMP<S_t, R_t, C_t>(A, wantP, wantQ, sort);
			}
		}
	}

	///Returns the unsorted Smith Normal Form of a matrix, with the coefficient matrices if desired
	template <typename S_t, typename R_t, typename C_t>
	Smith<S_t, R_t, C_t> diagonalize(const S_t& A, bool wantP, bool wantQ) {
		return diagonalize<S_t, R_t, C_t>(A, wantP, wantQ, 0);
	}

}
