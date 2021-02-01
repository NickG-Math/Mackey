#pragma once
#include "Smith_Base.h"
///@file
///@brief Contains the Smith normal form algorithm for dense matrices using the minimum as pivot

namespace Mackey{

	///The dense Smith Normal Form using the minimum (in absolute value) nonzero element as pivot
	template <typename S_t, typename R_t, typename C_t>
	class SmithMP : Smith<S_t, R_t, C_t> {

		template<typename sS_t, typename sR_t, typename sC_t>
		friend Smith<sS_t, sR_t, sC_t> diagonalize(const sS_t&, bool, bool, bool);

		SmithMP(const S_t& A, bool wantP, bool wantQ, bool sort) : Smith<S_t, R_t, C_t>(A, wantP, wantQ, sort) {
			SmithIt();
			if (sort)
				this->sorter();
		}
		//so we don't have to write this-> all the time
		using Smith<S_t, R_t, C_t>::S;
		using Smith<S_t, R_t, C_t>::P;
		using Smith<S_t, R_t, C_t>::Q;
		using Smith<S_t, R_t, C_t>::Pi;
		using Smith<S_t, R_t, C_t>::Qi;
		using Smith<S_t, R_t, C_t>::M;
		using Smith<S_t, R_t, C_t>::N;
		using Smith<S_t, R_t, C_t>::diagonal;
		using Smith<S_t, R_t, C_t>::wantP;
		using Smith<S_t, R_t, C_t>::wantQ;

		void SmithIt();
		void workRow(int);
		void workCol(int);
		void pivoting(int);
		void find_pivot(int, int&, int&);
	};

	template <typename S_t, typename R_t, typename C_t>
	void SmithMP<S_t, R_t, C_t>::SmithIt() {
		bool rowZero, columnZero;
		for (int start = 0; start < std::min(M, N); start++) {
			rowZero = S.row(start).tail(N - start - 1).isZero();
			do {
				if (!rowZero) {
					pivoting(start);
					workRow(start);
					rowZero = S.row(start).tail(N - start - 1).isZero();
				}
				columnZero = S.col(start).tail(M - start - 1).isZero();
				if (!columnZero) {
					pivoting(start);
					workCol(start);
					rowZero = S.row(start).tail(N - start - 1).isZero();
					columnZero = S.col(start).tail(M - start - 1).isZero();
				}
			} while (!rowZero || !columnZero);
			this->diagonal[start] = S(start, start);
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithMP<S_t, R_t, C_t>::workRow(int start) {
		for (int j = start + 1; j < N; j++) {
			if (S(start, j) != 0) {
				auto thequotient = floor_division(S(start, j), S(start,start));
				S.col(j).tail(M - start) -= thequotient * S.col(start).tail(M - start);
				if (wantQ) {
					Q.col(j) -= static_cast<Scalar_t<C_t>>(thequotient)* Q.col(start);
					Qi.row(start) += static_cast<Scalar_t<R_t>>(thequotient)* Qi.row(j);
				}
			}
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithMP<S_t, R_t, C_t>::workCol(int start) {
		for (int i = start + 1; i < M; i++) {
			if (S(i, start) != 0) {
				auto thequotient = floor_division(S(i, start), S(start, start));
				S.row(i).tail(N - start) -= thequotient * S.row(start).tail(N - start);
				if (wantP) {
					P.row(i) -= static_cast<Scalar_t<R_t>>(thequotient)* P.row(start);
					Pi.col(start) += static_cast<Scalar_t<C_t>>(thequotient)* Pi.col(i);
				}
			}
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithMP<S_t, R_t, C_t>::pivoting(int start)
	{
		if (abs(S(start, start)) == 1)
			return;
		int s, t;
		s = t = start;
		find_pivot(start, s, t);
		//bring pivot to the start
		if (t != start) {
			S.col(start).tail(M - start).swap(S.col(t).tail(M - start));
			if (wantQ) {
				Q.col(start).swap(Q.col(t));
				Qi.row(start).swap(Qi.row(t));
			}
		}
		if (s != start) {
			S.row(start).tail(N - start).swap(S.row(s).tail(N - start));
			if (wantP) {
				P.row(start).swap(P.row(s));
				Pi.col(start).swap(Pi.col(s));
			}
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithMP<S_t, R_t, C_t>::find_pivot(int start, int& s, int& t)
	{
		Scalar_t<S_t> min = 0;
		for (int j = start; j < S.cols();j++) {
			for (int i = start; i < S.rows(); i++) {
				if (S(i, j) != 0 && (min == 0 || abs(S(i, j)) < min)) {
					min = abs(S(i, j));
					s = i;
					t = j;
					if (min == 1)
						return;
				}
			}
		}
	}

}
