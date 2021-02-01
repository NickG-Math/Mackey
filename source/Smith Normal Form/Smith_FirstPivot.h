#pragma once
#include "Smith_Base.h"
///@file
///@brief Contains the Smith normal form algorithm for dense matrices using the first nonzero as pivot

namespace Mackey{

	///The dense Smith Normal Form using the first nonzero element as pivot
	template <typename S_t, typename R_t, typename C_t>
	class SmithFP : Smith<S_t, R_t, C_t> {

		template<typename sS_t, typename sR_t, typename sC_t>
		friend Smith<sS_t, sR_t, sC_t> diagonalize(const sS_t&, bool, bool, bool);

		SmithFP(const S_t& A, bool wantP, bool wantQ, bool sort) : Smith<S_t, R_t, C_t>(A, wantP, wantQ, sort) {
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
		using Smith<S_t, R_t, C_t>::diagonal;
		using Smith<S_t, R_t, C_t>::M;
		using Smith<S_t, R_t, C_t>::N;
		using Smith<S_t, R_t, C_t>::wantP;
		using Smith<S_t, R_t, C_t>::wantQ;

		void SmithIt();
		void workRow(int, int&, int&);
		void searchRow(int, int&, int&);
		void workCol(int, int&, int&);
		void searchCol(int, int&, int&);

	};


	template <typename S_t, typename R_t, typename C_t>
	void SmithFP<S_t, R_t, C_t>::SmithIt() {
		for (int start = 0; start < std::min(M, N); start++) {
			auto i = start + 1;
			auto j = start + 1;

			while (i < M || j < N) {
				searchCol(start, i, j);
				while (i < M) {
					workCol(start, i, j);
					searchCol(start, i, j);
				}
				searchRow(start, i, j);
				while (j < N) {
					workRow(start, i, j);
					searchRow(start, i, j);
				}
			}
			this->diagonal[start] = S(start, start);
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithFP<S_t, R_t, C_t>::workRow(int start, int& i, int& j) {
		if (S(start, j) % S(start, start) == 0) {
			auto thequotient = S(start, j) / S(start, start);
			S.col(j).tail(M - start) += -thequotient * S.col(start).tail(M - start);
			if (wantQ) {
				Q.col(j) += -static_cast<Scalar_t<C_t>>(thequotient)* Q.col(start);
				Qi.row(start) += static_cast<Scalar_t<R_t>>(thequotient)* Qi.row(j);
			}
			j++;
		}
		else if (S(start, start) % S(start, j) == 0) {
			S.col(start).tail(M - start).swap(S.col(j).tail(M - start));
			auto thequotient = S(start, j) / S(start, start);
			S.col(j).tail(M - start) += -thequotient * S.col(start).tail(M - start);
			if (wantQ) {
				Q.col(start).swap(Q.col(j));
				Q.col(j) += -static_cast<Scalar_t<C_t>>(thequotient)* Q.col(start);
				Qi.row(start).swap(Qi.row(j));
				Qi.row(start) += static_cast<Scalar_t<R_t>>(thequotient)* Qi.row(j);
			}
			j++;
			i = start + 1;
		}
		else {
			auto thequotient = floor_division(S(start, j), S(start, start));
			S.col(j).tail(M - start) += -thequotient * S.col(start).tail(M - start);
			S.col(start).swap(S.col(j));
			if (wantQ) {
				Q.col(j) += -static_cast<Scalar_t<C_t>>(thequotient)* Q.col(start);
				Q.col(start).swap(Q.col(j));
				Qi.row(start) += static_cast<Scalar_t<R_t>>(thequotient)* Qi.row(j);
				Qi.row(start).swap(Qi.row(j));
			}
			i = start + 1;
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithFP<S_t, R_t, C_t>::workCol(int start, int& i, int& j) {
		if (S(i, start) % S(start, start) == 0) {
			auto thequotient = S(i, start) / S(start, start);
			S.row(i).tail(N - start) += -thequotient * S.row(start).tail(N - start);
			if (wantP) {
				P.row(i) += -static_cast<Scalar_t<R_t>>(thequotient)* P.row(start);
				Pi.col(start) += static_cast<Scalar_t<C_t>>(thequotient)* Pi.col(i);
			}
			i++;
		}
		else if (S(start, start) % S(i, start) == 0) {
			S.row(start).tail(N - start).swap(S.row(i).tail(N - start));
			auto thequotient = S(i, start) / S(start, start);
			S.row(i).tail(N - start) += -thequotient * S.row(start).tail(N - start);
			if (wantP) {
				P.row(start).swap(P.row(i));
				P.row(i) += -static_cast<Scalar_t<R_t>>(thequotient)* P.row(start);
				Pi.col(start).swap(Pi.col(i));
				Pi.col(start) += static_cast<Scalar_t<C_t>>(thequotient)* Pi.col(i);
			}
			i++;
			j = start + 1;
		}
		else {
			auto thequotient = floor_division(S(i, start), S(start, start));
			S.row(i).tail(N - start) += -thequotient * S.row(start).tail(N - start);
			S.row(start).tail(N - start).swap(S.row(i).tail(N - start));
			if (wantP) {
				P.row(i) += -static_cast<Scalar_t<R_t>>(thequotient)* P.row(start);
				P.row(start).swap(P.row(i));
				Pi.col(start) += static_cast<Scalar_t<C_t>>(thequotient)* Pi.col(i);
				Pi.col(start).swap(Pi.col(i));
			}
			j = start + 1;
		}
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithFP<S_t, R_t, C_t>::searchRow(int start, int& i, int& j)
	{
		if (S(start, start) == 0) {
			bool found = 0;
			for (int k = start + 1; k < N; k++) {
				if (S(start, k) != 0) {
					S.col(start).tail(M - start).swap(S.col(k).tail(M - start));
					if (wantQ) {
						Q.col(start).swap(Q.col(k));
						Qi.row(start).swap(Qi.row(k));
					}
					j = k + 1;
					i = start + 1;
					found = 1;
					break;
				}
			}
			if (!found)
				j = N;
		}
		if (j < N && S(start, j) == 0){
			for (int k = j + 1; k < N; k++){
				if (S(start, k) != 0){
					j = k;
					return;
				}
			}
			j = N; 	/////If you made it here then row is zero outside of start
		}
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithFP<S_t, R_t, C_t>::searchCol(int start, int& i, int& j)
	{
		if (S(start, start) == 0) {
			bool found = 0;
			for (int k = start + 1; k < M; k++) {
				if (S(k, start) != 0) {
					S.row(start).tail(N - start).swap(S.row(k).tail(N - start));
					if (wantP) {
						P.row(start).swap(P.row(k));
						Pi.col(start).swap(Pi.col(k));
					}
					i = k + 1;
					j = start + 1;
					found = 1;
					break;
				}
			}
			if (!found)
				i = M;
		}
		if (i < M && S(i, start) == 0) {
			for (int k = i + 1; k < M; k++) {
				if (S(k, start) != 0) {
					i = k;
					return;
				}
			}
			i = M; /////////If you made it here then column is zero outside of start
		}
	}
}
