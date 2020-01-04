#pragma once
#include "Smith_Base.h"
///@file
///@brief Contains the Smith normal form algorithm for dense matrices where pivots are selected by a Markowitz metric 


namespace {
	template<typename T>
	long absum(const T& a) {
		long sum = abs(a[0]);
		for (int i = 1;i < a.size(); i++) {
			sum += (long)abs(a[i]);
		}
		return sum;
	}
}


namespace Mackey {

	///The dense Smith Normal Form selecting the pivot via a Markowitz metric
	template <typename S_t, typename R_t, typename C_t>
	class SmithSP : Smith<S_t, R_t, C_t> {

		template<typename sS_t, typename sR_t, typename sC_t>
		friend Smith<sS_t, sR_t, sC_t> diagonalize(const sS_t&, bool, bool, bool);

		SmithSP(const S_t& A, bool wantP, bool wantQ, bool sort) : Smith<S_t, R_t, C_t>(A, wantP, wantQ, sort) {
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

		std::vector<long> Rnorm, Cnorm;	///< Norms

		void SmithIt();
		long metric(int, int);
		void initialize_norms();
		void workRow(int);
		void workCol(int);
		void pivoting(int);
		void pivotingRow(int);
		void pivotingCol(int);
		void find_pivot(int, int&, int&);
		void find_pivotRow(int, int&);
		void find_pivotCol(int, int&);
	};

	///The Markowitz metric
	template <typename S_t, typename R_t, typename C_t>
	long SmithSP<S_t, R_t, C_t>::metric(int i, int j) {
		return Rnorm[i] + Cnorm[j]; //Can be changed
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSP<S_t, R_t, C_t>::initialize_norms() {
		Rnorm.reserve(M);
		Cnorm.reserve(N);
		for (int i = 0; i < M; i++)
			Rnorm.push_back(absum(S.row(i)));
		for (int j = 0; j < N; j++)
			Cnorm.push_back(absum(S.col(j)));
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithSP<S_t, R_t, C_t>::SmithIt() {
		initialize_norms();
		for (int start = 0; start < std::min(M,N); start++) {
			pivoting(start);
			bool flagR, flagC;
			flagR = flagC = 0;
			while (Rnorm[start] != abs(S(start, start)) || Cnorm[start] != abs(S(start, start))) {
				while (Rnorm[start] != abs(S(start, start))) {
					if (flagR) {
						flagC = 1;
						pivotingRow(start);
					}
					workRow(start);
					flagR = 1;
				}
				while (Cnorm[start] != abs(S(start, start))) {
					if (flagC) {
						flagR = 1;
						pivotingCol(start);
					}
					workCol(start);
					flagC = 1;
				}
			}
			diagonal[start] = S(start, start);
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSP<S_t, R_t, C_t>::workRow(int start) {
		for (int j = start + 1; j < N; j++) {
			if (S(start, j) != 0) {
				auto thequotient = floor_division(S(start, j), S(start, start));
				for (int i = start; i < M; i++)
					Rnorm[i] -= abs(S(i, j));
				S.col(j).tail(M - start) -= thequotient * S.col(start).tail(M - start);
				Cnorm[j] = absum(S.col(j).tail(M - start));
				for (int i = start; i < M; i++)
					Rnorm[i] += abs(S(i, j));
				if (wantQ) {
					Q.col(j) -= static_cast<Scalar_t<C_t>>(thequotient)* Q.col(start);
					Qi.row(start) += static_cast<Scalar_t<R_t>>(thequotient)* Qi.row(j);
				}
			}
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSP<S_t, R_t, C_t>::workCol(int start) {
		for (int i = start + 1; i < M; i++) {
			if (S(i, start) != 0) {
				auto thequotient = floor_division(S(i, start), S(start, start));
				for (int j = start; j < N; j++)
					Cnorm[j] -= abs(S(i, j));
				S.row(i).tail(N - start) -= thequotient * S.row(start).tail(N - start);
				Rnorm[i] = absum(S.row(i).tail(N - start));
				for (int j = start; j < N; j++)
					Cnorm[j] += abs(S(i, j));
				if (wantP) {
					P.row(i) -= static_cast<Scalar_t<R_t>>(thequotient)* P.row(start);
					Pi.col(start) += static_cast<Scalar_t<C_t>>(thequotient)* Pi.col(i);
				}
			}
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSP<S_t, R_t, C_t>::pivoting(int start)
	{
		int s, t;
		s = t = start;
		find_pivot(start, s, t);
		//bring minimum to the start
		if (t != start) {
			S.col(start).tail(M - start).swap(S.col(t).tail(M - start));
			swap(Cnorm[start], Cnorm[t]);
			if (wantQ) {
				Q.col(start).swap(Q.col(t));
				Qi.row(start).swap(Qi.row(t));
			}
		}
		if (s != start) {
			S.row(start).tail(N - start).swap(S.row(s).tail(N - start));
			swap(Rnorm[start], Rnorm[s]);
			if (wantP) {
				P.row(start).swap(P.row(s));
				Pi.col(start).swap(Pi.col(s));
			}
		}
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithSP<S_t, R_t, C_t>::pivotingRow(int start)
	{
		int t = start;
		find_pivotRow(start, t);
		//bring minimum to the start
		if (t != start) {
			S.col(start).tail(M - start).swap(S.col(t).tail(M - start));
			swap(Cnorm[start], Cnorm[t]);
			if (wantQ) {
				Q.col(start).swap(Q.col(t));
				Qi.row(start).swap(Qi.row(t));
			}
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSP<S_t, R_t, C_t>::pivotingCol(int start)
	{
		int s = start;
		find_pivotCol(start, s);
		//bring minimum to the start
		if (s != start) {
			S.row(start).tail(N - start).swap(S.row(s).tail(N - start));
			swap(Rnorm[start], Rnorm[s]);
			if (wantP) {
				P.row(start).swap(P.row(s));
				Pi.col(start).swap(Pi.col(s));
			}
		}
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithSP<S_t, R_t, C_t>::find_pivot(int start, int& s, int& t)
	{
		typename S_t::Scalar min, minNorm;
		min = minNorm = 0;
		for (int j = start; j < N;j++) {
			for (int i = start; i < M; i++) {
				if (S(i, j) != 0 && (min == 0 || min > abs(S(i, j)) || (min == abs(S(i, j)) && minNorm > metric(i, j)))) {
					min = abs(S(i, j));
					minNorm = metric(i, j);
					s = i;
					t = j;
					if (min == 1 && minNorm == 1)
						return;
				}
			}
		}
	}



	template <typename S_t, typename R_t, typename C_t>
	void SmithSP<S_t, R_t, C_t>::find_pivotRow(int start, int& t)
	{
		typename S_t::Scalar min, minNorm;
		min = minNorm = 0;
		for (int j = start; j < N;j++) {
			if (S(start, j) != 0 && (min == 0 || min > abs(S(start, j)) || (min == abs(S(start, j)) && minNorm > metric(start, j)))) {
				min = abs(S(start, j));
				minNorm = metric(start, j);
				t = j;
				if (min == 1 && minNorm == 1)
					return;
			}
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSP<S_t, R_t, C_t>::find_pivotCol(int start, int& s)
	{
		typename S_t::Scalar min, minNorm;
		min = minNorm = 0;
		for (int i = start; i < M; i++) {
			if (S(i, start) != 0 && (min == 0 || min > abs(S(i, start)) || (min == abs(S(i, start)) && minNorm > metric(i,start)))) {
				min = abs(S(i, start));
				minNorm = metric(i, start);
				s = i;
				if (min == 1 && minNorm == 1)
					return;
			}
		}
	}
}
