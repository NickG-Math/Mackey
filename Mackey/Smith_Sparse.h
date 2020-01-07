#pragma once
#include "Smith_Base.h"

///@file
///@brief Contains the Smith normal form algorithm for sparse matrices 

namespace Mackey {
	///The sparse Smith Normal Form, choosing the pivot via a Markowitz metric 
	template <typename S_t, typename R_t, typename C_t>
	class SmithSparse : Smith<S_t,R_t,C_t> {

		template<typename sS_t, typename sR_t, typename sC_t>
		friend Smith<sS_t, sR_t, sC_t> diagonalize(const sS_t&, bool, bool, bool);

		SmithSparse(const S_t& A, bool wantP, bool wantQ, bool sort) : Smith<S_t, R_t, C_t>(A, wantP, wantQ, sort) {
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

		///same as S_t but row major
		typedef Eigen::SparseMatrix<Scalar_t<S_t>, 1> S_t_r;

		Scalar_t<S_t> pivot;						///< The pivot used each turn
		std::vector<int> Rnorm, Cnorm;				
		int current;								 
		Eigen::SparseMatrix<Scalar_t<S_t>, 1> S_row; 

		int metric(int, int);
		void SmithIt();
		void workRow(int);
		void workCol(int);
		void pivoting(int);
		void pivotingRow(int);
		void pivotingCol(int);
		void initialize_norms();
		void find_pivot(int, int&, int&);
		void find_pivotRow(int, int&);
		void find_pivotCol(int, int&);
		void update(int a);
		void update();
	};

	///The Markowitz metric
	template <typename S_t, typename R_t, typename C_t>
	int SmithSparse<S_t, R_t, C_t>::metric(int i, int j) {
		return Rnorm[i]*Cnorm[j]; //Can be changed
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::initialize_norms() {
		Rnorm.resize(M);
		Cnorm.resize(N);
		for (int i = 0; i < S_row.outerSize(); i++)
			Rnorm[i] = S_row.innerVector(i).nonZeros();
		for (int j = 0; j < S.outerSize(); j++)
			Cnorm[j] = S.innerVector(j).nonZeros();
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::SmithIt() {
		S_row = S;
		initialize_norms();
		current = 0;
		for (int start = 0; start < std::min(M, N); start++) {
			if (S.nonZeros() == start) {
				for (int i = start; i < std::min(M, N); i++)
					diagonal[i] = 0;
				return;
			}
			pivoting(start);
			bool flagR, flagC;
			flagR = flagC = 0;
			while (Rnorm[start] != 1 || Cnorm[start] != 1) {
				while (Rnorm[start] != 1) {
					if (flagR) {
						flagC = 1;
						pivotingRow(start);
					}
					workRow(start);
					flagR = 1;
				}
				while (Cnorm[start] != 1) {
					if (flagC) {
						flagR = 1;
						pivotingCol(start);
					}
					workCol(start);
					flagC = 1;
				}
			}
			update();
			this->diagonal[start] = pivot;
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::workRow(int start) {
		update();
		for (typename S_t_r::InnerIterator it(S_row, start); it; ++it) {
			if (it.col() != start) {
				auto j = it.col();
				Scalar_t<S_t> thequotient = floor_division(it.value(), pivot);

				for (typename S_t::InnerIterator it2(S, j); it2; ++it2)
					Rnorm[it2.row()] -= 1;

				S.col(j) = (S.col(j) - thequotient * S.col(start)).pruned();
				Cnorm[j] = S.innerVector(j).nonZeros();
				for (typename S_t::InnerIterator it2(S, j); it2; ++it2) 
					Rnorm[it2.row()] += 1;
				if (wantQ) {
					Q.col(j) = (Q.col(j) - static_cast<Scalar_t<C_t>>(thequotient)* Q.col(start)).pruned();
					Qi.row(start) = (Qi.row(start) + static_cast<Scalar_t<R_t>>(thequotient)* Qi.row(j)).pruned();
				}
			}
		}
		current = -1;
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::workCol(int start) {
		update();
		for (typename S_t::InnerIterator it(S, start); it; ++it) {
			if (it.row() != start) {
				auto i = it.row();
				Scalar_t<S_t> thequotient = floor_division(it.value(), pivot);

				for (typename S_t_r::InnerIterator it2(S_row, i); it2; ++it2)
					Cnorm[it2.col()] -= 1;

				S_row.row(i) = (S_row.row(i) - thequotient * S_row.row(start)).pruned();
				Rnorm[i] = S_row.innerVector(i).nonZeros();				
				for (typename S_t_r::InnerIterator it2(S_row, i); it2; ++it2)
					Cnorm[it2.col()] += 1;
				if (wantP) {
					P.row(i) = (P.row(i) - static_cast<Scalar_t<R_t>>(thequotient) * P.row(start)).pruned();
					Pi.col(start) = (Pi.col(start) + static_cast<Scalar_t<C_t>>(thequotient)* Pi.col(i)).pruned();
				}
			}
		}
		current = 1;
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::pivoting(int start)
	{
		int s, t;
		s = t = start;
		find_pivot(start, s, t);
		if (t != start) {
			update(-1);
			swapCol(S, start, t);
			if (wantQ) {
				swapCol(Q, start, t);
				swapRow(Qi, start, t);
			}
			current = -1;
			swap(Cnorm[start], Cnorm[t]);
		}
		if (s != start) {
			update(1);
			swapRow(S_row, start, s);
			current = 1;
			swap(Rnorm[start], Rnorm[s]);
			if (wantP) {
				swapRow(P, start, s);
				swapCol(Pi, start, s);
			}
		}
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::pivotingRow(int start)
	{
		int t = start;
		find_pivotRow(start, t);
		if (t != start) {
			update(-1);
			swapCol(S, start, t);
			current = -1;
			swap(Cnorm[start], Cnorm[t]);
			if (wantQ) {
				swapCol(Q, start, t);
				swapRow(Qi, start, t);
			}
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::pivotingCol(int start)
	{
		int s = start;
		find_pivotCol(start, s);
		if (s != start) {
			update(1);
			swapRow(S_row, start, s);
			current = 1;
			swap(Rnorm[start], Rnorm[s]);
			if (wantP) {
				swapRow(P, start, s);
				swapCol(Pi, start, s);
			}
		}
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::find_pivot(int start, int& s, int& t)
	{
		update(-1);
		Scalar_t<S_t> min;
		long minNorm;
		min = minNorm = 0;
		for (int k = start; k < S.outerSize(); k++) {
			for (typename S_t::InnerIterator it(S, k); it; ++it) {
				if (it.row() >= start && (((min == 0 || min > abs(it.value())) || (min == abs(it.value()) && minNorm > metric(it.row(),it.col()) )))) {
					min = abs(it.value());
					pivot = it.value();
					s = it.row();
					t = it.col();
					minNorm = metric(s,t);
					if (min == 1 && minNorm == 1)
						return;
				}
			}
		}
	}



	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::find_pivotRow(int start, int& t)
	{
		update(1);
		Scalar_t<S_t> min;
		long minNorm;
		min = minNorm = 0;
		for (typename S_t_r::InnerIterator it(S_row, start); it; ++it) {
			if (it.col() >= start && (((min == 0 || min > abs(it.value())) || (min == abs(it.value()) && minNorm > metric(start, it.col()))))) {
				min = abs(it.value());
				pivot = it.value();
				t = it.col();
				minNorm = metric(start, t);
				if (min == 1 && minNorm == 1)
					return;
			}
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::find_pivotCol(int start, int& s)
	{
		update(-1);
		Scalar_t<S_t> min;
		long minNorm;
		min = minNorm = 0;
		for (typename S_t::InnerIterator it(S, start); it; ++it) {
			if (it.row() >= start && (((min == 0 || min > abs(it.value())) || (min == abs(it.value()) && minNorm > metric(it.row(), start))))) {
				min = abs(it.value());
				pivot = it.value();
				s = it.row();
				minNorm = metric(s, start);
				if (min == 1 && minNorm == 1)
					return;
			}
		}
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::update(int a) {
		if (a == 1 && current == -1) {
			S_row = S;
			current = 0;
		}
		else if (a == -1 && current == 1) {
			S = S_row;
			current = 0;
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::update() {
		update(1);
		update(-1);
	}
}
