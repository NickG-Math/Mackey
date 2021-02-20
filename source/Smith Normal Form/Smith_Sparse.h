#pragma once
#include "Smith_Base.h"

///@file
///@brief Contains the Smith normal form algorithm for sparse matrices 

namespace {
	template<typename T, typename S>
	struct Row_Column_Operation {
		//if scalar=0 then this indicates a swap of start,end
		//otherwise this indicates an operation end->end-scalar*start
		T scalar;
		S start;
		S end;
		Row_Column_Operation() {}
		Row_Column_Operation(T scalar, S start, S end) :scalar(scalar), start(start), end(end) {}
	};
}

namespace Mackey {
	///The sparse Smith Normal Form, choosing the pivot via a Markowitz metric 
	template <typename S_t, typename R_t, typename C_t>
	class SmithSparse : Smith<S_t, R_t, C_t> {

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

		using typename Smith<S_t, R_t, C_t>::ind;


		///same as S_t but row major
		typedef typename Eigen::SparseMatrix<scalar_t<S_t>, 1, typename S_t::StorageIndex> S_row_t;

		scalar_t<S_t> pivot;						///< The pivot used each turn
		std::vector<int> Rnorm, Cnorm;
		int current;
		S_row_t S_row;

		std::vector<Row_Column_Operation<scalar_t<S_t>, ind>> Pops, Qops;
		long metric(ind, ind) const;
		void SmithIt();
		void workRow(ind);
		void workCol(ind);
		void pivoting(ind);
		void pivotingRow(ind);
		void pivotingCol(ind);
		void initialize_norms();
		void find_pivot(ind, ind&, ind&);
		void find_pivotRow(ind, ind&);
		void find_pivotCol(ind, ind&);
		void update(int a);
		void update();
		void renderP();
		void renderQ();
	};

	///The Markowitz metric
	template <typename S_t, typename R_t, typename C_t>
	long SmithSparse<S_t, R_t, C_t>::metric(ind i, ind j) const {
		return Rnorm[i] * Cnorm[j]; //Can be changed
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::initialize_norms() {
		Rnorm.resize(M);
		Cnorm.resize(N);
		for (ind i = 0; i < S_row.outerSize(); i++)
			Rnorm[i] = S_row.innerVector(i).nonZeros();
		for (ind j = 0; j < S.outerSize(); j++)
			Cnorm[j] = S.innerVector(j).nonZeros();
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::SmithIt() {

		S_row = S;
		initialize_norms();
		if (wantP)
			Pops.reserve(2 * M);
		if (wantQ)
			Qops.reserve(2 * N);

		current = 0;
		for (ind start = 0; start < std::min(M, N); start++) {
			if (S.nonZeros() == start) {
				for (ind i = start; i < std::min(M, N); i++)
					diagonal[i] = 0;
				break;
			}
			pivoting(start);
			bool flagR, flagC;
			flagR = flagC = 0;
			while (Rnorm[start] > 1 || Cnorm[start] > 1) {
				while (Rnorm[start] > 1) {
					if (flagR) {
						flagC = 1;
						pivotingRow(start);
					}
					workRow(start);
					flagR = 1;
				}
				while (Cnorm[start] > 1) {
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
		if (wantP)
			renderP();
		if (wantQ)
			renderQ();
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::workRow(ind start) {
		update(-1);
		for (ind j = start + 1; j < S.outerSize(); j++) {
			typename S_t::InnerIterator it(S, j);
			if (it && it.row() == start) {
				scalar_t<S_t> thequotient = floor_division(it.value(), pivot);
				for (; it; ++it)
					Rnorm[it.row()]--;
				S.col(j) = (S.col(j) - thequotient * S.col(start)).pruned();
				Cnorm[j] = S.col(j).nonZeros();
				for (typename S_t::InnerIterator it2(S, j); it2; ++it2)
					Rnorm[it2.row()]++;
				if (wantQ)
					Qops.emplace_back(thequotient, start, j);
			}
		}
		current = -1;
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::workCol(ind start) {
		update(1);
		for (ind i = start + 1; i < S_row.outerSize(); i++) {
			typename S_row_t::InnerIterator it(S_row, i);
			if (it && it.col() == start) {
				scalar_t<S_t> thequotient = floor_division(it.value(), pivot);
				for (; it; ++it)
					Cnorm[it.col()] --;
				S_row.row(i) = (S_row.row(i) - thequotient * S_row.row(start)).pruned();
				Rnorm[i] = S_row.row(i).nonZeros();
				for (typename S_row_t::InnerIterator it2(S_row, i); it2; ++it2)
					Cnorm[it2.col()] ++;
				if (wantP)
					Pops.emplace_back(thequotient, start, i);
			}
		}
		current = 1;
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::pivoting(ind start)
	{
		ind s, t;
		s = t = start;
		find_pivot(start, s, t);
		if (t != start) {
			update(-1);
			swapCol(S, start, t);
			if (wantQ)
				Qops.emplace_back(0, start, t);
			current = -1;
			std::swap(Cnorm[start], Cnorm[t]);
		}
		if (s != start) {
			update(1);
			swapRow(S_row, start, s);
			current = 1;
			std::swap(Rnorm[start], Rnorm[s]);
			if (wantP)
				Pops.emplace_back(0, start, s);
		}
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::pivotingRow(ind start)
	{
		ind t = start;
		find_pivotRow(start, t);
		if (t != start) {
			update(-1);
			swapCol(S, start, t);
			current = -1;
			std::swap(Cnorm[start], Cnorm[t]);
			if (wantQ)
				Qops.emplace_back(0, start, t);
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::pivotingCol(ind start)
	{
		ind s = start;
		find_pivotCol(start, s);
		if (s != start) {
			update(1);
			swapRow(S_row, start, s);
			current = 1;
			std::swap(Rnorm[start], Rnorm[s]);
			if (wantP)
				Pops.emplace_back(0, start, s);
		}
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::find_pivot(ind start, ind& s, ind& t)
	{
		update(-1);
		scalar_t<S_t> min;
		long minNorm;
		min = minNorm = 0;
		for (ind k = start; k < S.outerSize(); k++) {
			for (typename S_t::InnerIterator it(S, k); it; ++it) {
				if (it.row() >= start && (((min == 0 || min > abs(it.value())) || (min == abs(it.value()) && minNorm > metric(it.row(), it.col()))))) {
					min = abs(it.value());
					pivot = it.value();
					s = it.row();
					t = it.col();
					minNorm = metric(s, t);
					if (min == 1 && minNorm == 1)
						return;
				}
			}
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::find_pivotRow(ind start, ind& t)
	{
		update(1);
		scalar_t<S_t> min;
		long minNorm;
		min = minNorm = 0;
		for (typename S_row_t::InnerIterator it(S_row, start); it; ++it) {
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
	void SmithSparse<S_t, R_t, C_t>::find_pivotCol(ind start, ind& s)
	{
		update(-1);
		scalar_t<S_t> min;
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


	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::renderQ() {

		std::vector<ind> counterQ(N, 1);
		for (const auto& op : Qops) {
			if (op.scalar == 0)
				std::swap(counterQ[op.start], counterQ[op.end]);
			else
				counterQ[op.end] = std::min(N, counterQ[op.end] + counterQ[op.start]);
		}
		Q.resize(N, N);
		Qi.resize(N, N);
		Q.reserve(counterQ);
		Qi.reserve(counterQ);
		Q.setIdentity();
		Qi.setIdentity();
		for (const auto& op : Qops) {
			if (op.scalar == 0) {
				swapCol(Q, op.start, op.end);
				swapRow(Qi, op.start, op.end);
			}
			else {
				Q.col(op.end) = (Q.col(op.end) - op.scalar * Q.col(op.start)).pruned();
				Qi.row(op.start) = (Qi.row(op.start) + op.scalar * Qi.row(op.end)).pruned();
			}
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::renderP() {
		std::vector<ind> counterP(M, 1);
		for (const auto& op : Pops) {
			if (op.scalar == 0)
				std::swap(counterP[op.start], counterP[op.end]);
			else
				counterP[op.end] = std::min(M, counterP[op.end] + counterP[op.start]);
		}
		P.resize(M, M);
		Pi.resize(M, M);
		P.reserve(counterP);
		Pi.reserve(counterP);
		P.setIdentity();
		Pi.setIdentity();
		for (const auto& op : Pops) {
			if (op.scalar == 0) {
				swapRow(P, op.start, op.end);
				swapCol(Pi, op.start, op.end);
			}
			else {
				P.row(op.end) = (P.row(op.end) - op.scalar * P.row(op.start)).pruned();
				Pi.col(op.start) = (Pi.col(op.start) + op.scalar * Pi.col(op.end)).pruned();
			}
		}
	}
}
