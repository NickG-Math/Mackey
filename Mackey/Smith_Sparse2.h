#pragma once
#include "Smith_Base.h"

///@file
///@brief Contains the Smith normal form algorithm for sparse matrices 

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

		///same as S_t but row major
		typedef spm_t_r<S_t> S_t_r;
		Scalar_t<S_t> pivot;						///< The pivot used each turn
		std::vector<int> Rnorm, Cnorm;
		int current;
		S_t_r S_row;

		std::vector<std::pair<Scalar_t<S_t>, std::array<ind, 2>>> Pops, Qops;

		bool isover;
		bool pivotcolunit(ind);
		void do_both(ind);
		int metric(ind, ind);
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
	int SmithSparse<S_t, R_t, C_t>::metric(ind i, ind j) {
		return Rnorm[i] * Cnorm[j]; //Can be changed
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::initialize_norms() {
		Rnorm.resize(M);
		Cnorm.resize(N);
		for (int j = 0; j < S.outerSize(); j++) {
			Cnorm[j] = S.innerVector(j).nonZeros();
			for (typename S_t::InnerIterator it(S, j); it; ++it)
				Rnorm[it.row()]++;
		}
	}


	template <typename S_t, typename R_t, typename C_t>
	bool SmithSparse<S_t, R_t, C_t>::pivotcolunit(ind start)
	{
		int min = 0;
		ind t = start;
		for (ind k = start; k < S.outerSize(); k++) {
			for (typename S_t::InnerIterator it(S, k); it; ++it) {
				if (it.row() == start && abs(it.value()) == 1 && (min == 0 || min > Cnorm[it.col()])) {
					pivot = it.value();
					min = Cnorm[it.col()];
					t = it.col();
				}
			}
		}
		if (min == 0)
			return 0;
		if (start != t) {
			swapCol(S, start, t);
			swap(Cnorm[start], Cnorm[t]);
			if (wantQ)
				Qops.push_back(std::make_pair<Scalar_t<S_t>, std::array<ind, 2>>(0, { start, t }));
		}
		return 1;
	}







	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::SmithIt() {
		initialize_norms();
		if (wantP)
			Pops.reserve(2 * M);
		if (wantQ)
			Qops.reserve(2 * N);
		current = -1;
		for (int start = 0; start < std::min(M, N); start++) {
			pivot = S.coeff(start, start);
			if (Rnorm[start] > 1) {
				auto u = pivotcolunit(start);
				if (u) {	//this clears the entire row and the pivot is pm1 so we can proceed to next iteration. But remember to update Rnorm!
					workRow(start);
					for (typename S_t::InnerIterator it2(S, start); it2; ++it2)
						if (it2.row() > start)
							Rnorm[it2.row()] -= 1;
				}
				else
					do_both(start);
			}
			if ((pivot == 0 && Cnorm[start] > 0) || (Cnorm[start] > 1 && abs(pivot) != 1))
				do_both(start);
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
		std::vector<std::pair<int, typename S_t::Scalar>> nonZeros;
		nonZeros.reserve(N - start);
		for (int k = start + 1; k < S.outerSize(); k++) {
			for (typename S_t::InnerIterator it(S, k); it; ++it) {
				if (it.row() == start)
					nonZeros.push_back(std::make_pair(it.col(), it.value()));
			}
		}
		for (const auto& k : nonZeros) {
			auto j = k.first;
			Scalar_t<S_t> thequotient = floor_division(k.second, pivot);
			for (typename S_t::InnerIterator it2(S, j); it2; ++it2)
				Rnorm[it2.row()] -= 1;
			S.col(j) = (S.col(j) - thequotient * S.col(start)).pruned();
			Cnorm[j] = S.innerVector(j).nonZeros();
			for (typename S_t::InnerIterator it2(S, j); it2; ++it2)
				Rnorm[it2.row()] += 1;
			if (wantQ) {
				std::array<ind, 2> u = { start,j };
				Qops.push_back(std::make_pair(thequotient, u));
			}
		}
		current = -1;
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::do_both(ind start) {
		update(1);
		for (ind i = 0; i < start; i++) {
			if (Cnorm[i] > 1) {
				pivot = this->diagonal[i];
				workCol(i);
			}
		}
		update();
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
	}



	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::workCol(ind start) {
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
					std::array<ind, 2> u = { start,i };
					Pops.push_back(std::make_pair(thequotient, u));
				}
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
				Qops.push_back(std::make_pair<Scalar_t<S_t>, std::array<ind, 2>>(0, { start, t }));
			current = -1;
			swap(Cnorm[start], Cnorm[t]);
		}
		if (s != start) {
			update(1);
			swapRow(S_row, start, s);
			current = 1;
			swap(Rnorm[start], Rnorm[s]);
			if (wantP)
				Pops.push_back(std::make_pair<Scalar_t<S_t>, std::array<ind, 2>>(0, { start, s }));
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
			swap(Cnorm[start], Cnorm[t]);
			if (wantQ)
				Qops.push_back(std::make_pair<Scalar_t<S_t>, std::array<ind, 2>>(0, { start, t }));
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
			swap(Rnorm[start], Rnorm[s]);
			if (wantP)
				Pops.push_back(std::make_pair<Scalar_t<S_t>, std::array<ind, 2>>(0, { start, s }));
		}
	}


	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::find_pivot(ind start, ind& s, ind& t)
	{
		update(-1);
		Scalar_t<S_t> min;
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
	void SmithSparse<S_t, R_t, C_t>::find_pivotCol(ind start, ind& s)
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


	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::renderQ() {

		std::vector<ind> counterQ(N, 1);
		for (const auto& op : Qops) {
			if (op.first == 0)
				swap(counterQ[op.second[0]], counterQ[op.second[1]]);
			else
				counterQ[op.second[1]] = std::min(N, counterQ[op.second[1]] + counterQ[op.second[0]]);
		}
		Q.resize(N, N);
		Qi.resize(N, N);
		Q.reserve(counterQ);
		Qi.reserve(counterQ);
		Q.setIdentity();
		Qi.setIdentity();
		for (const auto& op : Qops) {
			if (op.first == 0) {
				swapCol(Q, op.second[0], op.second[1]);
				swapRow(Qi, op.second[0], op.second[1]);
			}
			else {
				Q.col(op.second[1]) = (Q.col(op.second[1]) - op.first * Q.col(op.second[0])).pruned();
				Qi.row(op.second[0]) = (Qi.row(op.second[0]) + op.first * Qi.row(op.second[1])).pruned();
			}
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void SmithSparse<S_t, R_t, C_t>::renderP() {
		for (ind i = 0; i < std::min(M, N); i++) {
			if (Cnorm[i] > 1) {
				pivot = this->diagonal[i];
				workCol(i);
			}
		}
		std::vector<ind> counterP(M, 1);
		for (const auto& op : Pops) {
			if (op.first == 0)
				swap(counterP[op.second[0]], counterP[op.second[1]]);
			else
				counterP[op.second[1]] = std::min(M, counterP[op.second[1]] + counterP[op.second[0]]);
		}
		P.resize(M, M);
		Pi.resize(M, M);
		P.reserve(counterP);
		Pi.reserve(counterP);
		P.setIdentity();
		Pi.setIdentity();
		for (const auto& op : Pops) {
			if (op.first == 0) {
				swapRow(P, op.second[0], op.second[1]);
				swapCol(Pi, op.second[0], op.second[1]);
			}
			else {
				P.row(op.second[1]) = (P.row(op.second[1]) - op.first * P.row(op.second[0])).pruned();
				Pi.col(op.second[0]) = (Pi.col(op.second[0]) + op.first * Pi.col(op.second[1])).pruned();
			}
		}
	}

}