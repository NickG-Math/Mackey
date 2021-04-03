#pragma once
#include "../Smith.hpp"
///@file
///@brief Contains the base class for all Smith implementations. 

namespace mackey {


	namespace implementation_details { //lots of helper functions

		//Because C++ doesn't floor the division of integers...
		template<typename Scalar>
		inline Scalar floor_division(Scalar x, Scalar y) {
			if constexpr (SFINAE::is_finite_cyclic<Scalar>::value)
				return (x / y);
			else if constexpr (std::is_integral_v<Scalar>) { //C++ rounds integer division but we need floor!
				auto q = x / y;
				auto r = x % y;
				return (r == 0) ? q : q - ((x < 0) && !(y < 0));
			}
		}

		///Produces the transpositions that sort a given vector of pairs by the first index only
		template<typename T, typename ind>
		std::vector<ind> transpositions(std::vector<std::pair<T, ind>> a) {
			std::vector<ind> transpositions;
			if (a.size() <= 1)
				return transpositions;
			transpositions.reserve((a.size() - 1) * (a.size() - 1));
			for (size_t j = 0; j < a.size() - 1; j++) {
				for (size_t i = 0; i < a.size() - j - 1; i++) {
					if (a[i].first < a[i + 1].first) {
						transpositions.push_back(i);
						std::swap(a[i], a[i + 1]);
					}
				}
			}
			return transpositions;
		}

		template<typename T, typename storage>
		void swapCol(Eigen::SparseMatrix<T, 0, storage>& A, storage i, storage j) {
			Eigen::SparseMatrix<T, 0, storage> Ai = A.col(i);
			Eigen::SparseMatrix<T, 0, storage> Aj = A.col(j);
			A.col(i) = Aj;
			A.col(j) = Ai;
		}

		template<typename T, typename storage>
		void swapRow(Eigen::SparseMatrix<T, 1, storage>& A, storage i, storage j) {
			Eigen::SparseMatrix<T, 1, storage> Ai = A.row(i);
			Eigen::SparseMatrix<T, 1, storage> Aj = A.row(j);
			A.row(i) = Aj;
			A.row(j) = Ai;
		}

		template<typename T, typename ind>
		void swapCol(Eigen::MatrixBase<T>& A, ind i, ind j, bool tail = 0) {
			if (tail)
				A.col(i).tail(A.rows() - i).swap(A.col(j).tail(A.rows() - i));
			else
				A.col(i).swap(A.col(j));
		}

		template<typename T, typename ind>
		void swapRow(Eigen::MatrixBase<T>& A, ind i, ind j, bool tail = 0) {
			if (tail)
				A.row(i).tail(A.cols() - i).swap(A.row(j).tail(A.cols() - i));
			else
				A.row(i).swap(A.row(j));
		}

		template<typename RCO, typename ind>
		std::vector<ind> expected_size(const RCO& ops, ind size) {
			std::vector<ind> counter(size, 1);
			for (const auto& op : ops) {
				if (op.scalar == 0)
					std::swap(counter[op.start], counter[op.end]);
				else
					counter[op.end] = std::min(size, counter[op.end] + counter[op.start]);
			}
			return counter;
		}

		template<bool row, bool inverse_operation, typename T, typename RCO, typename S>
		void render(T& mat, const RCO& ops, const S& counter) {
			mat.reserve(counter);
			mat.setIdentity();
			for (const auto& op : ops) {
				auto start = op.start;
				auto end = op.end;
				auto scalar = op.scalar;
				if constexpr (inverse_operation) {
					std::swap(start, end);
					scalar = -scalar;
				}
				if (op.scalar == 0) {
					if constexpr (row)
						swapRow(mat, start, end);
					else
						swapCol(mat, start, end);
				}
				else {
					if constexpr (row)
						mat.row(start) = (mat.row(start) + scalar * mat.row(end));
					else
						mat.col(start) = (mat.col(start) + scalar * mat.col(end));
				}
			}
		}

	}


	template<typename _S, typename _R, typename _C>
	SmithNormalForm<_S, _R, _C>::SmithNormalForm(const _S& A, bool do_P, bool do_Q, bool do_sort, bool do_verify)
		: S(A), M(A.rows()), N(A.cols()), doP(do_verify ? 1 : do_P), doQ(do_verify ? 1 : do_Q) {
		diagonal.resize(std::min(M, N));
		if (doP) {
			P.resize(M, M);
			Pi.resize(M, M);
		}
		if (doQ) {
			Q.resize(N, N);
			Qi.resize(N, N);
		}
		compute();
		if constexpr (!partial_pivoting) {//sort for non fields: for a field, all values are invertible or 0
			if (do_sort)
				sorter();
		}
		if (do_verify)
			verify(A);
	}

	template<typename _S, typename _R, typename _C>
	void SmithNormalForm<_S, _R, _C>::compute() {
		initialize();
		for (ind start = 0; start < std::min(M, N); start++) {
			bool pivot_found = initial_pivoting_per_loop(start);
			if (!pivot_found) {
				for (auto i = start; i < std::min(M, N); i++)
					diagonal[i] = 0;
				break;
			}
			if constexpr (partial_pivoting) {
				if (!doneRow(start))
					eliminateRow(start);
			}
			else {
				bool bestpivot = 1;
				while (!doneRow(start) || !doneCol(start)) {
					while (!doneRow(start)) {
						if (!bestpivot) {
							find_and_set_pivot<1, 0>(start); //find pivot in row start
							col_pivot(start); //move pivot to start,start via column operation
						}
						eliminateRow(start);
						bestpivot = 0;
					}
					while (!doneCol(start)) {
						if (!bestpivot) {
							find_and_set_pivot<0, 1>(start); //find pivot in column start
							row_pivot(start); //move pivot to start,start via row operation
						}
						eliminateCol(start);
						bestpivot = 0;
					}
				}
			}
			this->diagonal[start] = piv.value;
		}
		renderPQ();
	}

	template<typename _S, typename _R, typename _C>
	void SmithNormalForm<_S, _R, _C>::sorter() {
		std::vector<std::pair<scalar_t<S_t>, ind>> toBeSorted;
		toBeSorted.reserve(diagonal.size());
		for (ind i = 0; i < diagonal.size(); i++) {
			if (abs(diagonal[i]) != 1)
				toBeSorted.push_back(std::make_pair(abs(diagonal[i]), i));
		}
		auto transp = implementation_details::transpositions(toBeSorted);
		for (const auto& i : transp) {
			auto first = toBeSorted[i].second;
			auto second = toBeSorted[i + 1].second;
			std::swap(diagonal[first], diagonal[second]);
			if (doP) {
				implementation_details::swapRow(P, first, second);
				implementation_details::swapCol(Pi, first, second);
			}
			if (doQ) {
				implementation_details::swapCol(Q, first, second);
				implementation_details::swapRow(Qi, first, second);
			}
		}
	}

	template<typename _S, typename _R, typename _C>
	void SmithNormalForm<_S, _R, _C>::verify(const S_t& A) const {
		S_t R = P * A * Q;
		if constexpr (SFINAE::is_Sparse<S_t>::value)
			R.prune(0, 0);
		bool err = 0;
		for (IteratorNNZ<S_t, 0, 0> it(R, 0); it; ++it) {
			if (it.row() != it.col()) {
				err = 1;
				std::cerr << "Verification Failed: Not Diagonal:\n" << R << "\n\n";
			}
			if (it.value() != diagonal[it.col()]) {
				err = 1;
				std::cerr << "Verification Failed: Wrong pivot at: " << it.col() << " Compare\n" << R << "\n\n" << "and: \n" << diagonal << "\n\n";
			}
			if (err) {
				std::cerr << "Original matrix:\n" << A << "\n\n";
				std::cerr << "P matrix:\n" << P << "\n\n";
				std::cerr << "Q matrix:\n" << Q << "\n\n";
				abort();
			}
		}
	}

	template <typename _S, typename _R, typename _C>
	int64_t SmithNormalForm<_S, _R, _C>::metric(ind i, ind j) const {
		if constexpr (R_exists && C_exists)
			return this->Rnorm[i] * this->Cnorm[j];
		else if constexpr (C_exists)
			return this->Cnorm[j];
		else
			return S.col(j).nonZeros();
	}

	template <typename _S, typename _R, typename _C>
	void SmithNormalForm<_S, _R, _C>::initialize() {
		if constexpr (!sparse) {
			if (doP) {
				P.setIdentity();
				Pi.setIdentity();
			}
			if (doQ) {
				Q.setIdentity();
				Qi.setIdentity();
			}
		}
		else {
			if (doP)
				this->Pops.reserve(2 * M);
			if (doQ)
				this->Qops.reserve(2 * N);
			if constexpr (!partial_pivoting) {
				this->S_row = S;
				this->current = 0;
			}
		}
		initialize_norms();
	}


	template <typename _S, typename _R, typename _C>
	void SmithNormalForm<_S, _R, _C>::row_pivot(ind start, bool tail) {
		if (piv.row == start)
			return;
		if constexpr (R_exists)
			std::swap(this->Rnorm[start], this->Rnorm[piv.row]);
		if constexpr (!sparse) {
			implementation_details::swapRow(S, start, piv.row, tail);
			if (doP) {
				implementation_details::swapRow(P, start, piv.row);
				implementation_details::swapCol(Pi, start, piv.row);
			}
		}
		else {
			if constexpr (!partial_pivoting) {
				update(1);
				implementation_details::swapRow(this->S_row, start, piv.row);
				this->current = 1;
			}
			else {
				Eigen::SparseMatrix<scalar_t<_S>, 1, ind> S_row = S;
				implementation_details::swapRow(S_row, start, piv.row);
				S = S_row;
			}
			if (doP)
				this->Pops.emplace_back(0, start, piv.row);
		}
	}

	template <typename _S, typename _R, typename _C>
	void SmithNormalForm<_S, _R, _C>::col_pivot(ind start, bool tail) {
		if (piv.col == start)
			return;
		if constexpr (C_exists)
			std::swap(this->Cnorm[start], this->Cnorm[piv.col]);
		if constexpr (!sparse) {
			implementation_details::swapCol(S, start, piv.col, tail);
			if (doQ) {
				implementation_details::swapCol(Q, start, piv.col);
				implementation_details::swapRow(Qi, start, piv.col);
			}
		}
		else {
			if constexpr (!partial_pivoting) {
				update(-1);
				implementation_details::swapCol(S, start, piv.col);
				this->current = -1;
			}
			else
				implementation_details::swapCol(S, start, piv.col);
			if (doQ)
				this->Qops.emplace_back(0, start, piv.col);
		}
	}


	template <typename _S, typename _R, typename _C>
	bool SmithNormalForm<_S, _R, _C>::initial_pivoting_per_loop(ind start) {
		if constexpr (!partial_pivoting && sparse) {
			update(-1);
			if (S.nonZeros() == start)
				return 0;
		}
		piv.value = 0;
		piv.row = start;
		piv.col = start;
		find_and_set_pivot<0, 0>(start);
		if (!partial_pivoting && !sparse && piv.value == 0) //when full pivoting, a zero pivot means that the entire matrix is 0.
			return 0;
		col_pivot(start, 1); //do tail optimization (in dense case)
		row_pivot(start, !partial_pivoting); //only do tail optimization (in dense case) when full pivoting
		return 1;
	}

	template <typename _S, typename _R, typename _C>
	bool SmithNormalForm<_S, _R, _C>::doneRow(ind start) {
		if constexpr (!partial_pivoting && R_exists)
			return (this->Rnorm[start] <= 1);
		else if constexpr (!partial_pivoting)
			return (this->S_row.row(start).nonZeros() <= 1);
		else //when partial pivoting use the variable donerow
			return this->donerow;
	}

	template <typename _S, typename _R, typename _C>
	bool SmithNormalForm<_S, _R, _C>::doneCol(ind start) {
		if constexpr (C_exists)
			return (this->Cnorm[start] <= 1);
		else
			return (S.col(start).nonZeros() <= 1);
	}

	template <typename _S, typename _R, typename _C>
	template<bool onlyrow, bool onlycol>
	void SmithNormalForm<_S, _R, _C>::find_and_set_pivot(ind start) {
		int64_t minNorm = 0;
		if constexpr (!partial_pivoting) { //full pivoting. We follow the instructions in the template parameters onlyrow and onlycolumn as to find the pivot either in the whole matrix, or in the row/column start.
			typedef std::conditional_t<!(sparse&& onlyrow), S_t, S_row_t> matrix_t;
			const matrix_t* A; //is S or S_row depending on whether we want to traverse the whole matrix (S), the column start (S) or the row start (S_row)
			if constexpr (!sparse)
				A = &S;
			else {
				if constexpr (onlyrow) {
					update(1);
					A = &this->S_row;
				}
				else {
					update(-1);
					A = &S;
				}
			}
			scalar_t<S_t> min = 0;
			for (IteratorNNZ<matrix_t, onlyrow, onlycol> it(*A, start); it; ++it)
				if (min == 0 || min > abs(it.value()) || (min == abs(it.value()) && minNorm > metric(it.row(), it.col()))) { //best means = nonzero, smallest absolute value, smallest norm
					min = abs(it.value());
					minNorm = metric(it.row(), it.col());
					piv = it;
					if (min == 1 && minNorm == 1)
						break;
				}
		}
		else { //partial pivoting, we ignore the instructions and try to find the pivot in the row if possible
			this->donerow = 0;
			piv.value = 0;
			for (IteratorNNZ<S_t, !sparse, 0> it(S, start); it; ++it)
				if (it.row() == start && (piv.value == 0 || minNorm > metric(it.row(), it.col()))) { //best means = nonzero, smallest absolute value, smallest norm
					minNorm = metric(it.row(), it.col());
					piv = it;
					if (minNorm == 1)
						break;
				}
			if (piv.value == 0) { //pivot not found as the whole row is 0
				if (metric(start, start) != 0) { //the column start isn't so we must search in that column. Since we do not have the Rnorm available, we just find the first nonzero
					IteratorNNZ<S_t, 0, 1> it(S, start);
					piv = it;
				}
				else
					this->donerow = 1; //proceed to next iteration			
			}

		}
	}

	template <typename _S, typename _R, typename _C>
	void SmithNormalForm<_S, _R, _C>::initialize_norms() { //only in full pivoting
		if constexpr (R_exists)
			this->Rnorm.resize(M);
		if constexpr (C_exists)
			this->Cnorm.resize(N);
		if constexpr (!sparse)
			update_norms<1, 1, 1>(IteratorNNZ<S_t, 0, 0>(S, 0, 0));
		else {
			for (ind i = 0; i < M; i++)
				update_norms<1, 0, 0>(i);
			for (ind j = 0; j < N; j++)
				update_norms<0, 1, 0>(j);
		}
	}

	template <typename _S, typename _R, typename _C>
	template<bool row, bool col, char increase_decrease_zero, typename iter>
	void SmithNormalForm<_S, _R, _C>::update_norms(iter it) { //only in full pivoting
		constexpr bool doR = row && R_exists;
		constexpr bool doC = col && C_exists;
		if constexpr (std::is_integral_v<iter>) {
			if constexpr (doR)
				this->Rnorm[it] = this->S_row.row(it).nonZeros();
			if constexpr (doC)
				this->Cnorm[it] = S.col(it).nonZeros();
		}
		else {
			if constexpr (increase_decrease_zero == 0) {
				if constexpr (doR)
					this->Rnorm[it.row()] = 0;
				if constexpr (doC)
					this->Cnorm[it.col()] = 0;
			}
			else if constexpr (doR || doC) {
				for (;it;++it) {
					if constexpr (doR && increase_decrease_zero == 1)
						this->Rnorm[it.row()]++;
					else if constexpr (doR && increase_decrease_zero == -1)
						this->Rnorm[it.row()]--;
					if constexpr (doC && increase_decrease_zero == 1)
						this->Cnorm[it.col()]++;
					else if constexpr (doC && increase_decrease_zero == -1)
						this->Cnorm[it.col()]--;
				}
			}
		}
	}

	template <typename _S, typename _R, typename _C>
	void SmithNormalForm<_S, _R, _C>::eliminateCol(ind start) { //only in full pivoting
		if constexpr (!sparse) {
			for (IteratorNNZ<S_t, 0, 1> it(S, start + 1, start); it; ++it) { //iterate through nonzeros in column=start, beginning with row=start+1
				auto thequotient = implementation_details::floor_division(it.value(), piv.value);
				update_norms<0, 1, -1>(IteratorNNZ<S_t, 1, 0>(S, it.row(), start)); //decreases the relevant Cnorms by 1 in anticipation of what comes next
				S.row(it.row()).tail(N - start) -= thequotient * S.row(start).tail(N - start);
				update_norms<1, 0, 0>(it); //sets the relevant Rnorm to 0 in anticipation of what comes next
				update_norms<1, 1, 1>(IteratorNNZ<S_t, 1, 0>(S, it.row(), start)); //fixes the norms
				if (doP) {
					P.row(it.row()) -= thequotient * P.row(start);
					Pi.col(start) += thequotient * Pi.col(it.row());
				}
			}
		}
		else {
			update(1);
			for (ind i = start + 1; i < this->S_row.outerSize(); i++) {
				typename S_row_t::InnerIterator it(this->S_row, i);
				if (it && it.col() == start) {
					auto thequotient = implementation_details::floor_division(it.value(), piv.value);
					update_norms<0, 1, -1>(it);
					this->S_row.row(i) = (this->S_row.row(i) - thequotient * this->S_row.row(start)).pruned();
					update_norms<1, 0, 0>(i);
					update_norms<0, 1, 1>(typename S_row_t::InnerIterator(this->S_row, i));
					if (doP)
						this->Pops.emplace_back(thequotient, start, i);
				}
			}
			this->current = 1;
		}
	}

	template <typename _S, typename _R, typename _C>
	void SmithNormalForm<_S, _R, _C>::eliminateRow(ind start) {
		if constexpr (!sparse) {
			for (IteratorNNZ<S_t, 1, 0> it(S, start, start + 1); it; ++it) { //iterate through nonzeros in row=start, beginning with col=start+1
				auto thequotient = implementation_details::floor_division(it.value(), piv.value);
				update_norms<1, 0, -1>(IteratorNNZ<S_t, 0, 1>(S, start, it.col()));
				S.col(it.col()).tail(M - start) -= thequotient * S.col(start).tail(M - start);
				update_norms<0, 1, 0>(it);
				update_norms<1, 1, 1>(IteratorNNZ<S_t, 0, 1>(S, start, it.col()));
				if (doQ) {
					Q.col(it.col()) -= thequotient * Q.col(start);
					Qi.row(start) += thequotient * Qi.row(it.col());
				}
			}
		}
		else {
			if constexpr (!partial_pivoting)
				update(-1);
			for (ind j = start + 1; j < S.outerSize(); j++) {
				typename S_t::InnerIterator it(S, j);
				if (it && it.row() == start) {
					scalar_t<S_t> thequotient = implementation_details::floor_division(it.value(), piv.value);
					update_norms<1, 0, -1>(it);
					S.col(j) = (S.col(j) - thequotient * S.col(start)).pruned();
					update_norms<0, 1, 0>(j);
					update_norms<1, 0, 1>(typename S_t::InnerIterator(S, j));
					if (doQ)
						this->Qops.emplace_back(thequotient, start, j);
				}
			}
			if constexpr (!partial_pivoting)
				this->current = -1;
		}
	}



	template <typename _S, typename _R, typename _C>
	void SmithNormalForm<_S, _R, _C>::renderPQ() { //only in partial pivoting + dense or any pivoting +sparse
		if (doQ)
			if constexpr (sparse) {
				auto counter = implementation_details::expected_size(this->Qops, N);
				implementation_details::render<0, 1>(Q, this->Qops, counter);
				implementation_details::render<1, 0>(Qi, this->Qops, counter);
			}
		if (doP) {
			if constexpr (partial_pivoting) { //S is only lower triangular.	We virtually perform row operations on S to get P and Pi
				for (ind start = 0; start < std::min(M, N); start++)
					if (!doneCol(start))
						for (IteratorNNZ<S_t, 0, 1> it(S, start + 1, start); it; ++it) {
							auto thequotient = it.value() / diagonal[start];
							if constexpr (!sparse) {
								P.row(it.row()) -= thequotient * P.row(start);
								Pi.col(start) += thequotient * Pi.col(it.col());
							}
							else
								this->Pops.emplace_back(thequotient, start, it.row());
						}
			}
			if constexpr (sparse) {
				auto counter = implementation_details::expected_size(this->Pops, M);
				implementation_details::render<1, 1>(P, this->Pops, counter);
				implementation_details::render<0, 0>(Pi, this->Pops, counter);
			}
		}
	}

	template <typename _S, typename _R, typename _C>
	void SmithNormalForm<_S, _R, _C>::update(int a) {
		if (a == 0) {
			update(1);
			update(-1);
		}
		else if (a == 1 && this->current == -1) {
			this->S_row = S;
			this->current = 0;
		}
		else if (a == -1 && this->current == 1) {
			S = this->S_row;
			this->current = 0;
		}
	}
}

//alternate eliminateCol for full pivot sparse when both are updated all the time (in total it's slower)
//update();
//typename S_t::InnerIterator it(S, start);
//for (++it;it;++it) {
//	auto thequotient = floor_division(it.value(), piv.value);
//	S_row.row(it.row()) -= thequotient * S_row.row(start);
//	if (doP)
//		Pops.emplace_back(thequotient, start, it.row());
//}
//S_row.prune(0, 0);
//S = S_row;
//current = 0;

//same for eliminateRow

//else if constexpr (!partial_pivoting) {
//update();
//typename S_row_t::InnerIterator it(S_row, start);
//for (++it;it;++it) {
//	auto thequotient = floor_division(it.value(), piv.value);
//	S.col(it.col()) -= thequotient * S.col(start);
//	if (doQ)
//		Qops.emplace_back(thequotient, start, it.col());
//}
//S.prune(0, 0);
//it = typename S_row_t::InnerIterator(S_row, start);
//for (++it;it;++it)
//	Cnorm[it.col()] = S.col(it.col()).nonZeros();
//S_row = S;
//for (ind i = 0; i < M; i++)
//	Rnorm[i] = S_row.row(i).nonZeros();
//current = 0;
//}