#pragma once
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <type_traits>

///@file
///@brief Contains the Smith normal form algorithm

namespace {

	template<typename Scalar>
	inline Scalar mod(Scalar x, Scalar y) {
		if constexpr (std::is_floating_point_v<Scalar>) {
			return std::fmod(x, y);
		}
		else {
			return x % y;
		}
	}

	template<typename Scalar>
	inline Scalar customfloor(Scalar x) {
		if constexpr (std::is_floating_point_v<Scalar>) {
			return floor(x);
		}
		else {
			return x;
		}
	}

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


}

namespace Mackey {
	/// The Smith Normal Form
	template <typename smithS_t, typename smithR_t, typename smithC_t>
	class Smith {
	public:
		smithS_t S;		///< The Smith Normal Form.
		smithR_t P;		///< One of the coefficient matrices (S=P*A*Q). Should be row major for best performance
		smithR_t Qi;	///< One of the coefficient matrices (A=Pi*S*Qi). Should be row major for best performance
		smithC_t Q;		///< One of the coefficient matrices (S=P*A*Q). Should be column major for best performance
		smithC_t Pi;	///< One of the coefficient matrices (A=Pi*S*Qi). Should be column major for best performance
		Eigen::Matrix<typename smithS_t::Scalar, 1, -1> diagonal; ///<The diagonal of the Smith normal form
		const int L;	///< The length of the diagonal of the Smith normal form.

		/////////////////////////////////////////////////
		/// Given a matrix A compute a diagonal matrix S and invertible matrices P,Q,Pi,Qi such that S=P*A*Q and A=Pi*S*Qi

		///The Pi,Qi are the inverses of P,Q. If wantP=1 then P,Pi are computed and if wantQ=1 then Q,Qi are computed.
		///If sort=1 the diagonal matrix S has decreasing entries, otherwise the order is not specified
		/////////////////////////////////////////////////
		Smith(const smithS_t&, bool, bool, bool);

		///Same as the other constructor, but don't sort
		Smith(const smithS_t& A, bool wantP, bool wantQ) : Smith(A, wantP, wantQ, 0) {}

	private:
		const int M, N;
		const bool wantP, wantQ;
		void SmithIt();
		void workRow(int);
		void workCol(int);
		bool pivot(int);
		bool find(int, int&, int&);
		void sorter();

	};



	template<typename smithS_t, typename smithR_t, typename smithC_t>
	Smith<smithS_t, smithR_t, smithC_t> ::Smith(const smithS_t& Original, bool wantP, bool wantQ, bool sort)
		: S(Original), L(std::min(S.rows(), S.cols())), M(S.rows()), N(S.cols()), wantP(wantP), wantQ(wantQ) {
		if (wantP) {
			P = smithR_t::Identity(M, M);
			Pi = smithC_t::Identity(M, M);
		}
		if (wantQ) {
			Q = smithC_t::Identity(N, N);
			Qi = smithR_t::Identity(N, N);
		}
		diagonal.resize(L);
		SmithIt();
		if (sort)
			sorter();
	}



	template<typename smithS_t, typename smithR_t, typename smithC_t>
	void Smith<smithS_t, smithR_t, smithC_t> ::sorter() {
		std::vector<std::pair<typename smithS_t::Scalar, int>> toBeSorted;
		toBeSorted.reserve(diagonal.size());
		for (int i = 0; i < diagonal.size(); i++) {
			if (abs(diagonal[i]) != 1)
				toBeSorted.push_back(std::make_pair(abs(diagonal[i]), i));
		}
		auto transp = transpositions(toBeSorted);
		for (const auto& i : transp) {
			auto first = toBeSorted[i].second;
			auto second = toBeSorted[i + 1].second;

			auto temp = diagonal[first];
			diagonal[first] = diagonal[second];
			diagonal[second] = temp;

			S.row(first).swap(S.row(second));
			S.col(first).swap(S.col(second));
			if (wantP) {
				P.row(first).swap(P.row(second));
				Pi.col(first).swap(Pi.col(second));
			}
			if (wantQ) {
				Q.col(first).swap(Q.col(second));
				Qi.row(first).swap(Qi.row(second));
			}
		}
	}


	template <typename smithS_t, typename smithR_t, typename smithC_t>
	void Smith<smithS_t, smithR_t, smithC_t>::SmithIt() {
		bool rowZero, columnZero;
		for (int start = 0; start < L; start++) {
			rowZero = S.row(start).tail(N - start - 1).isZero();
			do {
				if (!rowZero) {
					if (!pivot(start))
						return;	
					workRow(start);
					rowZero = S.row(start).tail(N - start - 1).isZero();
				}
				columnZero = S.col(start).tail(M - start - 1).isZero();
				if (!columnZero) {
					if (!pivot(start))
						return;	
					workCol(start);
					rowZero = S.row(start).tail(N - start - 1).isZero();
					columnZero = S.col(start).tail(M - start - 1).isZero();
				}
			} while (!rowZero || !columnZero);
			diagonal[start] = S(start, start);
		}
	}

	template <typename smithS_t, typename smithR_t, typename smithC_t>
	void Smith<smithS_t, smithR_t, smithC_t>::workRow(int start) {
		for (int j = start + 1; j < N; j++) {
			if (S(start, j) != 0) {
				auto thequotient = customfloor(S(start, j) / S(start, start));
				S.col(j).tail(M - start) -= thequotient * S.col(start).tail(M - start);
				if (wantQ) {
					Q.col(j) -= static_cast<typename smithC_t::Scalar>(thequotient)* Q.col(start);
					Qi.row(start) += static_cast<typename smithR_t::Scalar>(thequotient)* Qi.row(j);
				}
			}
		}
	}

	template <typename smithS_t, typename smithR_t, typename smithC_t>
	void Smith<smithS_t, smithR_t, smithC_t>::workCol(int start) {
		for (int i = start + 1; i < M; i++) {
			if (S(i, start) != 0) {
				auto thequotient = customfloor(S(i, start) / S(start, start));
				S.row(i).tail(N - start) -= thequotient * S.row(start).tail(N - start);
				if (wantP) {
					P.row(i) -= static_cast<typename smithR_t::Scalar>(thequotient)* P.row(start);
					Pi.col(start) += static_cast<typename smithC_t::Scalar>(thequotient)* Pi.col(i);
				}
			}
		}
	}

	template <typename smithS_t, typename smithR_t, typename smithC_t>
	bool Smith<smithS_t, smithR_t, smithC_t>::pivot(int start)
	{
		if (abs(S(start, start)) == 1)
			return 1;
		int s, t;
		s = t = start;
		if (!find(start, s, t))
			return 0;
		//bring minimum to the start
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
		return 1;
	}

	template <typename smithS_t, typename smithR_t, typename smithC_t>
	bool Smith<smithS_t, smithR_t, smithC_t>::find(int start, int& s, int& t)
	{
		typename smithS_t::Scalar min = 0;
		for (int j = start; j < S.cols();j++) {
			for (int i = start; i < S.rows(); i++) {
				if (S(i, j) != 0 && (min == 0 || abs(S(i, j)) < min)) {
					min = abs(S(i, j));
					s = i;
					t = j;
					if (min == 1)
						return 1;
				}
			}
		}
		if (min == 0)
			return 0;
		return 1;
	}
}
