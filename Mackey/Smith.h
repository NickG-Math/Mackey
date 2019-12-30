#pragma once
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <type_traits>

///@file
///@brief Contains the Smith normal form algorithms. 
///
///There's two of them, one fast but can overflow (good for finite coefficients) and one slow but reliable (good for infinite coefficients)
///We use SFINAE to detect coefficient type

namespace {

	template<typename Scalar>
	inline Scalar customfloor(Scalar x) {
		if constexpr (std::is_floating_point_v<Scalar>)
			return floor(x);
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
}

namespace Mackey {


	/// The common variables and methods for the two Smith Normal Form implementations
	template <typename S_t, typename R_t, typename C_t>
	class CommonSmith {
		S_t S;		///< The Smith Normal Form of a matrix A.
		R_t P;		///< One of the coefficient matrices (S=P*A*Q). Should be row major for best performance
		R_t Qi;	///< One of the coefficient matrices (A=Pi*S*Qi). Should be row major for best performance
		C_t Q;		///< One of the coefficient matrices (S=P*A*Q). Should be column major for best performance
		C_t Pi;	///< One of the coefficient matrices (A=Pi*S*Qi). Should be column major for best performance
		Eigen::Matrix<typename S_t::Scalar, 1, -1> diagonal; ///<The diagonal of the Smith normal form
		const int L;	///< The length of the diagonal of the Smith normal form.


		/// Initialization of the Smith Normal Form. Refer to the derived class Smith for the implementation
		CommonSmith(const S_t&, bool, bool);

		const int M, N;
		const bool wantP, wantQ;

		/// Sorts the Smith normal form so that it has decreasing entries (with the exception of entries 1,-1) 
		void sorter();
		
		///The implementation of the Smith Normal Form
		template <typename sS_t, typename sR_t, typename sC_t, typename FastOrSlow>
		friend class Smith;
	};

	template<typename S_t, typename R_t, typename C_t>
	CommonSmith<S_t, R_t, C_t> ::CommonSmith(const S_t& A, bool wantP, bool wantQ)
		: S(A), L(std::min(S.rows(), S.cols())), M(S.rows()), N(S.cols()), wantP(wantP), wantQ(wantQ) {
		if (wantP) {
			P = R_t::Identity(M, M);
			Pi = C_t::Identity(M, M);
		}
		if (wantQ) {
			Q = C_t::Identity(N, N);
			Qi = R_t::Identity(N, N);
		}
		diagonal.resize(L);
	}

	template<typename S_t, typename R_t, typename C_t>
	void CommonSmith<S_t, R_t, C_t> ::sorter() {
		std::vector<std::pair<typename S_t::Scalar, int>> toBeSorted;
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
}


namespace Mackey{

	///The safe Smith Normal Form implementation for general coefficients
	template <typename S_t, typename R_t, typename C_t, typename impl=void>
	class Smith : public CommonSmith<S_t, R_t, C_t> {
	public:

		///Computes the Smith Normal Form and sorts if desired
		Smith(const S_t& A, bool wantP, bool wantQ, bool sort) : CommonSmith<S_t, R_t, C_t>(A, wantP, wantQ) {
			SmithIt();
			if (sort)
				this->sorter();
		}
		///Same as the other constructor, but don't sort
		Smith(const S_t& A, bool wantP, bool wantQ) : Smith(A, wantP, wantQ, 0) {}

		//so we don't have to write this-> all the time (and to get these public)
		using CommonSmith<S_t, R_t, C_t>::S;
		using CommonSmith<S_t, R_t, C_t>::P;
		using CommonSmith<S_t, R_t, C_t>::Q;
		using CommonSmith<S_t, R_t, C_t>::Pi;
		using CommonSmith<S_t, R_t, C_t>::Qi;
		using CommonSmith<S_t, R_t, C_t>::diagonal;
		using CommonSmith<S_t, R_t, C_t>::L;

	private:
		void SmithIt();
		void workRow(int);
		void workCol(int);
		void pivot(int);
		void find(int, int&, int&);

		using CommonSmith<S_t, R_t, C_t>::M;
		using CommonSmith<S_t, R_t, C_t>::N;
		using CommonSmith<S_t, R_t, C_t>::wantP;
		using CommonSmith<S_t, R_t, C_t>::wantQ;
	};

	template <typename S_t, typename R_t, typename C_t, typename Slow>
	void Smith<S_t, R_t, C_t, Slow>::SmithIt() {
		bool rowZero, columnZero;
		for (int start = 0; start < this->L; start++) {
			rowZero = S.row(start).tail(N - start - 1).isZero();
			do {
				if (!rowZero) {
					pivot(start);
					workRow(start);
					rowZero = S.row(start).tail(N - start - 1).isZero();
				}
				columnZero = S.col(start).tail(M - start - 1).isZero();
				if (!columnZero) {
					pivot(start);
					workCol(start);
					rowZero = S.row(start).tail(N - start - 1).isZero();
					columnZero = S.col(start).tail(M - start - 1).isZero();
				}
			} while (!rowZero || !columnZero);
			this->diagonal[start] = S(start, start);
		}
	}

	template <typename S_t, typename R_t, typename C_t, typename Slow>
	void Smith<S_t, R_t, C_t,Slow>::workRow(int start) {
		for (int j = start + 1; j < N; j++) {
			if (S(start, j) != 0) {
				auto thequotient = customfloor(S(start, j) / S(start, start));
				S.col(j).tail(M - start) -= thequotient * S.col(start).tail(M - start);
				if (wantQ) {
					Q.col(j) -= static_cast<typename C_t::Scalar>(thequotient)* Q.col(start);
					Qi.row(start) += static_cast<typename R_t::Scalar>(thequotient)* Qi.row(j);
				}
			}
		}
	}

	template <typename S_t, typename R_t, typename C_t, typename Slow>
	void Smith<S_t, R_t, C_t,Slow>::workCol(int start) {
		for (int i = start + 1; i < M; i++) {
			if (S(i, start) != 0) {
				auto thequotient = customfloor(S(i, start) / S(start, start));
				S.row(i).tail(N - start) -= thequotient * S.row(start).tail(N - start);
				if (wantP) {
					P.row(i) -= static_cast<typename R_t::Scalar>(thequotient)* P.row(start);
					Pi.col(start) += static_cast<typename C_t::Scalar>(thequotient)* Pi.col(i);
				}
			}
		}
	}

	template <typename S_t, typename R_t, typename C_t, typename Slow>
	void Smith<S_t, R_t, C_t, Slow>::pivot(int start)
	{
		if (abs(S(start, start)) == 1)
			return;
		int s, t;
		s = t = start;
		find(start, s, t);
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
	}

	template <typename S_t, typename R_t, typename C_t, typename Slow>
	void Smith<S_t, R_t, C_t, Slow>::find(int start, int& s, int& t)
	{
		typename S_t::Scalar min = 0;
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



namespace Mackey {

	///The fast Smith Normal Form implementation for finite coefficients that have a member "order"
	template <typename S_t, typename R_t, typename C_t>
	class Smith<S_t, R_t, C_t, decltype(S_t::Scalar::order, void())> : public CommonSmith<S_t, R_t, C_t> {
	public:
		///Computes the Smith Normal Form and sorts if desired
		Smith(const S_t& A, bool wantP, bool wantQ, bool sort) : CommonSmith<S_t, R_t, C_t>(A, wantP, wantQ) {
			SmithIt();
			if (sort)
				this->sorter();
		}
		///Same as the other constructor, but don't sort
		Smith(const S_t& A, bool wantP, bool wantQ) : Smith(A, wantP, wantQ, 0) {}

		//so we don't have to write this-> all the time (and to get these public)
		using CommonSmith<S_t, R_t, C_t>::S;
		using CommonSmith<S_t, R_t, C_t>::P;
		using CommonSmith<S_t, R_t, C_t>::Q;
		using CommonSmith<S_t, R_t, C_t>::Pi;
		using CommonSmith<S_t, R_t, C_t>::Qi;
		using CommonSmith<S_t, R_t, C_t>::diagonal;
		using CommonSmith<S_t, R_t, C_t>::L;

	private:
		void SmithIt();
		void workRow(const int&, int&, int&);
		void searchRow(const int&, int&, int&);
		void workCol(const int&, int&, int&);
		void searchCol(const int&, int&, int&);

		using CommonSmith<S_t, R_t, C_t>::M;
		using CommonSmith<S_t, R_t, C_t>::N;
		using CommonSmith<S_t, R_t, C_t>::wantP;
		using CommonSmith<S_t, R_t, C_t>::wantQ;
	};


	template <typename S_t, typename R_t, typename C_t>
	void Smith<S_t, R_t, C_t, decltype(S_t::Scalar::order, void())>::SmithIt() {
		for (int start = 0; start < this->L; start++) {
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
	void Smith<S_t, R_t, C_t, decltype(S_t::Scalar::order, void())>::workRow(const int& start, int& i, int& j) {
		if (S(start, j) % S(start, start) == 0) {
			auto thequotient = S(start, j) / S(start, start);
			S.col(j).tail(M - start) += -thequotient * S.col(start).tail(M - start);
			if (wantQ) {
				Q.col(j) += -static_cast<typename C_t::Scalar>(thequotient)* Q.col(start);
				Qi.row(start) += static_cast<typename R_t::Scalar>(thequotient)* Qi.row(j);
			}
			j++;
		}
		else if (S(start, start) % S(start, j) == 0) {
			S.col(start).tail(M - start).swap(S.col(j).tail(M - start));
			auto thequotient = S(start, j) / S(start, start);
			S.col(j).tail(M - start) += -thequotient * S.col(start).tail(M - start);
			if (wantQ) {
				Q.col(start).swap(Q.col(j));
				Q.col(j) += -static_cast<typename C_t::Scalar>(thequotient)* Q.col(start);
				Qi.row(start).swap(Qi.row(j));
				Qi.row(start) += static_cast<typename R_t::Scalar>(thequotient)* Qi.row(j);
			}
			j++;
			i = start + 1;
		}
		else {
			auto thequotient = customfloor(S(start, j) / S(start, start));
			S.col(j).tail(M - start) += -thequotient * S.col(start).tail(M - start);
			S.col(start).swap(S.col(j));
			if (wantQ) {
				Q.col(j) += -static_cast<typename C_t::Scalar>(thequotient)* Q.col(start);
				Q.col(start).swap(Q.col(j));
				Qi.row(start) += static_cast<typename R_t::Scalar>(thequotient)* Qi.row(j);
				Qi.row(start).swap(Qi.row(j));
			}
			i = start + 1;
		}
	}

	template <typename S_t, typename R_t, typename C_t>
	void Smith<S_t, R_t, C_t, decltype(S_t::Scalar::order, void())>::workCol(const int& start, int& i, int& j) {
		if (S(i, start) % S(start, start) == 0) {
			auto thequotient = S(i, start) / S(start, start);
			S.row(i).tail(N - start) += -thequotient * S.row(start).tail(N - start);
			if (wantP) {
				P.row(i) += -static_cast<typename R_t::Scalar>(thequotient)* P.row(start);
				Pi.col(start) += static_cast<typename C_t::Scalar>(thequotient)* Pi.col(i);
			}
			i++;
		}
		else if (S(start, start) % S(i, start) == 0) {
			S.row(start).tail(N - start).swap(S.row(i).tail(N - start));
			auto thequotient = S(i, start) / S(start, start);
			S.row(i).tail(N - start) += -thequotient * S.row(start).tail(N - start);
			if (wantP) {
				P.row(start).swap(P.row(i));
				P.row(i) += -static_cast<typename R_t::Scalar>(thequotient)* P.row(start);
				Pi.col(start).swap(Pi.col(i));
				Pi.col(start) += static_cast<typename C_t::Scalar>(thequotient)* Pi.col(i);
			}
			i++;
			j = start + 1;
		}
		else {
			auto thequotient = customfloor(S(i, start) / S(start, start));
			S.row(i).tail(N - start) += -thequotient * S.row(start).tail(N - start);
			S.row(start).tail(N - start).swap(S.row(i).tail(N - start));
			if (wantP) {
				P.row(i) += -static_cast<typename R_t::Scalar>(thequotient)* P.row(start);
				P.row(start).swap(P.row(i));
				Pi.col(start) += static_cast<typename C_t::Scalar>(thequotient)* Pi.col(i);
				Pi.col(start).swap(Pi.col(i));
			}
			j = start + 1;
		}
	}


	template <typename S_t, typename R_t, typename C_t>
	void Smith<S_t, R_t, C_t, decltype(S_t::Scalar::order, void())>::searchRow(const int& start, int& i, int& j)
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
	void Smith<S_t, R_t, C_t, decltype(S_t::Scalar::order, void())>::searchCol(const int& start, int& i, int& j)
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
