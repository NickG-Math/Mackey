#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <type_traits>

///@file
///@brief Contains the Smith normal form algorithm

namespace {

	template<typename Scalar>
	inline Scalar mod(Scalar x, Scalar y) {
		if constexpr (std::is_integral_v<Scalar>) {
			return x % y;
		}
		else {
			return std::fmod(x, y);
		}
	}

	template<typename Scalar>
	inline Scalar customfloor(Scalar x) {
		if constexpr (std::is_integral_v<Scalar>) {
			return x;
		}
		else {
			return floor(x);
		}
	}

}

namespace Mackey{
	/// The Smith Normal Form
	template <typename smithS_t, typename smithR_t, typename smithC_t>
	class Smith {
	public:
		smithS_t S;		///< The Smith Normal Form.
		smithR_t P;		///< One of the coefficient matrices (S=P*A*Q). Should be row major for best performance
		smithR_t Qi;	///< One of the coefficient matrices (A=Pi*S*Qi). Should be row major for best performance
		smithC_t Q;		///< One of the coefficient matrices (S=P*A*Q). Should be column major for best performance
		smithC_t Pi;	///< One of the coefficient matrices (A=Pi*S*Qi). Should be column major for best performance
		const int L;	///< The length of the diagonal of S.


	/////////////////////////////////////////////////
	/// Given a matrix A compute a diagonal matrix S and invertible matrices P,Q,Pi,Qi such that S=P*A*Q and A=Pi*S*Qi

	///The Pi,Qi are the inverses of P,Q. If wantP=1 then P,Pi are computed and if wantQ=1 then Q,Qi are computed.
	/////////////////////////////////////////////////
		Smith(const smithS_t& Original, const bool& wantP, const bool& wantQ);

	private:
		const int M, N;
		const bool wantP, wantQ;
		void SmithIt();
		void workRow(const int&, int&, int&);
		void searchRow(const int&, int&, int&);
		void workCol(const int&, int&, int&);
		void searchCol(const int&, int&, int&);
	};


	template <typename smithS_t, typename smithR_t, typename smithC_t>
	Smith<smithS_t, smithR_t, smithC_t> ::Smith(const smithS_t& Original, const bool& wantP, const bool& wantQ)
		: S(Original), M(S.rows()), N(S.cols()), L(std::min(S.rows(), S.cols())), wantP(wantP), wantQ(wantQ) { //Initializer List
		if (wantP) {
			P.resize(M, M);
			P.setIdentity();
			Pi = P;
		}
		if (wantQ) {
			Q.resize(N, N);
			Q.setIdentity();
			Qi = Q;
		}
		SmithIt();
	}


	template <typename smithS_t, typename smithR_t, typename smithC_t>
	void Smith<smithS_t, smithR_t, smithC_t>::SmithIt() {
		for (int start = 0; start < L; start++) {
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
		}
	}

	template <typename smithS_t, typename smithR_t, typename smithC_t>
	void Smith<smithS_t, smithR_t, smithC_t>::workRow(const int& start, int& i, int& j) {
		if (mod(S(start, j), S(start, start)) == 0) {
			auto thequotient = S(start, j) / S(start, start);
			S.col(j).tail(M - start) += -thequotient * S.col(start).tail(M - start);
			if (wantQ) {
				Q.col(j) += -thequotient * Q.col(start);
				Qi.row(start) += thequotient * Qi.row(j);
			}
			j++;
		}
		else if (mod(S(start, start), S(start, j)) == 0) {
			S.col(start).tail(M - start).swap(S.col(j).tail(M - start));
			auto thequotient = S(start, j) / S(start, start);
			S.col(j).tail(M - start) += -thequotient * S.col(start).tail(M - start);
			if (wantQ) {
				Q.col(start).swap(Q.col(j));
				Q.col(j) += -thequotient * Q.col(start);
				Qi.row(start).swap(Qi.row(j));
				Qi.row(start) += thequotient * Qi.row(j);
			}
			j++;
			i = start + 1;
		}
		else {
			auto thequotient = customfloor(S(start, j) / S(start, start));
			S.col(j).tail(M - start) += -thequotient * S.col(start).tail(M - start);
			S.col(start).swap(S.col(j));
			if (wantQ) {
				Q.col(j) += -thequotient * Q.col(start);
				Q.col(start).swap(Q.col(j));
				Qi.row(start) += thequotient * Qi.row(j);
				Qi.row(start).swap(Qi.row(j));
			}
			i = start + 1;
		}
	}

	template <typename smithS_t, typename smithR_t, typename smithC_t>
	void Smith<smithS_t, smithR_t, smithC_t>::workCol(const int& start, int& i, int& j) {
		if (mod(S(i, start), S(start, start)) == 0) {
			auto thequotient = S(i, start) / S(start, start);
			S.row(i).tail(N - start) += -thequotient * S.row(start).tail(N - start);
			if (wantP) {
				P.row(i) += -thequotient * P.row(start);
				Pi.col(start) += thequotient * Pi.col(i);
			}
			i++;
		}
		else if (mod(S(start, start), S(i, start)) == 0) {
			S.row(start).tail(N - start).swap(S.row(i).tail(N - start));
			auto thequotient = S(i, start) / S(start, start);
			S.row(i).tail(N - start) += -thequotient * S.row(start).tail(N - start);
			if (wantP) {
				P.row(start).swap(P.row(i));
				P.row(i) += -thequotient * P.row(start);
				Pi.col(start).swap(Pi.col(i));
				Pi.col(start) += thequotient * Pi.col(i);
			}
			i++;
			j = start + 1;
		}
		else {
			auto thequotient = customfloor(S(i, start) / S(start, start));
			S.row(i).tail(N - start) += -thequotient * S.row(start).tail(N - start);
			S.row(start).tail(N - start).swap(S.row(i).tail(N - start));
			if (wantP) {
				P.row(i) += -thequotient * P.row(start);
				P.row(start).swap(P.row(i));
				Pi.col(start) += thequotient * Pi.col(i);
				Pi.col(start).swap(Pi.col(i));
			}
			j = start + 1;
		}
	}




	template <typename smithS_t, typename smithR_t, typename smithC_t>
	void Smith<smithS_t, smithR_t, smithC_t>::searchRow(const int& start, int& i, int& j)
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
			if (found == 0) {
				j = N;
			}
		}
		if (j < N && S(start, j) == 0)
		{
			for (int k = j + 1; k < N; k++)
			{
				if (S(start, k) != 0)
				{
					j = k;
					return;
				}
			}
			//If you made it here then row is zero outside of start
			j = N;
		}
	}


	template <typename smithS_t, typename smithR_t, typename smithC_t>
	void Smith<smithS_t, smithR_t, smithC_t>::searchCol(const int& start, int& i, int& j)
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
			if (found == 0) {
				i = M;
			}
		}
		if (i < M && S(i, start) == 0) {
			for (int k = i + 1; k < M; k++) {
				if (S(k, start) != 0) {
					i = k;
					return;
				}
			}
			//If you made it here then column is zero outside of start
			i = M;
		}
	}
}
