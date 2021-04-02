#pragma once
#include "../MackeyFunctor.hpp"
///@file
///@brief Contains the class describing Mackey Functors.

namespace mackey
{
	namespace implementation_details {

		template<typename Derived>
		inline bool genericcompare(const Eigen::MatrixBase<Derived>& A, const Eigen::MatrixBase<Derived>& B) {
			if (A.size() == 0 && B.size() == 1 && B(0, 0) == 0)
				return 1;
			if (B.size() == 0 && A.size() == 1 && A(0, 0) == 0)
				return 1;
			if (A.size() != B.size() || A != B)
				return 0;
			return 1;
		}

		template<typename T>
		inline bool genericcompare(const std::vector<T>& A, const std::vector<T>& B) {
			if (A.size() != B.size())
				return 0;
			else
				for (size_t i = 0; i < A.size(); i++)
					if (!genericcompare(A[i], B[i]))
						return 0;
			return 1;
		}

		///The sum of two linear maps. We account for a map being empty if domain/range is 0 hence why the dimensions of domains and ranges are given.
		template<typename Matrix_t>
		Matrix_t blockdiag(const Matrix_t& A, const Matrix_t& B, int domA, int ranA, int domB, int ranB) {
			if (A.size() == 0 && B.size() != 0) {
				Matrix_t C = Eigen::MatrixBase<Matrix_t>::Zero(ranA + B.rows(), domA + B.cols());
				C.block(ranA, domA, B.rows(), B.cols()) = B;
				return C;
			}
			else if (B.size() == 0 && A.size() != 0) {
				Matrix_t C = Eigen::MatrixBase<Matrix_t>::Zero(A.rows() + ranB, A.cols() + domB);
				C.block(0, 0, A.rows(), A.cols()) = A;
				return C;
			}
			Matrix_t C = Eigen::MatrixBase<Matrix_t>::Zero(A.rows() + B.rows(), A.cols() + B.cols());
			C.block(0, 0, A.rows(), A.cols()) = A;
			C.block(A.rows(), A.cols(), B.rows(), B.cols()) = B;
			return C;
		}

		///The sum of two linear maps.
		template<typename Matrix_t>
		Matrix_t blockdiag(const Matrix_t& A, const Matrix_t& B) {
			return blockdiag(A, B, A.cols(), A.rows(), B.cols(), B.rows());
		}

		///Given matrix and change of basis matrices, produces the matrix w.r.t. the new basis and normalizes it for the given group in the image
		template<typename T, typename S>
		T change_base_matrix_normalize(const T& mat, const T& P, const T& Q, const S& Group_image) {
			if (mat.size() == 0)
				return T();
			//Cast for extended precision in matrix multiplication before normalizing. Usefull as T may well have Scalar type char
			Eigen::Matrix<float, -1, -1> a, p, q, d;
			Eigen::Matrix<int, -1, -1> e;
			a = mat.template cast<float>();
			p = P.template cast<float>();
			q = Q.template cast<float>();
			d = p * a * q;
			e = d.template cast<int>();
			Group_image.normalize(e);
			//cast back
			return e.template cast<typename T::Scalar>();
		}
	}


	//Resize all member variables.
	template<typename rank_t>
	void MackeyFunctor<rank_t>::resize(int levels) {
		Groups.resize(levels); Tr.resize(levels - 1); Res.resize(levels - 1); Weyl.resize(levels - 1);
	}


	//Compares two Mackey and returns 1 if equal.
	template<typename rank_t>
	bool MackeyFunctor<rank_t>::operator==(const MackeyFunctor<rank_t>& M) const {
		if (!name.empty() && !M.name.empty())
			return (name == M.name);
		else
			return (Groups == M.Groups && implementation_details::genericcompare(Tr, M.Tr) && implementation_details::genericcompare(Res, M.Res) && implementation_details::genericcompare(Weyl, M.Weyl));
	}

	template<typename rank_t>
	void MackeyFunctor<rank_t>::notation()
	{
		if (!name.empty())
			return;
		size_t number_of_zeros = 0;
		for (const auto& i : Groups) {
			if (i.istrivial()) {
				name.append("0");
				number_of_zeros++;
			}
			else if (i.iscyclic())
				name.append(std::to_string(i[0]));
			else {
				name.clear();
				return;
			}
		}
		if (number_of_zeros == Groups.size()) {
			name = "0";
			return;
		}
		bool sharp = 1;
		for (size_t i = 0; i < Tr.size(); i++) {
			if (!Groups[i].istrivial() && !Groups[i + 1].istrivial()) {
				if ((Tr[i](0, 0) == 1 || Tr[i](0, 0) * Groups[i][0] == Groups[i + 1][0]) && ((Groups[i][0] == 1 && Res[i](0, 0) == 2) || ((Res[i](0, 0) - 2) % Groups[i][0] == 0))) {
					if (sharp)
						name.append(" # ");
					name.append(std::to_string(i));
					sharp = 0;
				}
				else if (!((Res[i](0, 0) == 1 || Res[i](0, 0) * Groups[i + 1][0] == Groups[i][0]) && ((Groups[i + 1][0] == 1 && Tr[i](0, 0) == 2) || ((Tr[i](0, 0) - 2) % Groups[i + 1][0] == 0)))) {
					name.clear();
					return;
				}

			}
		}
	}

	template<typename rank_t>
	MackeyFunctor<rank_t>::MackeyFunctor(const rank_t& bigarray, int length) {
		resize(length);
		for (int i = 0; i < length; i++) {
			if (bigarray(i) != 0)
				Groups[i] = bigarray.segment(i, 1);
		}
		for (int i = 0; i < length - 1; i++) {
			if (bigarray(i) != 0 && bigarray(i + 1) != 0) {
				Tr[i] = bigarray.segment(i + length, 1);
				Res[i] = bigarray.segment(i + 2 * length - 1, 1);
			}
			if (bigarray(i) != 0)
				Weyl[i] = bigarray.segment(i + 3 * length - 2, 1);
		}
	}

	///The sum of two Mackey functors
	template<typename rank_t>
	MackeyFunctor<rank_t> operator +(const MackeyFunctor<rank_t>& M, const MackeyFunctor<rank_t>& N) {
		MackeyFunctor<rank_t> K;
		K.resize(M.Groups.size());
		for (size_t i = 0; i < M.Groups.size(); i++) {
			K.Groups[i] = M.Groups[i] + N.Groups[i];
		}
		for (size_t i = 0; i < M.Groups.size() - 1; i++) {
			K.Tr[i] = implementation_details::blockdiag(M.Tr[i], N.Tr[i], M.Groups[i].number_of_summands(), M.Groups[i + 1].number_of_summands(), N.Groups[i].number_of_summands(), N.Groups[i + 1].number_of_summands());
			K.Res[i] = implementation_details::blockdiag(M.Res[i], N.Res[i], M.Groups[i + 1].number_of_summands(), M.Groups[i].number_of_summands(), N.Groups[i + 1].number_of_summands(), N.Groups[i].number_of_summands());
			K.Weyl[i] = implementation_details::blockdiag(M.Weyl[i], N.Weyl[i]);
		}
		if (!M.name.empty() && !N.name.empty()) {
			if (M.name == "0")
				K.name = N.name;
			else if (N.name == "0")
				K.name = M.name;
			else
				K.name = M.name + " + " + N.name;
		}
		return K;
	}



	//Printing a Mackey functor
	template<typename rank_t>
	std::string MackeyFunctor<rank_t>::print() const {
		std::stringstream out;
		out << "Groups:";
		for (const auto& i : Groups) {
			if (i.istrivial())
				out << "0,";
			else
				out << i << ",";
		}
		out << "\n Tr:";
		for (const auto& i : Tr) {
			if (i.size() == 0)
				out << "0 \n";
			else
				out << i << "\n";
		}
		out << "\n Res:";
		for (const auto& i : Res) {
			if (i.size() == 0)
				out << "0 \n";
			else
				out << i << "\n";
		}
		out << "\n Weyl:";
		for (const auto& i : Weyl) {
			if (i.size() == 0)
				out << "0 \n";
			else
				out << i << "\n";
		}
		return out.str();
	}


	//Print Mackey Functor to output stream
	template<typename rank_t>
	std::ostream& operator<<(std::ostream& out, const MackeyFunctor<rank_t>& M) {
		if (!M.name.empty())
			out << M.name;
		else
			out << M.print();
		return out;
	}



	template<typename rank_t>
	MackeyFunctor<rank_t> MackeyFunctor<rank_t>::apply(const iso_t<rank_t>& iso, const iso_t<rank_t>& isoinv) const {
		MackeyFunctor<rank_t> N = *this;
		for (size_t i = 0; i < Tr.size(); i++)
			N.Tr[i] = implementation_details::change_base_matrix_normalize(Tr[i], iso[i + 1], isoinv[i], Groups[i + 1]);
		for (size_t i = 0; i < Res.size(); i++)
			N.Res[i] = implementation_details::change_base_matrix_normalize(Res[i], iso[i], isoinv[i + 1], Groups[i]);
		for (size_t i = 0; i < Weyl.size(); i++)
			N.Weyl[i] = implementation_details::change_base_matrix_normalize(Weyl[i], iso[i], isoinv[i], Groups[i]);
		return N;
	}


	template<typename rank_t>
	std::pair<std::vector<iso_t<rank_t>>, std::vector<iso_t<rank_t>>> MackeyFunctor<rank_t>::automorphisms() const {
		std::vector<std::vector<dense_t<rank_t>>> isosfirst, isossecond;
		isosfirst.reserve(Groups.size());
		isossecond.reserve(Groups.size());
		for (const auto& i : Groups) {
			auto pair = i.all_automorphisms();
			isosfirst.push_back(pair.first);
			isossecond.push_back(pair.second);
		}
		auto comb = combinations(isosfirst);
		auto combinv = combinations(isossecond);
		return std::make_pair(comb, combinv);
	}

	template<typename rank_t>
	std::vector<MackeyFunctor<rank_t>> MackeyFunctor<rank_t>::isomorphism_class() const {
		auto autom = automorphisms();
		std::vector<MackeyFunctor<rank_t>> result;
		result.reserve(autom.first.size());
		for (size_t i = 0; i < autom.first.size(); i++) {
			result.push_back(apply(autom.first[i], autom.second[i]));
		}
		return result;
	}

	template<typename rank_t>
	bool MackeyFunctor<rank_t>::isomorphic(const std::vector<MackeyFunctor<rank_t>>& isoclass) const {
		if (isoclass[0].Groups != Groups)
			return 0;
		for (const auto& i : isoclass) {
			if (*this == i)
				return 1;
		}
		return 0;
	}

	template<typename rank_t>
	bool MackeyFunctor<rank_t>::isomorphic(const MackeyFunctor<rank_t>& M) {
		if (Groups != M.Groups)
			return 0;
		auto isoclass = isomorphism_class();
		return isomorphic(isoclass);
	}
}
