#pragma once
#include "../Box.hpp"

///@file
///@brief Contains the construction of box products of Chains.


namespace mackey {

	//Returns the tensor product
	template<typename _output_t, typename _optional_t>
	_output_t& Tensor<_output_t, _optional_t>::tensor() {
		return _tensor;
	}

	//Returns the optional parameter aka the detailed rank
	template<typename _output_t, typename _optional_t>
	auto& Tensor<_output_t, _optional_t>::optional() {
		static_assert(optional_exists, "You need to provide a non void type for optional in order for it to be computed");
		return this->_optional;
	}

	//Returns the tensor product
	template<typename _output_t, typename _optional_t>
	const _output_t& Tensor<_output_t, _optional_t>::tensor() const {
        return _tensor;
	}

	//Returns the optional parameter aka the detailed rank
	template<typename _output_t, typename _optional_t>
	const auto& Tensor<_output_t, _optional_t>::optional() const {
		static_assert(optional_exists, "You need to provide a non void type for optional in order for it to be computed");
		return this->_optional;
	}

	//Constructor tensors A,B optionally at given location only
	template<typename _output_t, typename _optional_t>
	template<typename _input_t>
	Tensor<_output_t, _optional_t>::Tensor(const _input_t& A, const _input_t& B, int location) {
		tensor(A, B, location);
	}


	//_input_t is rank_t, _output_t is rank_t, _optional_t is void
	template<typename _output_t, typename _optional_t>
	template<typename _input_t, typename implementation_details::is_rank<_input_t>>
	void Tensor<_output_t, _optional_t>::tensor(const _input_t& A, const _input_t& B, int) {
		int sum = 0;
		for (int i = 0; i < A.size(); i++) {
			for (int j = 0; j < B.size(); j++) {
				sum += std::min(A[i], B[j]);
			}
		}
		_tensor = _output_t(sum);
		int current = 0;
		for (int j = 0; j < B.size(); j++) {
			for (int i = 0; i < A.size(); i++) {
				_tensor.segment(current, std::min(A[i], B[j])).setConstant(std::max(A[i], B[j]));
				current += std::min(A[i], B[j]);
			}
		}
	}

	//_input_t is vector<rank_t>, _output_t is rank_t, _optional_t can be std::vector<int64_t>
	template<typename _output_t, typename _optional_t>
	template<typename _input_t, typename implementation_details::is_vector_rank<_input_t>, typename S, typename implementation_details::is_rank<S>>
	void Tensor<_output_t, _optional_t>::tensor(const _input_t& A, const _input_t& B, int location_of_rank) {
		_input_t verydetailedrank;
		if constexpr (optional_exists)
			this->_optional.reserve(location_of_rank + 1);
		verydetailedrank.resize(location_of_rank + 1);
		auto lowlimit = std::max(0, location_of_rank - (int)B.size() + 1);
		auto highlimit = std::min(location_of_rank, (int)A.size() - 1);
		int64_t sum = 0;
		for (int j = lowlimit; j <= highlimit; j++) {
			verydetailedrank[j] = Tensor<_output_t>(A[j], B[location_of_rank - j]).tensor();
			sum += verydetailedrank[j].size();
		}
		_tensor.resize(sum);
		int64_t track = 0;
		for (const auto& i : verydetailedrank) {
			_tensor.segment(track, i.size()) = i;
			if constexpr (optional_exists)
				this->_optional.push_back(summation(i));
			track += i.size();
		}
	}

	//_input_t is vector<rank_t>, _output_t is vector<rank_t>, _optional_t is void
	template<typename _output_t, typename _optional_t>
	template<typename _input_t, typename T, typename S, typename implementation_details::is_vector_rank<T>, typename implementation_details::is_vector_rank<S>>
	void Tensor<_output_t, _optional_t>::tensor(const _input_t& A, const _input_t& B, int) {
		_tensor.reserve(A.size() + B.size() - 1);
		for (int i = 0; i < A.size() + B.size() - 1; i++)
			_tensor.push_back(Tensor<typename _output_t::value_type>(A, B, i).tensor());
	}

	//Computing the differential of the tensor product
	namespace implementation_details {

		///Forms the block diagonal matrix with k many copies of a, and then applies the permutations left^{-1} * block *right 
		///Significant performance improvement over the naive implementation by forming left^{-1}*block in one step
		template<typename Derived, typename S>
		Derived permutation_block(const S& left, const S& right, const Eigen::MatrixBase<Derived>& a, int k) {
			int n = a.cols();
			int m = a.rows();
			Derived b = Eigen::MatrixBase<Derived>::Zero(m * k, n * k);
			for (int i = 0; i < b.rows(); i++) {
				auto u = left[i] % m;
				auto v = left[i] / m;
				b.row(i).segment(n * v, n) = a.row(u).segment(0, n);
			}
			//now apply the permutation using the Eigen API
			typedef typename S::value_type pScalar;
			return b * Eigen::PermutationMatrix<-1, -1, pScalar>(Eigen::Map<const Eigen::Matrix<pScalar, -1, 1>>(right.data(), right.size()));
		}



		///Mixes the left and right differentials to produce the differential for the box product
		///The left differentials are stored in a vector L and the right differentials in a vector R.
		///The mixing is L[0], then R[0] directly to the right, then L[1] directly below then...
		///It's assumed that L,R have consistent dimensions for the blocks to fit coherently.
		template<typename T>
		T MatrixMixer(std::vector<T>& L, std::vector<T>& R, int rows, int cols) {
			T mixed;
			mixed.setZero(rows, cols);
			int horz = 0, vert = 0;
			for (size_t i = 0;i < L.size();i++) {
				auto rowsL = L[i].rows();
				auto colsL = L[i].cols();
				if (rowsL > 0) {
					mixed.block(vert, horz, rowsL, colsL) = L[i];
					vert += rowsL;
				}
				auto rowsR = R[i].rows();
				auto colsR = R[i].cols();

				if (colsR > 0) {
					mixed.block(vert, horz, rowsR, colsR) = R[i];
					horz += colsR;
				}
			}
			return mixed;
		}


		///Does permutation_block (using right^{-1} !!!) and matrix mixing in one step for sparse matrices, using the triplets format.
		template<typename T, typename S, typename storage, typename U, typename V>
		void box_sparse(T& trip, const S& left, const S& right, const T& a, U rows, U cols, V copies, storage h_offset, storage v_offset) {
			for (size_t k = 0; k < a.size(); k++) {
				auto i = a[k].row();
				auto j = a[k].col();
				while (i < rows * copies && j < cols * copies) {
					trip.emplace_back(left[i] + v_offset, right[j] + h_offset, a[k].value());
					i += rows;
					j += cols;
				}
			}
		}

		///Computes the differential (A_* tensor B_*)_i -> (A_* tensor B_*)_{i-1}. Dense implementation. Please note that we use the Canon_to_Convenient perms here due to how permutation_block works
		template<typename rank_t, typename diff_t, typename std::enable_if_t<mackey::SFINAE::is_Dense<diff_t>::value, int> = 0>
		diff_t diff_of_Tensor(const Chains<rank_t, diff_t>& A, const Chains<rank_t, diff_t>& B, int i, int64_t rows, int64_t cols) {
			auto lowlimit = std::max(0, i - B.maxindex());
			auto highlimit = std::min(i, A.maxindex());
			std::vector<diff_t> LeftDiff(i + 1);
			std::vector<diff_t> RightDiff(i + 1);
			for (int j = lowlimit; j <= highlimit; j++) {
				if (A.rank[j].size() == 0 || B.rank[i - j].size() == 0)
					continue;
				auto Domain = ChangeBasis<storage_t<diff_t>>(A.rank[j], B.rank[i - j], 0, 1);
				if (j >= 1) //We have a LeftDiff
				{
					auto RangeL = ChangeBasis<storage_t<diff_t>>(A.rank[j - 1], B.rank[i - j], 0, 1);
					//LeftDiff= R_CO_CA * CO_CO * D_CA_CO and permutation_block does ran^{-1} * diff * dom hence why we provide Range_Canon_to_Left
					LeftDiff[j] = permutation_block(RangeL.canon_to_left, Domain.canon_to_left, A.diff[j], summation(B.rank[i - j]));
				}
				if (i - j >= 1) //We have a RightDiff
				{
					scalar_t<diff_t> sign(1 - 2 * (j % 2)); //(1 - 2 * (j % 2)) is (-1) ^ j
					auto RangeR = ChangeBasis<storage_t<diff_t>>(A.rank[j], B.rank[i - j - 1], 0, 1);
					RightDiff[j] = permutation_block(RangeR.canon_to_right, Domain.canon_to_right, static_cast<diff_t>(sign * B.diff[i - j]), summation(A.rank[j]));
				}
			}
			return MatrixMixer(LeftDiff, RightDiff, rows, cols);
		}

		///Computes the differential (A_* tensor B_*)_i -> (A_* tensor B_*)_{i-1}. Sparse implementation. Please note that we use the Conv_to_Canon perms here due to how box_sparse works
		template<typename rank_t, typename diff_t, typename std::enable_if_t<mackey::SFINAE::is_Sparse<diff_t>::value, int> = 0>
		diff_t diff_of_Tensor(const Chains<rank_t, diff_t>& A, const Chains<rank_t, diff_t>& B, int i, int64_t rows, int64_t cols) {
			auto lowlimit = std::max(0, i - B.maxindex());
			auto highlimit = std::min(i, A.maxindex());
			storage_t<diff_t> totalsize, v_offset, h_offset;
			totalsize = v_offset = h_offset = 0;

			//first find the size of nonzeros
			for (int j = lowlimit; j <= highlimit; j++) {
				if (j >= 1)
					totalsize += A.diff[j].nonZeros() * summation(B.rank[i - j]);
				if (i - j >= 1)
					totalsize += B.diff[i - j].nonZeros() * summation(A.rank[j]);
			}

			std::vector<triplet_t<diff_t>> mat;
			mat.reserve(totalsize);

			//put the nonzeros in mat
			for (int j = lowlimit; j <= highlimit; j++) {
				if (A.rank[j].size() == 0 || B.rank[i - j].size() == 0)
					continue;
				auto Domain = ChangeBasis<storage_t<diff_t>>(A.rank[j], B.rank[i - j], 1, 0);
				if (j >= 1 && A.rank[j - 1].size() != 0) //We have a LeftDiff
				{
					auto RangeL = ChangeBasis<storage_t<diff_t>>(A.rank[j - 1], B.rank[i - j], 1, 0);
					//LeftDiff= R_CO_CA * CO_CO * D_CA_CO and box_sparse does ran * diff * dom^{-1} hence why we provide Domain.left_to_canon
					box_sparse(mat, RangeL.left_to_canon, Domain.left_to_canon, make_triplets(A.diff[j]), A.diff[j].rows(), A.diff[j].cols(), summation(B.rank[i - j]), h_offset, v_offset);
					v_offset += A.diff[j].rows() * summation(B.rank[i - j]);
				}
				if (i - j >= 1 && B.rank[i - j - 1].size() != 0) //We have a RightDiff
				{
					scalar_t<diff_t> sign(1 - 2 * (j % 2)); //(1 - 2 * (j % 2)) is(-1) ^ j
					auto RangeR = ChangeBasis<storage_t<diff_t>>(A.rank[j], B.rank[i - j - 1], 1, 0);
					box_sparse(mat, RangeR.right_to_canon, Domain.right_to_canon, make_triplets(static_cast<diff_t>(sign * B.diff[i - j])), B.diff[i - j].rows(), B.diff[i - j].cols(), summation(A.rank[j]), h_offset, v_offset);
					h_offset += B.diff[i - j].cols() * summation(A.rank[j]);
				}
			}
			//turn mat to a matrix
			diff_t diff;
			diff.resize(rows, cols);
			diff.setFromTriplets(mat.begin(), mat.end());
			return diff;
		}


	}

	template<typename _output_t, typename _optional_t>
	template<typename _input_t, typename T, typename implementation_details::is_arrow<T>>
	void Tensor<_output_t, _optional_t>::tensor(const _input_t& A, const _input_t& B, int i) {
		typedef typename T::rank_t rank_t;
		typedef typename T::diff_t diff_t;
		auto Ranks = Tensor<rank_t, _optional_t>(A.rank, B.rank, i);
		auto rank_domain = Ranks.tensor();
		rank_t rank_range;
		if (i > 0)
			rank_range = Tensor<rank_t, _optional_t>(A.rank, B.rank, i - 1).tensor();
		diff_t diff;
		if (rank_domain.size() != 0 && rank_range.size() != 0)
			diff = implementation_details::diff_of_Tensor<rank_t, diff_t>(A, B, i, summation(rank_range), summation(rank_domain));
		_tensor = Arrow<rank_t, diff_t>(rank_domain, rank_range, diff);
		if constexpr (optional_exists)
			this->_optional = Ranks.optional();
	}

	template<typename _output_t, typename _optional_t>
	template<typename _input_t, typename T, typename implementation_details::is_chains<T>>
	void Tensor<_output_t, _optional_t>::tensor(const _input_t& A, const _input_t& B, int i) {
		if (i == -1)
			i = A.maxindex() + B.maxindex();
		else
			i = std::min(i, A.maxindex() + B.maxindex());
		_tensor.reserve(i + 1);
		if constexpr (optional_exists) {
			this->_optional.reserve(i + 1);
			for (int j = 0; j <= i; j++) {
				Tensor<arrow_t<T>, typename _optional_t::value_type> BoxedAtPointj(A, B, j);
				_tensor.push_back(BoxedAtPointj.tensor());
				this->_optional.push_back(BoxedAtPointj.optional());
			}
		}
		else
			for (int j = 0; j <= i;j++)
				_tensor.push_back(Tensor<arrow_t<T>>(A, B, j).tensor());
	}


	template<typename _output_t, typename _optional_t>
	template<typename _input_t, typename T, typename implementation_details::is_junction<T>>
	void Tensor<_output_t, _optional_t>::tensor(const _input_t& A, const _input_t& B, int i) {
		Tensor<arrow_t<T>, _optional_t> OUT(A, B, i);
		_tensor.setOut(OUT.tensor());
		if (i < A.maxindex() + B.maxindex())
			_tensor.setIn(Tensor<arrow_t<T>, _optional_t>(A, B, i + 1).tensor(), 0);
		if constexpr (optional_exists)
			this->_optional = OUT.optional();
	}
}
