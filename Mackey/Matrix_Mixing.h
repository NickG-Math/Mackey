#pragma once
#include <Eigen/Dense>
#include <vector>

///@file
///@brief Contains block operations on matrices.


///Forms the block diagonal matrix with n blocks all equal to A.
template<typename Derived>
Derived blkdiag(const Eigen::MatrixBase<Derived>& A, const int n) {
	Derived B=Eigen::MatrixBase<Derived>::Zero(n * A.rows(), n * A.cols());
	for (int i = 0;i < n;i++) {
		B.block(i * A.rows(), i * A.cols(), A.rows(), A.cols()) = A;
	}
	return B;
}

/////////////////////////////////////////////////////////////////////////////
///Mixes the left and right differentials to produce the differential for the box product

///The left differentials are stored in a vector L and the right differentials in a vector R.
///The mixing is L[0], then R[0] directly to the right, then L[1] directly below then...
///It's assumed that L,R have consistent dimensions for the blocks to fit coherently.
////////////////////////////////////////////////////////////////////////////
template<typename Derived>
Derived MatrixMixer(std::vector<Derived>& L, std::vector<Derived>& R) {
	auto rows = L[0].rows();
	auto cols = std::max(L[0].cols(),R[0].cols());
	for (size_t i = 0; i < L.size()-1; i++) {
		rows += std::max(L[i+1].rows(),R[i].rows());
		cols += std::max(L[i + 1].cols(), R[i+1].cols());
	}
	rows += R[R.size()-1].rows();

	Derived mixed;
	mixed.setZero(rows, cols);
	int horz = 0, vert=0;
	for (size_t i = 0;i < L.size();i++) {
		auto rowsL = L[i].rows();
		auto colsL = L[i].cols();
		if (rowsL > 0) {
			mixed.block(vert, horz, rowsL, colsL) =std::move(L[i]);
			vert+= rowsL;
		}
		auto rowsR = R[i].rows();
		auto colsR = R[i].cols();

		if (colsR > 0) {
			mixed.block(vert, horz, rowsR, colsR) = std::move(R[i]);
			horz+= colsR;
		}
	}
	return mixed;
}
