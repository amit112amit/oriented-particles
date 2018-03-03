#include "HelperFunctions.h"

namespace OPS{
    // Given two sets of 3D points, find the rotation + translation + scale
    // which best maps the first set to the second.
    // Source: http://en.wikipedia.org/wiki/Kabsch_algorithm

    // The input 3D points are stored as columns.
    Eigen::Affine3d find3DAffineTransform(Eigen::Ref<Eigen::Matrix3Xd> in,
	    Eigen::Ref<Eigen::Matrix3Xd> out){
	// Default output
	Eigen::Affine3d A;
	A.linear() = Eigen::Matrix3d::Identity(3, 3);
	A.translation() = Eigen::Vector3d::Zero();

	if (in.cols() != out.cols())
	    throw "Find3DAffineTransform(): input data mis-match";

	// Find the centroids then shift to the origin
	Eigen::Vector3d in_ctr = Eigen::Vector3d::Zero();
	Eigen::Vector3d out_ctr = Eigen::Vector3d::Zero();
	for (int col = 0; col < in.cols(); col++) {
	    in_ctr += in.col(col);
	    out_ctr += out.col(col);
	}
	in_ctr /= in.cols();
	out_ctr /= out.cols();
	for (int col = 0; col < in.cols(); col++) {
	    in.col(col) -= in_ctr;
	    out.col(col) -= out_ctr;
	}
	Eigen::MatrixXd Cov = in * out.transpose();
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov,
		Eigen::ComputeThinU | Eigen::ComputeThinV);

	// Find the rotation
	double_t d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
	if (d > 0)
	    d = 1.0;
	else
	    d = -1.0;
	Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
	I(2, 2) = d;
	Eigen::Matrix3d R = svd.matrixV() * I * svd.matrixU().transpose();

	// The final transform
	A.linear() = R;
	A.translation() = out_ctr - R * in_ctr;

	return A;
    }

    // A function to test Find3DAffineTransform()

    void TestFind3DAffineTransform() {

	// Create datasets with known transform
	Eigen::Matrix3Xd in(3, 100), out(3, 100);
	Eigen::Quaternion<double> Q(1, 3, 5, 2);
	Q.normalize();
	Eigen::Matrix3d R = Q.toRotationMatrix();
	for (int row = 0; row < in.rows(); row++) {
	    for (int col = 0; col < in.cols(); col++) {
		in(row, col) = log(2.0 * row + 10.0) / sqrt(1.0 * col + 4.0)
		    + sqrt(col * 1.0) / (row + 1.0);
	    }
	}
	Eigen::Vector3d S;
	S << -5, 6, -27;
	for (int col = 0; col < in.cols(); col++)
	    out.col(col) = R * in.col(col) + S;

	Eigen::Affine3d A = find3DAffineTransform(in, out);

	// See if we got the transform we expected
	if ((R - A.linear()).cwiseAbs().maxCoeff() > 1e-13
		|| (S - A.translation()).cwiseAbs().maxCoeff() > 1e-13){
	    //throw "Could not determine the affine transform accurately enough";
	}
    }

    // Find 2d affine transformation
    Eigen::Affine2d find2DAffineTransform(Eigen::Ref<Eigen::Matrix2Xd> in,
	    Eigen::Ref<Eigen::Matrix2Xd> out){
	// Default output
	Eigen::Affine2d A;
	A.linear() = Eigen::Matrix2d::Identity(2, 2);
	A.translation() = Eigen::Vector2d::Zero();

	if (in.cols() != out.cols())
	    throw "Find2DAffineTransform(): input data mis-match";

	// Find the centroids then shift to the origin
	Eigen::Vector2d in_ctr = Eigen::Vector2d::Zero();
	Eigen::Vector2d out_ctr = Eigen::Vector2d::Zero();
	for (int col = 0; col < in.cols(); col++) {
	    in_ctr += in.col(col);
	    out_ctr += out.col(col);
	}
	in_ctr /= in.cols();
	out_ctr /= out.cols();
	for (int col = 0; col < in.cols(); col++) {
	    in.col(col) -= in_ctr;
	    out.col(col) -= out_ctr;
	}

	// SVD
	Eigen::MatrixXd Cov = in * out.transpose();
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(Cov,
		Eigen::ComputeThinU | Eigen::ComputeThinV);

	// Find the rotation
	double_t d = (svd.matrixV() * svd.matrixU().transpose()).determinant();
	if (d > 0)
	    d = 1.0;
	else
	    d = -1.0;
	Eigen::Matrix2d I = Eigen::Matrix2d::Identity(2, 2);
	I(1, 1) = d;
	Eigen::Matrix2d R = svd.matrixV() * I * svd.matrixU().transpose();

	// The final transform
	A.linear() = R;
	A.translation() = out_ctr - R * in_ctr;

	return A;
    }

    // A function to test Find3DAffineTransform()
    void TestFind2DAffineTransform() {
	// Create datasets with known transform
	Eigen::Matrix2Xd in(2, 100), out(2, 100);
	Eigen::Matrix2d R;
	R << 0.70761, -0.70761, 0.70761, 0.70761;
	Eigen::Vector2d S;
	S << -5,6;
	for (int row = 0; row < in.rows(); row++) {
	    for (int col = 0; col < in.cols(); col++) {
		in(row, col) = log(2.0 * row + 10.0) / sqrt(1.0 * col + 4.0)
		    + sqrt(col * 1.0) / (row + 1.0);
	    }
	}
	for (int col = 0; col < in.cols(); col++)
	    out.col(col) = R * in.col(col) + S;
	Eigen::Affine2d A = find2DAffineTransform(in, out);
	// See if we got the transform we expected
	if ((R - A.linear()).cwiseAbs().maxCoeff() > 1e-13
		|| (S - A.translation()).cwiseAbs().maxCoeff() > 1e-13){
	    throw "Could not determine the affine transform accurately enough";
	}
    }

}
