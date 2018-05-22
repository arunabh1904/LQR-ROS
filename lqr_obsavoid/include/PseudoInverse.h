/*
 * PseudoInverse.h
 *
 */

#ifndef PSEUDOINVERSE_HPP_
#define PSEUDOINVERSE_HPP_

#include <eigen3/Eigen/SVD>

namespace Math {
    // Eigen Pseudo Inverse Using SVD and setting un-observable (as defined by the singular min to max ratio) singular values to inf
    Eigen::VectorXd pseudoInverse(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, double singularMinToMaxRatio = .001, double dampingFactor = 0);
    Eigen::MatrixXd pseudoInverseMatrix(const Eigen::MatrixXd& A, double singularMinToMaxRatio = .001, double dampingFactor = 0);
}


#endif /* PSEUDOINVERSE_HPP_ */
