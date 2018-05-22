#ifndef ARE_SOLVER
#define ARE_SOLVER

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <assert.h>
#include <complex>


using Eigen::MatrixXd;


struct ARE_Solution
{
    MatrixXd X;
    MatrixXd K;
};

ARE_Solution dare(const MatrixXd& Ad, const MatrixXd& Bd, const MatrixXd& Qd, const MatrixXd& Rd);
ARE_Solution care(const MatrixXd& A, const MatrixXd& B, const MatrixXd& Q, const MatrixXd& R);

void discretizatize(double dt, const MatrixXd& A, const MatrixXd& B, const MatrixXd& Q, const MatrixXd& R, MatrixXd& Ad, MatrixXd& Bd, MatrixXd& Qd, MatrixXd& Rd);
MatrixXd expm(const MatrixXd& A);

#endif // ARE_SOLVER

