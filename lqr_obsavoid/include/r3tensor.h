#ifndef R3TENSOR_H
#define R3TENSOR_H

#include <Eigen/Dense>
#include <assert.h>
#include <ostream>

/**
 *  R3Tensor implements the functionality necessary to implement the Rank 3 Tensor Jacobians for MPC control
 */
template<int R, int C, int D>
class R3Tensor
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW // need to do this because we dynamically allociate this object and it has a fixed sized Eigen Matrix

protected:
    Eigen::Matrix<double,R*D,C> tensor_storage;

public:
    R3Tensor();
    R3Tensor(const R3Tensor& other);
    R3Tensor(const Eigen::Matrix<double,R*D,C>& mat);
    //R3Tensor(const Eigen::MatrixXd& mat);

    static R3Tensor Zero();

    void setZero(); ///<- Sets all cells to zero

    Eigen::Matrix<double,R*D,C> operator() ()const;
    Eigen::Matrix<double,R*D,C>& operator() ();

    Eigen::MatrixXd toMatrixXd() const; /// Returns a (R*D)xC matrix where each depth slice is appended as additional rows


    Eigen::Block<Eigen::Matrix<double,R*D,C>,R,C> rc_slice(int depth);       ///<- Returns a reference to the RxC matrix slice that can be modified
    Eigen::Block<const Eigen::Matrix<double,R*D,C>,R,C>  rc_slice(int depth) const; ///<- Returns a copy of the RxC matrix slice that cannot be directly modified

    Eigen::Map< Eigen::Matrix<double,R,D> ,0, Eigen::Stride<0,0> > rd_slice(int col) ;   ///<- Returns a copy of the RxD matrix slice that can be directly modified
    Eigen::Map< const Eigen::Matrix<double,R,D> ,0, Eigen::Stride<0,0> > rd_slice(int col) const;   ///<- Returns a copy of the RxD matrix slice that cannot be directly modified

    Eigen::Map< Eigen::Matrix<double,D,C> ,0, Eigen::Stride<R*D,R> > dc_slice(int col) ;   ///<- Returns a copy of the DxC matrix slice that can be directly modified
    Eigen::Map< const Eigen::Matrix<double,D,C> ,0, Eigen::Stride<R*D,R> > dc_slice(int col) const;    ///<- Returns a copy of the DxC matrix slice that cannot be directly modified

    R3Tensor<R,D,C> cd_transpose() const; ///<- Returns a transpose where the depth and column indices are swapped.
    R3Tensor<D,C,R> rd_transpose() const; ///<- Returns a transpose where the depth and row indices are swapped.

    double& operator()(int r, int c, int d); ///<- provides element access and allows modification
    double operator()(int r, int c, int d) const; ///<- provides element access without modification

    R3Tensor<R,C,D>& depthMultiplyInPlace(const Eigen::Matrix<double,D,D>& RHS); ///<- returns the current matrix multiplied by an matrix in the depth direction

    template<int RHS_C>
    R3Tensor<R,C,RHS_C> depthMultiply(const Eigen::Matrix<double,D,RHS_C>& RHS) const; ///<- returns the a matrix multiplied by an matrix in the depth direction


    template<int LHS_C>
    R3Tensor<R,LHS_C,D> operator*(const Eigen::Matrix<double,C,LHS_C>& RHS ) const; ///<- returns the tensor multiplied on the right by a matrix

    R3Tensor& operator*=(const Eigen::Matrix<double,C,R>& RHS ); ///<- returns the tensor multiplied on the right by a matrix

    R3Tensor& operator*=(double);
    R3Tensor& operator/=(double);
    R3Tensor& operator+=(const R3Tensor&);
    R3Tensor& operator-=(const R3Tensor&);
    R3Tensor& operator=(const R3Tensor& );


    R3Tensor operator*(double) const;
    R3Tensor operator/(double) const;
    R3Tensor operator+(const R3Tensor& ) const;
    R3Tensor operator-(const R3Tensor& ) const;

    double norm() const; /// The Frobenius norm

};

template<int R, int C, int D>
R3Tensor<R,C,D> operator*(double val, R3Tensor<R,C,D> tensor);

template<int R, int C, int D>
Eigen::Matrix<double,D,C> operator*(const Eigen::RowVectorXd& LHS, const R3Tensor<R,C,D>& RHS );

template<int R, int C, int D, int LHS_R>
R3Tensor<LHS_R,C,D> operator*(const Eigen::Matrix<double,LHS_R,R>& LHS, const R3Tensor<R,C,D>& RHS );

template<int R, int C, int D>
std::ostream& operator<<(std::ostream& os, const R3Tensor<R,C,D> & tensor)
{
    for( int i=0; i<D; i++)
        os << tensor.rc_slice(i) << std::endl << std::endl;

    return os;
}

#include "r3tensor.hpp"

#endif // R3TENSOR_H
