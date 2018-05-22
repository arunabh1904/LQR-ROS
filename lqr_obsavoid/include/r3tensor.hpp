#ifndef R3TENSOR_HPP
#define R3TENSOR_HPP
#include "r3tensor.h"

#include "iostream"
using namespace std;

using namespace Eigen;


template<int R, int C, int D>
R3Tensor<R,C,D>::R3Tensor()
{
    return;
}

template<int R, int C, int D>
R3Tensor<R,C,D>::R3Tensor(const R3Tensor<R,C,D>& other)
{
    tensor_storage=other.tensor_storage;
}

//template<int R, int C, int D>
//R3Tensor<R,C,D>::R3Tensor(const Eigen::MatrixXd& mat)
//{
//    assert(mat.rows() == R*D && mat.cols()==C && "Error in R3Tensor MatrixXd constructor, dims dont match");

//    tensor_storage = mat;//.block(0,0,R*D,C);
//}

template<int R, int C, int D>
R3Tensor<R,C,D>::R3Tensor(const Eigen::Matrix<double,R*D,C>& mat)
{
    tensor_storage = mat;
}

template<int R, int C, int D>
Eigen::Matrix<double,R*D,C>& R3Tensor<R,C,D>::operator() ()
{
    return tensor_storage;
}

template<int R, int C, int D>
Eigen::Matrix<double,R*D,C> R3Tensor<R,C,D>::operator() () const
{
    return tensor_storage;
}

template<int R, int C, int D>
Eigen::MatrixXd R3Tensor<R,C,D>::toMatrixXd() const
{
    return tensor_storage;
}

template<int R, int C, int D>
R3Tensor<R,C,D> R3Tensor<R,C,D>::Zero()
{
    R3Tensor<R,C,D> retTen;
    retTen.setZero();
    return retTen;
}


template<int R, int C, int D>
void R3Tensor<R,C,D>::setZero()
{
    tensor_storage.setZero(R*D,C);

    return;
}


template<int R, int C, int D>
Eigen::Block<Eigen::Matrix<double,R*D,C>,R,C> R3Tensor<R,C,D>::rc_slice(int depth)
{
    assert( depth >=0 && depth < D && "R3Tensor::rc_slice depth must be in range");

    return Eigen::Block< Eigen::Matrix<double,R*D,C>, R, C>(tensor_storage.derived(), depth*R, 0);
}

template<int R, int C, int D>
Eigen::Block<const Eigen::Matrix<double,R*D,C>,R,C> R3Tensor<R,C,D>::rc_slice(int depth) const
{
    assert( depth >=0 && depth < D && "R3Tensor::rc_slice depth must be in range");

    return Eigen::Block< const Eigen::Matrix<double,R*D,C>, R, C>(tensor_storage.derived(), depth*R, 0);
}

template<int R, int C, int D>
Eigen::Map< Eigen::Matrix<double,R,D> ,0, Stride<0,0> > R3Tensor<R,C,D>::rd_slice(int col)
{
    assert( col >=0 && col < C && "R3Tensor::rd_slice column must be in range");
    Eigen::Map< Eigen::Matrix<double,R,D> ,0, Stride<0,0> > rd_Map(tensor_storage.data()+col*R*D,R,D);
    return rd_Map;
}

template<int R, int C, int D>
Eigen::Map< const Eigen::Matrix<double,R,D> ,0, Stride<0,0> > R3Tensor<R,C,D>::rd_slice(int col) const
{
    assert( col >=0 && col < C && "R3Tensor::rd_slice column must be in range");
    //Eigen::Matrix<double,R,D> rd_Map;
    Eigen::Map< const Eigen::Matrix<double,R,D> ,0, Stride<0,0> > rd_Map(tensor_storage.data()+col*R*D,R,D);

    return rd_Map;
}

template<int R, int C, int D>
Eigen::Map< Eigen::Matrix<double,D,C> ,0, Eigen::Stride<R*D,R> > R3Tensor<R,C,D>::dc_slice(int row)   ///<- Returns a copy of the RxD matrix slice that cannot be directly modified
{
    assert( row >=0 && row < R && "R3Tensor::cd_slice row must be in range");

    Eigen::Map< Eigen::Matrix<double,D,C> ,0, Eigen::Stride<R*D,R> > dc_Map(tensor_storage.data()+row,D,C);
    return dc_Map;

}

template<int R, int C, int D>
Eigen::Map<const Eigen::Matrix<double,D,C> ,0, Eigen::Stride<R*D,R> > R3Tensor<R,C,D>::dc_slice(int row) const   ///<- Returns a copy of the RxD matrix slice that cannot be directly modified
{
    assert( row >=0 && row < R && "R3Tensor::cd_slice row must be in range");

    Eigen::Map< const Eigen::Matrix<double,D,C> ,0, Eigen::Stride<R*D,R> > dc_Map(tensor_storage.data()+row,D,C);
    return dc_Map;

}



template<int R, int C, int D>
double& R3Tensor<R,C,D>::operator()(int r, int c, int d)
{
    assert( c >=0 && c < C && "R3Tensor::operator() column must be in range");
    assert( r >=0 && r < R && "R3Tensor::operator() row must be in range");
    assert( d >=0 && d < D && "R3Tensor::operator() depth must be in range");
    return tensor_storage(r+R*d,c);
}

template<int R, int C, int D>
double R3Tensor<R,C,D>::operator()(int r, int c, int d) const
{
    assert( c >=0 && c < C && "R3Tensor::operator() column must be in range");
    assert( r >=0 && r < R && "R3Tensor::operator() row must be in range");
    assert( d >=0 && d < D && "R3Tensor::operator() depth must be in range");
    return tensor_storage(r+R*d,c);
}

template<int R, int C, int D>
R3Tensor<R,C,D>& R3Tensor<R,C,D>::depthMultiplyInPlace(const Eigen::Matrix<double,D,D>& RHS)
{
    for(int i=0; i<C;i++)
        rd_slice(i) *= RHS;

    return *this;
}

template<int R, int C, int D>
template<int RHS_C>
R3Tensor<R,C,RHS_C> R3Tensor<R,C,D>::depthMultiply(const Eigen::Matrix<double,D,RHS_C>& RHS) const
{

    R3Tensor<R,C,RHS_C> copy;
    for(int i=0; i<C;i++)
        copy.rd_slice(i) = rd_slice(i)*RHS;

    return copy;

}


template<int R, int C, int D>
template<int LHS_C>
R3Tensor<R,LHS_C,D> R3Tensor<R,C,D>::operator*(const Eigen::Matrix<double,C,LHS_C>& RHS ) const
{
    R3Tensor<R,LHS_C,D> retTensor(tensor_storage*RHS);


    return retTensor;

}

template<int R, int C, int D, int LHS_R>
R3Tensor<LHS_R,C,D> operator*(const Eigen::Matrix<double,LHS_R,R>& LHS, const R3Tensor<R,C,D>& RHS )
{
    R3Tensor<LHS_R,C,D> soln;
    for(int i=0; i< D; i++)
        soln.rc_slice(i) = LHS*RHS.rc_slice(i);
    return soln;
}


template<int R, int C, int D>
R3Tensor<R,C,D>& R3Tensor<R,C,D>::operator*=(const Eigen::Matrix<double,C,R>& RHS )
{
    tensor_storage*=RHS;

    return *this;
}

template<int R, int C, int D>
R3Tensor<R,C,D> operator*(double val, R3Tensor<R,C,D> tensor)
{
    return tensor*val;
}

template<int R, int C, int D>
Eigen::Matrix<double,D,C> operator*(const Eigen::RowVectorXd& LHS, const R3Tensor<R,C,D>& RHS )
{
    assert(LHS.cols() == R && "Row Vector must have as many columns as tensor rows");

    Eigen::Matrix<double,D,C> retMat;

    for(int i=0; i<D; i++)
        retMat.block(i,0,1,C) = LHS*RHS.rc_slice(i);

    return retMat;
}

template<int R, int C, int D>
R3Tensor<R,C,D>& R3Tensor<R,C,D>::operator*=(double RHS)
{
    tensor_storage*=RHS;
    return *this;
}

template<int R, int C, int D>
R3Tensor<R,C,D>& R3Tensor<R,C,D>::operator/=(double RHS)
{
    tensor_storage/=RHS;
    return *this;
}

template<int R, int C, int D>
R3Tensor<R,C,D>& R3Tensor<R,C,D>::operator+=(const R3Tensor&RHS)
{

    tensor_storage+=RHS.tensor_storage;

    return *this;

}

template<int R, int C, int D>
R3Tensor<R,C,D>& R3Tensor<R,C,D>::operator-=(const R3Tensor&RHS)
{
    tensor_storage-=RHS.tensor_storage;

    return *this;

}

template<int R, int C, int D>
R3Tensor<R,C,D>& R3Tensor<R,C,D>::operator=(const R3Tensor<R,C,D>& other)
{
    tensor_storage = other.tensor_storage;

    return *this;
}


template<int R, int C, int D>
R3Tensor<R,C,D> R3Tensor<R,C,D>::operator*(double RHS) const
{
    R3Tensor<R,C,D> retTensor(*this);
    retTensor *=RHS;
    return retTensor;

}

template<int R, int C, int D>
R3Tensor<R,C,D> R3Tensor<R,C,D>::operator/(double RHS) const
{
    R3Tensor<R,C,D> retTensor(*this);
    retTensor /=RHS;
    return retTensor;

}

template<int R, int C, int D>
R3Tensor<R,C,D> R3Tensor<R,C,D>::operator+(const R3Tensor& RHS) const
{
    R3Tensor<R,C,D> retTensor(*this);
    retTensor +=RHS;
    return retTensor;

}

template<int R, int C, int D>
R3Tensor<R,C,D> R3Tensor<R,C,D>::operator-(const R3Tensor& RHS) const
{
    R3Tensor<R,C,D> retTensor(*this);
    retTensor -=RHS;
    return retTensor;

}

template<int R, int C, int D>
R3Tensor<R,D,C> R3Tensor<R,C,D>::cd_transpose() const
{
    R3Tensor<R,D,C> retTensor;

    for(int i = 0;i<C;i++)
        retTensor.rc_slice(i) = rd_slice(i);

    return retTensor;
}

template<int R, int C, int D>
R3Tensor<D,C,R> R3Tensor<R,C,D>::rd_transpose() const
{
    R3Tensor<D,C,R> retTensor;

    for(int i = 0;i<C;i++)
        retTensor.rd_slice(i) = rd_slice(i).transpose();

    return retTensor;
}


template<int R, int C, int D>
double R3Tensor<R,C,D>::norm() const
{
    return tensor_storage.norm();
}

#endif // R3TENSOR_HPP

