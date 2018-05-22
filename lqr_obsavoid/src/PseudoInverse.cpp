#include "PseudoInverse.h"
#include <iostream>
#include <cmath>

Eigen::VectorXd Math::pseudoInverse(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, double singularMinToMaxRatio, double dampingFactor)
{
    // uses a singular value decomposisition (SVD) to solve for the least squares solution

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV ); // computes the SVD


    //    cout << "SVD COMPUTED" << endl;
    Eigen::MatrixXd S_inv = Eigen::MatrixXd::Zero(svd.matrixV().cols(),svd.matrixU().cols());

   // std::cout << svd.singularValues() << std::endl;
    for( int i = 0; i < std::min(A.rows(),A.cols()); i++)
    {
        double val = 0;
        if( svd.singularValues()[i] > svd.singularValues()[0]*singularMinToMaxRatio )// threashold singular values anything less than 1/1000 of the max is set to 0
            //val = 1.0 / svd.singularValues()[i];
            val = svd.singularValues()[i] / (svd.singularValues()[i]*svd.singularValues()[i]+dampingFactor);

        S_inv(i,i) = val;
    }
    //    cout << "Singular Values Flipped" << endl;

    Eigen::VectorXd answer;

    //        cout << "Size of U:  " << svd.matrixU().rows() << " x " << svd.matrixU().cols() << endl;
    //        cout << "Size of V:  " << svd.matrixV().rows() << " x " << svd.matrixV().cols() << endl;
    //        cout << "Size of S:  " << S.rows() << " x " << S.cols() << endl;
    //        cout << "Size of b:  " << b.rows() << " x " << b.cols() << endl;


    answer = svd.matrixV()*(S_inv*(svd.matrixU().transpose()*b));

    //   cout << "Answer Computed " << endl;
    return answer;
}


Eigen::MatrixXd Math::pseudoInverseMatrix(const Eigen::MatrixXd& A, double singularMinToMaxRatio, double dampingFactor)
{
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV ); // computes the SVD


    //    cout << "SVD COMPUTED" << endl;
    Eigen::MatrixXd S_inv = Eigen::MatrixXd::Zero(svd.matrixV().cols(),svd.matrixU().cols());
    /*
    std:: cout << std::endl;
    std:: cout << A << std::endl << std::endl;
    std::cout << svd.matrixV() << std::endl << std::endl;
    std::cout << svd.matrixU() << std::endl << std::endl;
    std::cout << svd.singularValues() << std::endl;
    */

    double minAcceptableVal = svd.singularValues()[0]*singularMinToMaxRatio;
    for( int i = 0; i < std::min(A.rows(),A.cols()); i++)
    {
        double val = 0;
        if( svd.singularValues()[i] > minAcceptableVal )
        {   // threashold singular values anything less than singularMinToMaxRatio of the max is set to 0
            val = svd.singularValues()[i] / (svd.singularValues()[i]*svd.singularValues()[i]+dampingFactor);
        }

        S_inv(i,i) = val;
    }
    //    cout << "Singular Values Flipped" << endl;

    Eigen::MatrixXd answer;

    //    cout << "Size of A:  " << A.rows() << " x " << A.cols() << endl;
    //    cout << "Size of U:  " << svd.matrixU().rows() << " x " << svd.matrixU().cols() << endl;
    //    cout << "Size of V:  " << svd.matrixV().rows() << " x " << svd.matrixV().cols() << endl;
    //    cout << "Size of Sinv:  " << S_inv.rows() << " x " << S_inv.cols() << endl;
    //    cout << endl << S_inv << endl << endl;


    answer = svd.matrixV()*(S_inv*(svd.matrixU().transpose()));

    //std::cout << "Answer Computed " << std::endl << answer << "\n__________________\n"<<  A*answer << "\n__________________\n" << answer*A << std::endl<<std::endl;
    return answer;
}
