#include "are_solver.h"

ARE_Solution dare(const MatrixXd& A, const MatrixXd& B, const MatrixXd& Q, const MatrixXd& R)
{
    //https://en.wikipedia.org/wiki/Algebraic_Riccati_equation

    assert(A.rows() == A.cols() && A.rows() == B.rows() && A.rows() == Q.rows() && A.rows() == Q.cols() && B.cols() == R.rows() && B.cols() == R.cols() && "DARE Solver: Check input matrix sizes");

    int N = A.rows();
    int M = B.cols();

    MatrixXd AinvT = A.colPivHouseholderQr().solve(MatrixXd::Identity(N,N)).transpose();
    MatrixXd Rinv = R.ldlt().solve(MatrixXd::Identity(M,M));

    MatrixXd Z; // synaptic pencile

    Z.resize(2*N,2*N);

    // Z = [A+B*Rinv*B'*Ainv'*Q -B*Rinv*B'*Ainv'; -Ainv'*Q Ainv'];

    Z.block(0,N,N,N) = -B*Rinv*B.transpose()*AinvT;
    Z.block(0,0,N,N) = A-Z.block(0,N,N,N)*Q;
    Z.block(N,0,N,N) = -AinvT*Q;
    Z.block(N,N,N,N) = AinvT;


    // [V,L] = eig(Z)
    Eigen::EigenSolver<MatrixXd> eigZ(Z);
    //eigZ.compute(Z,true);

    Eigen::MatrixXcd U1;
    Eigen::MatrixXcd U2;
    U1.resize(N,N);
    U2.resize(N,N);

    Eigen::MatrixXcd eigenValues = eigZ.eigenvalues();


    //    U = zeros(size(Z,1),size(Z,2)/2);
    //    uind = 1;
    //    for i=1:size(L)
    //        if abs(L(i,i)) < 1
    //            U(:,uind) = V(:,i);
    //            uind = uind +1;
    //        end
    //    end
    int u_col = 0;
    for(int eigInd=0; eigInd < 2*N && u_col<N; eigInd++ )
    {
        if( std::abs(eigenValues(eigInd).real()) < 1 )
        {
            U1.block(0,u_col,N,1) = eigZ.eigenvectors().block(0,eigInd,N,1);
            U2.block(0,u_col,N,1) = eigZ.eigenvectors().block(N,eigInd,N,1);
            u_col ++;
        }
    }

    assert(u_col == N && "DARE ERROR: No Solution Found");



    //    X = U(floor(end/2+1):end,:)*U(1:end/2,:)^-1;
    //    K = (R+B'*X*B)^-1*(B'*X*A);
    ARE_Solution solution;
    //Eigen::MatrixXcd U1 = U.block(0,0,N,N);
    Eigen::MatrixXcd X_c = U2*U1.colPivHouseholderQr().solve(Eigen::MatrixXcd::Identity(N,N));

    solution.X.resize(N,N);
    for(int i=0; i<N;i++)
        for(int j=0;j<N;j++)
            solution.X(i,j) = X_c(i,j).real();

    solution.K = (R+B.transpose()*solution.X*B).colPivHouseholderQr().solve(B.transpose()*solution.X*A);

    return solution;

}

ARE_Solution care(const MatrixXd& A, const MatrixXd& B, const MatrixXd& Q, const MatrixXd& R)
{

    assert(A.rows() == A.cols() && A.rows() == B.rows() && A.rows() == Q.rows() && A.rows() == Q.cols() && B.cols() == R.rows() && B.cols() == R.cols() && "DARE Solver: Check input matrix sizes");

    int N = A.rows();
    int M = B.cols();

    // MatrixXd AinvT = A.colPivHouseholderQr().solve(MatrixXd::Identity(N,N)).transpose();
    MatrixXd Rinv = R.ldlt().solve(MatrixXd::Identity(M,M));

    MatrixXd Z; // synaptic pencile

    Z.resize(2*N,2*N);

    // Z = [A+B*Rinv*B'*Ainv'*Q -B*Rinv*B'*Ainv'; -Ainv'*Q Ainv']; https://en.wikipedia.org/wiki/Algebraic_Riccati_equation
    Z.block(0,0,N,N) = A;
    Z.block(0,N,N,N) = -B*Rinv*B.transpose();
    Z.block(N,0,N,N) = -Q;
    Z.block(N,N,N,N) = -A.transpose();


    // [V,L] = eig(Z)
    Eigen::EigenSolver<MatrixXd> eigZ(Z);
    //eigZ.compute(Z,true);

    Eigen::MatrixXcd U1;
    Eigen::MatrixXcd U2;
    U1.resize(N,N);
    U2.resize(N,N);

    Eigen::MatrixXcd eigenValues = eigZ.eigenvalues();


    //    U = zeros(size(Z,1),size(Z,2)/2);
    //    uind = 1;
    //    for i=1:size(L)
    //        if abs(L(i,i)) < 1
    //            U(:,uind) = V(:,i);
    //            uind = uind +1;
    //        end
    //    end
    int u_col = 0;
    for(int eigInd=0; eigInd < 2*N && u_col<N; eigInd++ )
    {
        if( eigenValues(eigInd).real() < 0 )
        {
            U1.block(0,u_col,N,1) = eigZ.eigenvectors().block(0,eigInd,N,1);
            U2.block(0,u_col,N,1) = eigZ.eigenvectors().block(N,eigInd,N,1);
            u_col ++;
        }
    }
    assert(u_col == N && "CARE ERROR: No Solution Found");


    //    X = U(floor(end/2+1):end,:)*U(1:end/2,:)^-1;
    //    K = (R+B'*X*B)^-1*(B'*X*A);
    ARE_Solution solution;
    //Eigen::MatrixXcd U1 = U.block(0,0,N,N);
    Eigen::MatrixXcd X_c = U2*U1.colPivHouseholderQr().solve(Eigen::MatrixXcd::Identity(N,N));

    solution.X.resize(N,N);
    for(int i=0; i<N;i++)
        for(int j=0;j<N;j++)
            solution.X(i,j) = X_c(i,j).real();

    solution.K = Rinv*B.transpose()*solution.X;

    return solution;

}

//#include <iostream>
//using namespace std;
MatrixXd expm(const MatrixXd& A)
{
    assert(A.rows()==A.cols() && "expm error: A must be square");
    int N = A.rows();
    MatrixXd soln;
    soln.setIdentity(A.rows(),A.cols());

    MatrixXd Amult = A;
    double nfact = 1;

    double s_minmax = 1;

    for( int i = 2; i < 100; i++ )
    {
        // check to see if we should terminate the series
        double a_minmax = 0;
        for( int j = 0; j < N ; j ++ )
            for( int k=0; k < N; k ++ )
                a_minmax = std::max(a_minmax,std::abs(Amult(j,k)));

        a_minmax /= nfact;
        if( a_minmax/s_minmax > 1e-10 )
        {
            s_minmax += a_minmax;
        } else
        {
            break;
        }

        soln += Amult/nfact;
        Amult = Amult*A;
        nfact *= i;

        // cout << Amult/nfact <<endl<<endl;
    }

    return soln;
}

//#include <iostream>
//using namespace std;
void discretizatize(double dt, const MatrixXd& A, const MatrixXd& B, const MatrixXd& Q, const MatrixXd& R, MatrixXd& Ad, MatrixXd& Bd, MatrixXd& Qd, MatrixXd& Rd)
{
    //https://en.wikipedia.org/wiki/Discretization

    assert(A.rows() == A.cols() && A.rows() == B.rows() && A.rows() == Q.rows() && A.rows() == Q.cols() && B.cols() == R.rows() && B.cols() == R.cols() && "Discretization: Check input matrix sizes");

    int N = A.rows();
    int M = B.cols();

    Eigen::MatrixXd AD;
    AD.setZero(N+M,N+M);
    AD.block(0,0,N,N) = A;
    AD.block(0,N,N,M) = B;

    Eigen::MatrixXd tmp;
    tmp = expm(AD*dt);

    Ad = tmp.block(0,0,N,N);
    Bd = tmp.block(0,N,N,M);

    // convert Q
    AD.setZero(2*N,2*N);
    AD.block(0,0,N,N) = -A;
    AD.block(0,N,N,N) = Q;
    AD.block(N,N,N,N) = A.transpose();

    tmp = expm(AD*dt);
    Qd = tmp.block(N,N,N,N).transpose()*tmp.block(0,N,N,N);
   // cout << Ad << endl << tmp.block(N,N,N,N).transpose() <<endl;

    // convert R
    Rd = R/dt;

    return;


}
