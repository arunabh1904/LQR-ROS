#include "NewtonRaphson.h"

template<>
NewtonRaphson<_nr_default>::NewtonRaphson(jacobian_update_method update_method, int updateJacobianEvery )
{
    iteration = 0;

    H.resize(0,0); /**< Teh scaling Matrix **/
    J.resize(0,0); /**< The current jacobian **/
    J_inv.resize(0,0); /**< The current inverse jacobian **/
    states.resize(0); /**< The current state vector **/
    delta_States.resize(0); /**< The last change in states **/
    error_last.resize(0); /**< The previous step's error vector **/
    error_this.resize(0); /**< The current error vector **/
    delta_error.resize(0); /**< The current change in error **/

    this->errorFcnPtr = 0;
    this->errorJacobianFcnPtr = 0;
    this->errorFcnPtr_2 = 0;
    this->errorJacobianFcnPtr_2 = 0;
    this->update_method = update_method;

    this->classObjectPointer = 0;
    this->class_errorFcnPtr = 0;
    this->class_errorJacobianFcnPtr = 0;

    this->modelReductionThreshold = -1;

    this->updateJacobianEvery = updateJacobianEvery;
    sqrtEps = std::sqrt(std::numeric_limits<double>::epsilon());

}

template<>
NewtonRaphson<_nr_default>::NewtonRaphson(  Eigen::VectorXd (*errorFcnPtr)(const Eigen::VectorXd&), Eigen::MatrixXd (*errorJacobianFcnPtr)( const Eigen::VectorXd& ), jacobian_update_method update_method, int updateJacobianEvery, double modelReductionThreshold_ )
{
    iteration = 0;

    H.resize(0,0); /**< Teh scaling Matrix **/
    J.resize(0,0); /**< The current jacobian **/
    J_inv.resize(0,0); /**< The current inverse jacobian **/
    states.resize(0); /**< The current state vector **/
    delta_States.resize(0); /**< The last change in states **/
    error_last.resize(0); /**< The previous step's error vector **/
    error_this.resize(0); /**< The current error vector **/
    delta_error.resize(0); /**< The current change in error **/

    this->errorFcnPtr = errorFcnPtr;
    this->errorJacobianFcnPtr = errorJacobianFcnPtr;
    this->errorFcnPtr_2 = 0;
    this->errorJacobianFcnPtr_2 = 0;
    this->update_method = update_method;

    this->classObjectPointer = 0;
    this->class_errorFcnPtr = 0;
    this->class_errorJacobianFcnPtr = 0;


    this->updateJacobianEvery = updateJacobianEvery;

    this->modelReductionThreshold = modelReductionThreshold_;

    sqrtEps = std::sqrt(std::numeric_limits<double>::epsilon());


    assert(errorFcnPtr != 0);
}

template<>
NewtonRaphson<_nr_default>::NewtonRaphson(  void (*errorFcnPtr)(Eigen::VectorXd&, Eigen::VectorXd&), void (*errorJacobianFcnPtr)( Eigen::VectorXd&, Eigen::MatrixXd& ), jacobian_update_method update_method, int updateJacobianEvery, double modelReductionThreshold_ )
{
    iteration = 0;
    H.resize(0,0); /**< The scaling Matrix **/
    J.resize(0,0); /**< The current jacobian **/
    J_inv.resize(0,0); /**< The current inverse jacobian **/
    states.resize(0); /**< The current state vector **/
    delta_States.resize(0); /**< The last change in states **/
    error_last.resize(0); /**< The previous step's error vector **/
    error_this.resize(0); /**< The current error vector **/
    delta_error.resize(0); /**< The current change in error **/

    this->errorFcnPtr = 0;
    this->errorJacobianFcnPtr = 0;
    this->errorFcnPtr_2 = errorFcnPtr;
    this->errorJacobianFcnPtr_2 = errorJacobianFcnPtr;

    this->classObjectPointer = 0;
    this->class_errorFcnPtr = 0;
    this->class_errorJacobianFcnPtr = 0;


    this->update_method = update_method;

    this->updateJacobianEvery = updateJacobianEvery;

    this->modelReductionThreshold = modelReductionThreshold_;
    sqrtEps = std::sqrt(std::numeric_limits<double>::epsilon());

    assert(errorJacobianFcnPtr_2 != 0 );
    assert(errorFcnPtr_2 != 0);
}

template<>
NewtonRaphson<_nr_default>::NewtonRaphson( _nr_default* classPointer, void (_nr_default::*errorFcnPtr)(Eigen::VectorXd&, Eigen::VectorXd&), void (_nr_default::*errorJacobianFcnPtr)( Eigen::VectorXd&, Eigen::MatrixXd& ), jacobian_update_method update_method, int updateJacobianEvery, double modelReductionThreshold_ )
{
    cout << "You are using the wrong NewtonRapson constructor.  Use the non-class specific one!" << endl;
    assert(false);
}
//void NewtonRaphson::updateJacobian(int count )
//{
//    // update inverse jacobian
//    jacobian_update_method method = update_method;

//    if( count == 0 || update_method == EVALUATE || updateJacobianEvery == 0 || count%updateJacobianEvery == 0 || J_inv.rows()==0 || J_inv.cols()==0)
//    {
//        method = EVALUATE; // recalculate the jacobian
//    }


//    switch( method )
//    {
//    case EVALUATE:
//    {
//        if( errorJacobianFcnPtr )
//        {
//            J = (*errorJacobianFcnPtr)(states);
//        }
//        else if( errorJacobianFcnPtr_2 )
//        {
//            (*errorJacobianFcnPtr_2)(states, J);
//        } else
//        {
//            J = CalculateJacobian(errorFcnPtr,states);
//        }

//        // Remove Zero rows
//        HtoJ_index.clear();
//        H.resize(J.rows(),J.cols());
//        unsigned int count = 0;
//        for( unsigned int i=0; i<J.cols(); i++ )
//            if( J.col(i).norm() != 0 )
//            {
//                H.col(count) = J.col(i);
//                HtoJ_index.push_back(i);
//                count ++;
//            }

//        J_inv.setZero(J.cols(),J.rows());
//        if( modelReductionThreshold <= 0 )
//        {
//            if( H.rows() >= count )
//                J_QR.compute(H.block(0,0,H.rows(),count).transpose()*H.block(0,0,H.rows(),count));
//            else
//                J_QR.compute(H.block(0,0,H.rows(),count)*H.block(0,0,H.rows(),count).transpose());
//        }

//        //        if( update_method != EVALUATE )
//        //        {
//        //            Eigen::MatrixXd ds_tmp;
//        //            if( modelReductionThreshold <= 0 )
//        //            {
//        //                if( H.rows() >= count )
//        //                    ds_tmp = J_QR.compute(H.block(0,0,H.rows(),count).transpose()*H.block(0,0,H.rows(),count)).solve(H.block(0,0,H.rows(),count).transpose());
//        //                else
//        //                    ds_tmp = H.block(0,0,H.rows(),count).transpose()*J_QR.compute(H.block(0,0,H.rows(),count)*H.block(0,0,H.rows(),count).transpose()).solve(Eigen::MatrixXd::Identity(H.rows(),H.rows()));
//        //            } else
//        //            {
//        //                //cout << H << endl;
//        //                ds_tmp = Math::pseudoInverseMatrix(H.block(0,0,H.rows(),count), modelReductionThreshold);
//        //            }

//        //            J_inv.setZero(J.cols(),J.rows());
//        //            for(unsigned int i=0; i<HtoJ_index.size(); i++ )
//        //            {
//        //                J_inv.row(HtoJ_index[i]) = ds_tmp.row(i);
//        //            }
//        //        }else
//        //        {
//        //            if( modelReductionThreshold <= 0 )
//        //            {
//        //                if( H.rows() >= HtoJ_index.size() )
//        //                    J_QR.compute(H.block(0,0,H.rows(),count).transpose()*H.block(0,0,H.rows(),count));
//        //                else
//        //                    J_QR.compute(H.block(0,0,H.rows(),count)*H.block(0,0,H.rows(),count).transpose());
//        //            }
//        //        }

//        break;
//    }
//    case APPROXIMATE_BROYDEN_1:
//        //https://en.wikipedia.org/wiki/Broyden%27s_method
//        if( count != 0 )
//        {
//            Eigen::VectorXd j_inv_delta_f = J_inv*(delta_error);
//            J_inv += (delta_States - j_inv_delta_f)/(delta_States.dot(j_inv_delta_f))*(delta_States.transpose()*J_inv);
//        }
//        break;

//    case APPROXIMATE_BROYDEN_2:
//        //https://en.wikipedia.org/wiki/Broyden%27s_method
//        if( count != 0 )
//        {
//            Eigen::VectorXd j_inv_delta_f = J_inv*delta_error;
//            J_inv += (delta_States - j_inv_delta_f)/(pow(delta_error.norm(),2.0))*(delta_error.transpose());
//        }
//        break;
//    }
//}

//void NewtonRaphson::calculateDeltaStates( )
//{
//    // Timing::stopWatch tic;

//    Eigen::VectorXd ds_tmp;
//    if( modelReductionThreshold <= 0 )
//    {
//        if( H.rows() >= (int) HtoJ_index.size() )
//            ds_tmp = -J_QR.solve(H.block(0,0,H.rows(),HtoJ_index.size()).transpose()*error_this);
//        else
//            ds_tmp = -H.block(0,0,H.rows(),HtoJ_index.size()).transpose()*J_QR.solve(error_this);

//    } else
//    {
//        ds_tmp = -Math::pseudoInverse(H.block(0,0,H.rows(),HtoJ_index.size()), error_this, modelReductionThreshold);
//    }

//    delta_States.setZero(J.cols());
//    for( unsigned int i=0; i<HtoJ_index.size(); i++ )
//    {
//        delta_States(HtoJ_index[i]) = ds_tmp(i);
//    }

//    switch( update_method )
//    {
//    case EVALUATE:
//    {
//        //        if( modelReductionThreshold <= 0 )
//        //        {
//        //            if( H.rows() >= HtoJ_index.size() )
//        //                ds_tmp = -J_QR.solve(H.block(0,0,H.rows(),HtoJ_index.size()).transpose()*error_this);
//        //            else
//        //                ds_tmp = -H.block(0,0,H.rows(),HtoJ_index.size()).transpose()*J_QR.solve(error_this);

//        //        } else
//        //        {
//        //            ds_tmp = -Math::pseudoInverse(H.block(0,0,H.rows(),HtoJ_index.size()), error_this, modelReductionThreshold);
//        //        }

//        //        delta_States.setZero();
//        //        for( unsigned int i=0; i<HtoJ_index.size(); i++ )
//        //        {
//        //            delta_States(HtoJ_index[i]) = ds_tmp(i);
//        //        }
//        break;
//    }

//    case APPROXIMATE_BROYDEN_1:
//        delta_States += -1.0*J_inv * error_this;
//        break;

//    case APPROXIMATE_BROYDEN_2:
//        delta_States += -1.0*J_inv * error_this;
//        break;
//    }

//    //std::cout << "Time for Inverse: " << tic.stop() << std::endl;
//}

//void NewtonRaphson::updateError()
//{
//    // get current error
//    error_last = error_this;

//    if( errorFcnPtr )
//        error_this = (*errorFcnPtr)(states);
//    else
//        (*errorFcnPtr_2)(states, error_this);

//}

//Eigen::VectorXd NewtonRaphson::solve(const Eigen::VectorXd& initial_guess, double tolerance, unsigned int maxAttempts, bool displayProgress)
//{
//    if( initial_guess.rows() == states.rows() )
//        iteration &= iteration && (initial_guess-states).norm()/states.norm() < .01;
//    else
//        iteration = 0;

//    if(iteration == 0 )
//    {
//        states = initial_guess;
//        updateError();

//        states_last = states;
//        delta_States.setZero(initial_guess.rows(),initial_guess.cols());
//        delta_error.setZero(error_this.rows());
//        error_last = error_this;

//        delta_States_last.setZero(delta_States.rows(),1);

//    } else
//    {
//        states = initial_guess;
//        updateError();
//        delta_States = states-states_last;
//        delta_error = error_this-error_last;

//        delta_States_last.setZero(delta_States.rows(),1);

//    }
//    bool oscilitory = false;
//    bool initialPass = true;
//    double percentErrorChange  =0;
//    do {
//        if( displayProgress )
//            std::cout << std::setprecision(3) << "Iteration: " << iteration << "\tCurrent Error: " << solutionRMSError() << "\tDelta Error: " <<  delta_error.norm()/delta_States.rows()  << "\tPercent Error Change: " << (delta_error.norm()/delta_States.rows()) / solutionRMSError() * 100<< "%" << std::endl;


//        updateJacobian(iteration);
//        removeNullSpace();
//        calculateDeltaStates();

//        if( !initialPass && (delta_States+delta_States_last).norm() < tolerance * 10.0  )
//        {
//            //cout << "d states: " << (delta_States+delta_States_last).norm() << "\t\t" << (delta_States+delta_States_last).transpose() << endl;

//            delta_States = (delta_States+delta_States_last)/2.0;
//            oscilitory = true;
//        }

//        states += delta_States;

//        updateError();

//        delta_States_last = delta_States;

//        delta_States = states-states_last; // just incase BC/Constraint inforcement changed something
//        states_last = states;

//        delta_error = error_this-error_last;

//        if( !initialPass && (error_this+error_last).norm() < tolerance*10.0 )
//        {
//            // cout << "e Error: " << (error_this+error_last).norm() << "\t\t" << (error_this+error_last).transpose() << endl;

//            //delta_States = (delta_States+delta_States_last)/2.0;
//            oscilitory = true;
//        }


//        iteration ++;

//        initialPass = false;


//        percentErrorChange = (delta_error.norm()/delta_States.rows()) / solutionRMSError() * 100;
//        if( isnan(percentErrorChange) )
//            percentErrorChange = 0;

//    }while( !oscilitory && solutionRMSError() > tolerance && (sqrt(delta_error.squaredNorm()/delta_error.rows())*10.0 > tolerance || delta_error.norm()==0) && iteration < maxAttempts && (percentErrorChange > 5e-2 || percentErrorChange == 0 ) );

//    /*if( H.rows() < HtoJ_index.size() && update_method == EVALUATE)
//    {
//        updateJacobian(iteration);
//        removeNullSpace();
//        updateError();
//    }*/


//    if( displayProgress )
//        std::cout << std::setprecision(3) << "Iteration: " << iteration << "\tCurrent Error: " << solutionRMSError() << "\tPercent Error Change: " <<  percentErrorChange <<"%"  <<"\tTolerance: " << tolerance << "\tDONE!"<<std::endl;


//    return states;
//}

//Eigen::VectorXd NewtonRaphson::solution() const
//{
//    return states;
//}

//double NewtonRaphson::solutionRMSError() const
//{
//    return sqrt(error_this.squaredNorm()/error_this.rows());
//}

//void NewtonRaphson::changeUpdateMethod(jacobian_update_method newMethod, int updateJacobianEvery )
//{
//    this->updateJacobianEvery = updateJacobianEvery;
//    this->update_method = newMethod;
//}

//void NewtonRaphson::removeNullSpace()
//{
//    if( H.rows() < HtoJ_index.size() && update_method == EVALUATE)
//    {
//        Eigen::VectorXd ds_tmp(HtoJ_index.size());
//        for( unsigned int i=0; i<HtoJ_index.size(); i++ )
//        {
//            ds_tmp(i) = states(HtoJ_index[i]);
//        }
//        if( modelReductionThreshold <= 0 )
//        {
//            ds_tmp = H.block(0,0,H.rows(),HtoJ_index.size()).transpose()*J_QR.solve(H.block(0,0,H.rows(),HtoJ_index.size())*ds_tmp);
//        } else
//        {
//            //cout << H.block(0,0,H.rows(),HtoJ_index.size())*H.block(0,0,H.rows(),HtoJ_index.size()).transpose() << endl << "____" <<endl;
//            //cout << H.block(0,0,H.rows(),HtoJ_index.size())*ds_tmp << endl<< "----" << endl;
//            J_QR.compute(H.block(0,0,H.rows(),HtoJ_index.size())*H.block(0,0,H.rows(),HtoJ_index.size()).transpose());
//            ds_tmp = H.block(0,0,H.rows(),HtoJ_index.size()).transpose()*J_QR.solve(H.block(0,0,H.rows(),HtoJ_index.size())*ds_tmp);
//            //ds_tmp = H.block(0,0,H.rows(),HtoJ_index.size()).transpose()*Math::pseudoInverse(H.block(0,0,H.rows(),HtoJ_index.size())*H.block(0,0,H.rows(),HtoJ_index.size()).transpose(), H.block(0,0,H.rows(),HtoJ_index.size())*ds_tmp, modelReductionThreshold);
//        }


//        for( unsigned int i=0; i<HtoJ_index.size(); i++ )
//        {
//            states(HtoJ_index[i]) = ds_tmp(i);
//        }
//    }

//}

//void NewtonRaphson::setModelReductionThreshold(double thresh )
//{
//    modelReductionThreshold = thresh;
//}

//void NewtonRaphson::changeJacobianFunction( Eigen::MatrixXd (*errorJacobianFcnPtr)( const Eigen::VectorXd& ) )
//{
//    this->errorJacobianFcnPtr = errorJacobianFcnPtr;
//}

//void NewtonRaphson::changeErrorFunction( Eigen::VectorXd (*errorFcnPtr)(const Eigen::VectorXd&) )
//{
//    this->errorFcnPtr = errorFcnPtr;
//}
