#ifndef IntegrationMethod_H
#define IntegrationMethod_H

#include <Eigen/Dense>
#include <string>
#include <vector>

enum IntegrationMethodName {
    CONSTRAINT_BASED,
    BACKWORDS_EULER,
    LOBBATO_IIIA_4TH,
    LOBBATO_IIIA_2ND,
    LOBBATO_IIIB_2ND, /**< Not Recomended **/
    LOBBATO_IIIB_4TH, /**< Not Recomended **/
    LOBBATO_IIIC_2ND,
    LOBBATO_IIIC_4TH,
    RANDAU_IA_3RD,    /**< Not Recomended **/
    RANDAU_IA_5TH,
    RANDAU_IIA_3RD,
    RANDAU_IIA_5TH,
    GAUSS_LEGENDRE_2ND,
    GAUSS_LEGENDRE_4TH,
    GAUSS_LEGENDRE_6TH,
    GAUSS_LEGENDRE_8TH,
    GAUSS_LEGENDRE_10TH,

    // NON IMPLICIT METHODS
    DORMAND_PRINCE_45,
    BOGACKI_SHAMPINE_23, /**< Not Recomended For Interpolation **/
    FEHLBERG_45,         /**< Not Recomended For Interpolation **/
    RUNGE_KUTTA_3RD,     /**< Not Recomended For Interpolation **/
    RUNGE_KUTTA_4TH,
    CASH_KARP_45,       /**< Not Recomended For Interpolation **/
    COOPER_VERNER_8TH,   /**< Not Recomended For Interpolation **/
    CURTIS_8TH,
    DORMAN_PRINCE_8TH,
    //DORMAN_PRINCE_12TH,/** doesnot work**/
    RK_10,
    RK_10_8_CURTIS,
    RK_12_10
};

/**
 * @brief The IntegrationMethod class Encodes the Bucher Table for the method. For more information see the <a href="http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods"> Runge Kutta Methods List</a>.
 *
 * This is intended to be used with CollocationNode and ODECollocationSolver.
 *
 */
class IntegrationMethod {

public:
    std::string name() const; /**< Returns the method name **/
    unsigned int numStages() const; /**< Returns the number of internode slopes **/
    unsigned int order() const; /**< Returns the order of the method **/

    Eigen::MatrixXd buildInterpolationMatrix(double h) const;

    double a(unsigned int i, unsigned int j) const; /**< Returnds the a(i,j) slope multipler index **/
    double b(unsigned int i) const; /**< Returns the b(i) index **/
    double b_errEst( unsigned int i) const; /**< Returns the _b2(i) index if the method is adaptive **/
    double c(unsigned int i) const; /**< Returns the c(i) index **/

    unsigned int c_duplication_index( unsigned int i) const; /** returns the most accurate state index for interpolation of predictor/corrector methods **/


   explicit IntegrationMethod(IntegrationMethodName name = LOBBATO_IIIA_4TH );

   enum integratorType{ EXPLICIT, IMPLICIT, NONE };
   integratorType type() const; /**< Returns the integration type **/

   bool first_same_as_last() const;
   bool isAdaptive() const;


protected:
    Eigen::MatrixXd _a;
    Eigen::VectorXd _b;
    Eigen::VectorXd _b2;
    Eigen::VectorXd _c;

    std::vector<unsigned int> c_duplicate;

    unsigned int _order;
    unsigned int _numStages;
    std::string _name;
    bool _adaptive_method;
    bool _first_same_as_last;

    integratorType _type; /**< The integrator type **/

};



// The predefined methods
// IMPLICIT METHODS
extern const IntegrationMethod BackwordEuler;
extern const IntegrationMethod Lobbato_IIIA_4th;
extern const IntegrationMethod Lobbato_IIIA_2nd;
extern const IntegrationMethod Lobbato_IIIB_4th;
extern const IntegrationMethod Lobbato_IIIB_2nd;
extern const IntegrationMethod Lobbato_IIIC_4th;
extern const IntegrationMethod Lobbato_IIIC_2nd;
extern const IntegrationMethod Randau_IA_3rd;
extern const IntegrationMethod Randau_IA_5th;
extern const IntegrationMethod Randau_IIA_3rd;
extern const IntegrationMethod Randau_IIA_5th;
extern const IntegrationMethod Gauss_Legendre_2nd;
extern const IntegrationMethod Gauss_Legendre_4th;
extern const IntegrationMethod Gauss_Legendre_6th;
extern const IntegrationMethod Gauss_Legendre_8th;
extern const IntegrationMethod Gauss_Legendre_10th;


// EXPLICIT METHODS
extern const IntegrationMethod Bogacki_Shampine_3rd;
extern const IntegrationMethod Runge_Kutta_3rd;
extern const IntegrationMethod Runge_Kutta_4th;
extern const IntegrationMethod Dormand_Prince_45; /**< Not Recomended For Interpolation **/
extern const IntegrationMethod Fehlberg_45;
extern const IntegrationMethod Cash_Karp_45;           /**< Not Recomended For Interpolation **/
extern const IntegrationMethod Cooper_Verner_8th; /**< Not Recomended For Interpolation **/
extern const IntegrationMethod Curtis_8th; /**< Not Great For Interpolation **/
extern const IntegrationMethod Dormand_Prince_8th;
extern const IntegrationMethod Rk_10;/**< Not Recomended For Interpolation **/
extern const IntegrationMethod Rk_10_8_Curtis;
//static const IntegrationMethod Rk_12_10(RK_12_10); /**< Not Recomended For Interpolation **/


// A CONSTAINT DRIVEN NODE NO INTEGRATION PERFORMED
extern const IntegrationMethod ConstraintBased;




#endif // IntegrationMethod_H
