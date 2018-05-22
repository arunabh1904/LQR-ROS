#ifndef SUPERELLIPSOID_KEEPOUT
#define SUPERELLIPSOID_KEEPOUT

#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>

/**
 * @brief The SuperEllipsoidKeepout class
 * This models an object as a bounding superellipsoid, which can be varried from a sphere to a rectange, and calculates
 * the translation and rotation derivatives associated with a set of points being internal to or external to the object.
 * The intent is to use this to get a better bounding box around our robots than just a sphere.  It should enable use to
 * do smarter obstacle avoidance.  It also can be used to keep controls within a certain range, but this might be better
 * handeled otherwise.
 *
 * For more information on superellipoid shapes see: https://en.wikipedia.org/wiki/Superellipsoid
 *
 */
class SuperEllipsoidKeepout
{
public:
    enum TYPE {
                KEEP_OUT,///<- mode returns zero cost if points are outside ellipsoid
                KEEP_IN, ///<- mode returns zero cost if points are inside ellipsoid
                EXACT,    ///<- mode returns zero cost if points are exactly on the boundry
                KEEP_OUT_POTENTIAL_FIELD
              };

    /**
     * @brief SuperEllipsoidKeepout
     * @param xLength The total x length (diameter of circle)
     * @param yLength The total y length (diameter of circle)
     * @param zLength The total z length (diameter of circle)
     * @param r Shape parameter for superellipsoid 2 is equiv to a normal ellipsoid, 4 is more boxy in the xy
     * @param t Shape parameter for superellipsoid 2 is equiv to a normal ellipsoid, 4 is more boxy in the xz,yz plane
     * @param constraint_type KEEP_OUT, KEEP_IN, or EXACT
     */
    SuperEllipsoidKeepout( double xLength, double yLength, double zLength, double r, double t, TYPE constraint_type );

    SuperEllipsoidKeepout( double p_field_mag, double p_field_pow, double xLength, double yLength, double zLength, double r, double t, TYPE constraint_type );

    /**
     * @brief setPointList sets the points to use as constraints
     * @param pointList a set of 3d points that act as targets (KEEP_IN) or obstacles (KEEP_OUT)
     */
    void setPointList( const std::vector<Eigen::Vector3d>& pointList ); ///<- Sets a new list of points for constraints

    /**
     * @brief cost calculates the total cost associated with the point list supplied
     * @param robot_pos_inPointFrame
     * @param rotation_robot_intoPointFrame
     * @return  Cost superellipsoid function -1
     */
    double  cost( const Eigen::Vector3d& robot_pos_inPointFrame, const Eigen::Matrix3d& rotation_robot_intoPointFrame ) const; ///<- The hard constraint. Should be zero when satisfied

    /**
     * @brief cost_jacobian calculates the first derivative of cost with respect to position and rotation
     * @param robot_pos_inPointFrame Position of the robot in the frame the points are described
     * @param rotation_robot_intoPointFrame Orientation of the robot in the frame the poinints are described
     * @return how the cost changes with displacement and angular-rate (think angular velocity)
     *
     * The first three rows are from displacement, the second three rows are from rotation
     */
    Eigen::Matrix<double,6,1> cost_jacobian( const Eigen::Vector3d& robot_pos_inPointFrame, const Eigen::Matrix3d& rotation_robot_intoPointFrame ) const; ///<- the derivative of the constraint with respect to position and orientation (returned as a linear and angular rate)

    /**
     * @brief cost_hessian calculates the second derivative of cost with respect to position and rotation
     * @param robot_pos_inPointFrame  Position of the robot in the frame the points are described
     * @param rotation_robot_intoPointFrame Orientation of the robot in the frame the poinints are described
     * @return how the costJacobian changes with displacement and angular-rate (think angular velocity)
     *
     * The first three columns are from displacement, the second three columns are from rotation
     */
    Eigen::Matrix<double,6,6> cost_hessian(  const Eigen::Vector3d& robot_pos_inPointFrame, const Eigen::Matrix3d& rotation_robot_intoPointFrame ) const; ///<- the second derivative of the constraint with respect to position and orientation (returned as a linear and angular rate)


    bool verify_numerically(bool print_to_screen);

    void setMagnitude( double newMag );
    void setExponent( double newExponent );
    void setSuperEllipseT(double newT );
    void setSuperEllipseR(double newR );

protected:

    double xLen;
    double yLen;
    double zLen;
    double r;
    double t;

    double p_field_mag;
    double n_pow;

    TYPE constraint_type;

    std::vector<Eigen::Vector3d> point_list;

    double calc_cost(double xdif, double ydif, double zdif) const;
    static int sgn(double val);
    static Eigen::Matrix3d cross(const Eigen::Vector3d&);

    double logistic(double x, double k=100) const;
    double logisticDer(double x, double k=100) const;
    double logisitcDer2(double x, double k=100) const;

    static const double logisitcDer2_kx_zero = -2.3993572805154676678327396972823;
    double max_threshold;
    double min_threshold;
    double logistic_K;

};


#endif // SUPERELLIPSOID_KEEPOUT

