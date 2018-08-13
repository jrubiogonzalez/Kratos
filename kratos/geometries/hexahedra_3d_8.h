//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Janosch Stascheit
//                   Felix Nagel
//  contributors:    Hoang Giang Bui
//                   Josep Maria Carbonell
//                   Bodhinanda Chandra
//

#if !defined(KRATOS_HEXAHEDRA_3D_8_H_INCLUDED )
#define  KRATOS_HEXAHEDRA_3D_8_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/quadrilateral_3d_4.h"
#include "integration/hexahedron_gauss_legendre_integration_points.h"

namespace Kratos
{
/**
 * @class Hexahedra3D8
 * @ingroup KratosCore
 * @brief An eight node hexahedra geometry with linear shape functions
 * @details The node ordering corresponds with: 
 *             v
 *      3----------2            
 *      |\     ^   |\          
 *      | \    |   | \        
 *      |  \   |   |  \       
 *      |   7------+---6        
 *      |   |  +-- |-- | -> u   
 *      0---+---\--1   |        
 *       \  |    \  \  |        
 *        \ |     \  \ |         
 *         \|      w  \|         
 *          4----------5   
 * @author Riccardo Rossi
 * @author Janosch Stascheit
 * @author Felix Nagel
 */
template<class TPointType> class Hexahedra3D8 : public Geometry<TPointType>
{
public:
    /**
     * Type Definitions
     */

    /**
     * Geometry as base class.
     */
    typedef Geometry<TPointType> BaseType;

    /**
     * Type of edge and face geometries
     */
    typedef Line3D2<TPointType> EdgeType;
    typedef Quadrilateral3D4<TPointType> FaceType;

    /**
     * Pointer definition of Hexahedra3D8
     */
    KRATOS_CLASS_POINTER_DEFINITION( Hexahedra3D8 );

    /**
     * Integration methods implemented in geometry.
     */
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /**
     * A Vector of counted pointers to Geometries. Used for
     * returning edges of the geometry.
     */
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    /**
     * Redefinition of template parameter TPointType.
     */
    typedef TPointType PointType;

    /**
     * Type used for indexing in geometry class.std::size_t used for indexing
     * point or integration point access methods and also all other
     * methods which need point or integration point index.
     */
    typedef typename BaseType::IndexType IndexType;


    /**
     * This typed used to return size or dimension in
     * geometry. Dimension, WorkingDimension, PointsNumber and
     * ... return this type as their results.
     */
    typedef typename BaseType::SizeType SizeType;

    /**
     * Array of counted pointers to point. This type used to hold
     * geometry's points.
     */
    typedef typename BaseType::PointsArrayType PointsArrayType;

    /**
     * This type used for representing an integration point in
     * geometry. This integration point is a point with an
     * additional weight component.
     */
    typedef typename BaseType::IntegrationPointType IntegrationPointType;

    /**
     * A Vector of IntegrationPointType which used to hold
     * integration points related to an integration
     * method. IntegrationPoints functions used this type to return
     * their results.
     */
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /**
     * A Vector of IntegrationPointsArrayType which used to hold
     * integration points related to different integration method
     * implemented in geometry.
     */
    typedef typename BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;

    /**
     * A third order tensor used as shape functions' values
     * container.
     */
    typedef typename BaseType::ShapeFunctionsValuesContainerType
    ShapeFunctionsValuesContainerType;

    /**
     * A fourth order tensor used as shape functions' local
     * gradients container in geometry.
     */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType
    ShapeFunctionsLocalGradientsContainerType;

    /**
     * A third order tensor to hold jacobian matrices evaluated at
     * integration points. Jacobian and InverseOfJacobian functions
     * return this type as their result.
     */
    typedef typename BaseType::JacobiansType JacobiansType;

    /**
     * A third order tensor to hold shape functions' local
     * gradients. ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /**
     * A third order tensor to hold shape functions' local second derivatives.
     * ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsSecondDerivativesType ShapeFunctionsSecondDerivativesType;

    /**
     * Type of the normal vector used for normal to edges in geomety.
     */
    typedef typename BaseType::NormalType NormalType;

    /**
     * Type of coordinates array
     */
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /**
     * Type of Matrix
     */
    typedef Matrix MatrixType;


    /**
     * Life Cycle
     */
    Hexahedra3D8( const PointType& Point1, const PointType& Point2,
                  const PointType& Point3, const PointType& Point4,
                  const PointType& Point5, const PointType& Point6,
                  const PointType& Point7, const PointType& Point8 )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        this->Points().push_back( typename PointType::Pointer( new PointType( Point1 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point2 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point3 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point4 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point5 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point6 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point7 ) ) );
        this->Points().push_back( typename PointType::Pointer( new PointType( Point8 ) ) );
    }

    Hexahedra3D8( typename PointType::Pointer pPoint1,
                  typename PointType::Pointer pPoint2,
                  typename PointType::Pointer pPoint3,
                  typename PointType::Pointer pPoint4,
                  typename PointType::Pointer pPoint5,
                  typename PointType::Pointer pPoint6,
                  typename PointType::Pointer pPoint7,
                  typename PointType::Pointer pPoint8 )
        : BaseType( PointsArrayType(), &msGeometryData )
    {
        this->Points().push_back( pPoint1 );
        this->Points().push_back( pPoint2 );
        this->Points().push_back( pPoint3 );
        this->Points().push_back( pPoint4 );
        this->Points().push_back( pPoint5 );
        this->Points().push_back( pPoint6 );
        this->Points().push_back( pPoint7 );
        this->Points().push_back( pPoint8 );
    }

    explicit Hexahedra3D8( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        if ( this->PointsNumber() != 8 )
            KRATOS_ERROR << "Invalid points number. Expected 8, given " << this->PointsNumber() << std::endl;
    }

    /**
     * Copy constructor.
     * Construct this geometry as a copy of given geometry.
     *
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    Hexahedra3D8( Hexahedra3D8 const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Copy constructor from a geometry with other point type.
     * Construct this geometry as a copy of given geometry which
     * has different type of points. The given goemetry's
     * TOtherPointType* must be implicity convertible to this
     * geometry PointType.
     *
     * @note This copy constructor don't copy the points and new
     * geometry shares points with given source geometry. It's
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> Hexahedra3D8( Hexahedra3D8<TOtherPointType> const& rOther )
        : BaseType( rOther )
    {
    }

    /// Destructor. Does nothing!!!
    ~Hexahedra3D8() override {}

    GeometryData::KratosGeometryFamily GetGeometryFamily() const override
    {
        return GeometryData::Kratos_Hexahedra;
    }

    GeometryData::KratosGeometryType GetGeometryType() const override
    {
        return GeometryData::Kratos_Hexahedra3D8;
    }

    /**
     * Operators
     */

    /**
     * Assignment operator.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    Hexahedra3D8& operator=( const Hexahedra3D8& rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    /**
     * Assignment operator for geometries with different point type.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    template<class TOtherPointType>
    Hexahedra3D8& operator=( Hexahedra3D8<TOtherPointType> const & rOther )
    {
        BaseType::operator=( rOther );

        return *this;
    }


    /**
     * Operations
     */

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const override
    {
        return typename BaseType::Pointer( new Hexahedra3D8( ThisPoints ) );
    }

    // Geometry< Point<3> >::Pointer Clone() const override
    // {
    //     Geometry< Point<3> >::PointsArrayType NewPoints;

    //     //making a copy of the nodes TO POINTS (not Nodes!!!)
    //     for ( IndexType i = 0 ; i < this->size() ; i++ )
    //     {
    //         NewPoints.push_back(Kratos::make_shared< Point<3> >((*this)[i]));
    //     }

    //     //creating a geometry with the new points
    //     Geometry< Point<3> >::Pointer p_clone( new Hexahedra3D8< Point<3> >( NewPoints ) );

    //     return p_clone;
    // }

    //lumping factors for the calculation of the lumped mass matrix
    Vector& LumpingFactors( Vector& rResult ) const override
    {
        if(rResult.size() != 8)
            rResult.resize( 8, false );
        std::fill( rResult.begin(), rResult.end(), 1.00 / 8.00 );
        return rResult;
    }


    /**
     * Information
     */

    /**
     * This method calculates and returns Length or charactereistic
     * length of this geometry depending on it's dimension. For one
     * dimensional geometry for example Line it returns length of it
     * and for the other geometries it gives Characteristic length
     * otherwise.
     *
     * @return double value contains length or Characteristic
     * length
     * @see Area()
     * @see Volume()
     * @see DomainSize()
     *
     * :TODO: might be necessary to reimplement
     */
    double Length() const override
    {
        return sqrt( fabs( this->DeterminantOfJacobian( PointType() ) ) );
    }

    /**
     * This method calculates and returns area or surface area of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns zero, for two dimensional it gives area
     * and for three dimensional geometries it gives surface area.
     *
     *
     * @return double value contains area or surface area.
     * @see Length()
     * @see Volume()
     * @see DomainSize()
     *
     * :TODO: might be necessary to reimplement
     */
    double Area() const override
    {
        return Volume();
    }

    double Volume() const override //Not a closed formula for a hexahedra
    {
        Vector temp;
        this->DeterminantOfJacobian( temp, msGeometryData.DefaultIntegrationMethod() );
        const IntegrationPointsArrayType& integration_points = this->IntegrationPoints( msGeometryData.DefaultIntegrationMethod() );
        double Volume = 0.00;

        for ( unsigned int i = 0; i < integration_points.size(); i++ )
        {
            Volume += temp[i] * integration_points[i].Weight();
        }

        return Volume;
    }


    /**
     * This method calculate and return length, area or volume of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns its length, for two dimensional it gives area
     * and for three dimensional geometries it gives its volume.
     *
     * @return double value contains length, area or volume.
     * @see Length()
     * @see Area()
     * @see Volume()
     *
     * :TODO: might be necessary to reimplement
     */
    double DomainSize() const override
    {
        return Volume();
    }

    /**
     * Returns a matrix of the local coordinates of all points
     * @param rResult a Matrix that will be overwritten by the results
     * @return the coordinates of all points of the current geometry
     */
    Matrix& PointsLocalCoordinates( Matrix& rResult ) const override
    {
        if ( rResult.size1() != 8 || rResult.size2() != 3 )
            rResult.resize( 8, 3, false );

        rResult( 0, 0 ) = -1.0;
        rResult( 0, 1 ) = -1.0;
        rResult( 0, 2 ) = -1.0;

        rResult( 1, 0 ) = 1.0;
        rResult( 1, 1 ) = -1.0;
        rResult( 1, 2 ) = -1.0;

        rResult( 2, 0 ) = 1.0;
        rResult( 2, 1 ) = 1.0;
        rResult( 2, 2 ) = -1.0;

        rResult( 3, 0 ) = -1.0;
        rResult( 3, 1 ) = 1.0;
        rResult( 3, 2 ) = -1.0;

        rResult( 4, 0 ) = -1.0;
        rResult( 4, 1 ) = -1.0;
        rResult( 4, 2 ) = 1.0;

        rResult( 5, 0 ) = 1.0;
        rResult( 5, 1 ) = -1.0;
        rResult( 5, 2 ) = 1.0;

        rResult( 6, 0 ) = 1.0;
        rResult( 6, 1 ) = 1.0;
        rResult( 6, 2 ) = 1.0;

        rResult( 7, 0 ) = -1.0;
        rResult( 7, 1 ) = 1.0;
        rResult( 7, 2 ) = 1.0;

        return rResult;
    }

    /**
     * Returns whether given arbitrary point is inside the Geometry and the respective 
     * local point for the given global point
     * @param rPoint The point to be checked if is inside o note in global coordinates
     * @param rResult The local coordinates of the point
     * @param Tolerance The  tolerance that will be considered to check if the point is inside or not
     * @return True if the point is inside, false otherwise
     */
    virtual bool IsInside( 
        const CoordinatesArrayType& rPoint, 
        CoordinatesArrayType& rResult, 
        const double Tolerance = std::numeric_limits<double>::epsilon()
        ) override
    {
        this->PointLocalCoordinates( rResult, rPoint );

        if ( std::abs( rResult[0] ) <= (1.0 + Tolerance) )
        {
            if ( std::abs( rResult[1] ) <= (1.0 + Tolerance) )
            {
                if ( std::abs( rResult[2] ) <= (1.0 + Tolerance) )
                {
                    return true;
                }
            }
        }

        return false;
    }


    /** This method gives you number of all edges of this
    geometry.
    @return SizeType containes number of this geometry edges.
    @see Edges()
    @see Edge()
     */
    // will be used by refinement algorithm, thus uncommented. janosch.
    SizeType EdgesNumber() const override
    {
        return 12;
    }

    SizeType FacesNumber() const override
    {
        return 6;
    }

    /** This method gives you all edges of this geometry.

    @return GeometriesArrayType containes this geometry edges.
    @see EdgesNumber()
    @see Edge()
     */
    GeometriesArrayType Edges( void ) override
    {
        GeometriesArrayType edges = GeometriesArrayType();
        typedef typename Geometry<TPointType>::Pointer EdgePointerType;
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 1 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 2 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 3 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 0 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 5 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 6 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 6 ),
                                              this->pGetPoint( 7 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 7 ),
                                              this->pGetPoint( 4 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 4 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 5 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 6 ) ) ) );
        edges.push_back( EdgePointerType( new EdgeType(
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 7 ) ) ) );
        return edges;
    }

    GeometriesArrayType Faces( void ) override
    {
        GeometriesArrayType faces = GeometriesArrayType();
        typedef typename Geometry<TPointType>::Pointer FacePointerType;
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 0 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 1 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 4 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 6 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 1 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 7 ),
                                              this->pGetPoint( 6 ),
                                              this->pGetPoint( 2 ),
                                              this->pGetPoint( 3 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 7 ),
                                              this->pGetPoint( 3 ),
                                              this->pGetPoint( 0 ),
                                              this->pGetPoint( 4 ) ) ) );
        faces.push_back( FacePointerType( new FaceType(
                                              this->pGetPoint( 4 ),
                                              this->pGetPoint( 5 ),
                                              this->pGetPoint( 6 ),
                                              this->pGetPoint( 7 ) ) ) );
        return faces;
    }


    bool HasIntersection( const Point& rLowPoint, const Point& rHighPoint ) override
    {
        using Quadrilateral3D4Type = Quadrilateral3D4<TPointType>;
        // Check if faces have intersection
        if(Quadrilateral3D4Type(this->pGetPoint(0),this->pGetPoint(1), this->pGetPoint(2), this->pGetPoint(3)).HasIntersection(rLowPoint, rHighPoint))
            return true;
        if(Quadrilateral3D4Type(this->pGetPoint(4),this->pGetPoint(5), this->pGetPoint(6), this->pGetPoint(7)).HasIntersection(rLowPoint, rHighPoint))
            return true;
        if(Quadrilateral3D4Type(this->pGetPoint(0),this->pGetPoint(1), this->pGetPoint(5), this->pGetPoint(4)).HasIntersection(rLowPoint, rHighPoint))
            return true;
        if(Quadrilateral3D4Type(this->pGetPoint(3),this->pGetPoint(2), this->pGetPoint(6), this->pGetPoint(7)).HasIntersection(rLowPoint, rHighPoint))
            return true;
        if(Quadrilateral3D4Type(this->pGetPoint(0),this->pGetPoint(4), this->pGetPoint(7), this->pGetPoint(3)).HasIntersection(rLowPoint, rHighPoint))
            return true;
        if(Quadrilateral3D4Type(this->pGetPoint(1),this->pGetPoint(5), this->pGetPoint(6), this->pGetPoint(2)).HasIntersection(rLowPoint, rHighPoint))
            return true;
        
        CoordinatesArrayType local_coordinates;
        // if there are no faces intersecting the box then or the box is inside the hexahedron or it does not have intersection
        if(IsInside(rLowPoint,local_coordinates))
            return true;

        return false;
    }


    /**
     * Shape Function
     */

    /**
     * Calculates the value of a given shape function at a given point.
     *
     * @param ShapeFunctionIndex The number of the desired shape function
     * @param rPoint the given point in local coordinates at which the
     * value of the shape function is calculated
     *
     * @return the value of the shape function at the given point
     * TODO: TO BE VERIFIED
     */
    double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                       const CoordinatesArrayType& rPoint ) const override
    {
        switch ( ShapeFunctionIndex )
        {
        case 0:
            return( 0.125*( 1.0 - rPoint[0] )*( 1.0 - rPoint[1] )*( 1.0 - rPoint[2] ) );
        case 1:
            return( 0.125*( 1.0 + rPoint[0] )*( 1.0 - rPoint[1] )*( 1.0 - rPoint[2] ) );
        case 2:
            return( 0.125*( 1.0 + rPoint[0] )*( 1.0 + rPoint[1] )*( 1.0 - rPoint[2] ) );
        case 3:
            return( 0.125*( 1.0 - rPoint[0] )*( 1.0 + rPoint[1] )*( 1.0 - rPoint[2] ) );
        case 4:
            return( 0.125*( 1.0 - rPoint[0] )*( 1.0 - rPoint[1] )*( 1.0 + rPoint[2] ) );
        case 5:
            return( 0.125*( 1.0 + rPoint[0] )*( 1.0 - rPoint[1] )*( 1.0 + rPoint[2] ) );
        case 6:
            return( 0.125*( 1.0 + rPoint[0] )*( 1.0 + rPoint[1] )*( 1.0 + rPoint[2] ) );
        case 7:
            return( 0.125*( 1.0 - rPoint[0] )*( 1.0 + rPoint[1] )*( 1.0 + rPoint[2] ) );
        default:
            KRATOS_ERROR << "Wrong index of shape function!" << *this << std::endl;
        }

        return 0;
    }

    /** This method gives all non-zero shape functions values
    evaluated at the rCoordinates provided

    @return Vector of values of shape functions \f$ F_{i} \f$
    where i is the shape function index (for NURBS it is the index
    of the local enumeration in the element).

    @see ShapeFunctionValue
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    Vector& ShapeFunctionsValues (Vector &rResult, const CoordinatesArrayType& rCoordinates) const override
    {
      if(rResult.size() != 8) rResult.resize(8,false);
      rResult[0] =  0.125*( 1.0 - rCoordinates[0] )*( 1.0 - rCoordinates[1] )*( 1.0 - rCoordinates[2] ) ;
      rResult[1] =  0.125*( 1.0 + rCoordinates[0] )*( 1.0 - rCoordinates[1] )*( 1.0 - rCoordinates[2] ) ;
      rResult[2] =  0.125*( 1.0 + rCoordinates[0] )*( 1.0 + rCoordinates[1] )*( 1.0 - rCoordinates[2] ) ;
      rResult[3] =  0.125*( 1.0 - rCoordinates[0] )*( 1.0 + rCoordinates[1] )*( 1.0 - rCoordinates[2] ) ;
      rResult[4] =  0.125*( 1.0 - rCoordinates[0] )*( 1.0 - rCoordinates[1] )*( 1.0 + rCoordinates[2] ) ;
      rResult[5] =  0.125*( 1.0 + rCoordinates[0] )*( 1.0 - rCoordinates[1] )*( 1.0 + rCoordinates[2] ) ;
      rResult[6] =  0.125*( 1.0 + rCoordinates[0] )*( 1.0 + rCoordinates[1] )*( 1.0 + rCoordinates[2] ) ;
      rResult[7] =  0.125*( 1.0 - rCoordinates[0] )*( 1.0 + rCoordinates[1] )*( 1.0 + rCoordinates[2] ) ;        
        return rResult;
    }

    /**
     * Calculates the gradients in terms of local coordinateds
     * of all shape functions in a given point.
     *
     * @param rPoint the current point at which the gradients are calculated
     * @return the gradients of all shape functions
     * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
     */
    Matrix& ShapeFunctionsLocalGradients( Matrix& rResult, const CoordinatesArrayType& rPoint ) const override
    {
        if ( rResult.size1() != 8 || rResult.size2() != 3 )
            rResult.resize( 8, 3, false );

        rResult( 0, 0 ) = -0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 - rPoint[2] );
        rResult( 0, 1 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[2] );
        rResult( 0, 2 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[1] );
        rResult( 1, 0 ) =  0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 - rPoint[2] );
        rResult( 1, 1 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[2] );
        rResult( 1, 2 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[1] );
        rResult( 2, 0 ) =  0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 - rPoint[2] );
        rResult( 2, 1 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[2] );
        rResult( 2, 2 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[1] );
        rResult( 3, 0 ) = -0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 - rPoint[2] );
        rResult( 3, 1 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[2] );
        rResult( 3, 2 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[1] );
        rResult( 4, 0 ) = -0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 + rPoint[2] );
        rResult( 4, 1 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[2] );
        rResult( 4, 2 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[1] );
        rResult( 5, 0 ) =  0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 + rPoint[2] );
        rResult( 5, 1 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[2] );
        rResult( 5, 2 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[1] );
        rResult( 6, 0 ) =  0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 + rPoint[2] );
        rResult( 6, 1 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[2] );
        rResult( 6, 2 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[1] );
        rResult( 7, 0 ) = -0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 + rPoint[2] );
        rResult( 7, 1 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[2] );
        rResult( 7, 2 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[1] );

        return rResult;
    }


    /**
     * Calculates the second derivatives in terms of local coordinateds
     * of all shape functions in a given point.
     *
     * @param rResult the given container will be overwritten by the results
     * @param rPoint the given local coordinates the derivatives will be evaluated for.
     * @return a third order tensor containing the second order derivatives for each shape function
     */
    ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives( ShapeFunctionsSecondDerivativesType& rResult, const CoordinatesArrayType& rPoint ) const override
    {
        if ( rResult.size() != this->PointsNumber() )
        {
            // KLUDGE: While there is a bug in
            // ublas vector resize, I have to put this beside resizing!!
            ShapeFunctionsGradientsType temp( this->PointsNumber() );
            rResult.swap( temp );
        }

        for ( unsigned int i = 0; i < this->PointsNumber(); ++i )
        {
            rResult[i].resize(3, 3, false);
        }

        rResult[0]( 0, 0 ) = 0.0;
        rResult[0]( 0, 1 ) = 0.125 * ( 1.0 - rPoint[2] );
        rResult[0]( 0, 2 ) = 0.125 * ( 1.0 - rPoint[1] );
        rResult[0]( 1, 0 ) = 0.125 * ( 1.0 - rPoint[2] );
        rResult[0]( 1, 1 ) = 0.0;
        rResult[0]( 1, 2 ) = 0.125 * ( 1.0 - rPoint[0] );
        rResult[0]( 2, 0 ) = 0.125 * ( 1.0 - rPoint[1] );
        rResult[0]( 2, 1 ) = 0.125 * ( 1.0 - rPoint[0] );
        rResult[0]( 2, 2 ) = 0.0;

        rResult[1]( 0, 0 ) = 0.0;
        rResult[1]( 0, 1 ) = -0.125 * ( 1.0 - rPoint[2] );
        rResult[1]( 0, 2 ) = -0.125 * ( 1.0 - rPoint[1] );
        rResult[1]( 1, 0 ) = -0.125 * ( 1.0 - rPoint[2] );
        rResult[1]( 1, 1 ) = 0.0;
        rResult[1]( 1, 2 ) = 0.125 * ( 1.0 + rPoint[0] );
        rResult[1]( 2, 0 ) = -0.125 * ( 1.0 - rPoint[1] );
        rResult[1]( 2, 1 ) = 0.125 * ( 1.0 + rPoint[0] );
        rResult[1]( 2, 2 ) = 0.0;

        rResult[2]( 0, 0 ) = 0.0;
        rResult[2]( 0, 1 ) = 0.125 * ( 1.0 - rPoint[2] );
        rResult[2]( 0, 2 ) = -0.125 * ( 1.0 + rPoint[1] );
        rResult[2]( 1, 0 ) = 0.125 * ( 1.0 - rPoint[2] );
        rResult[2]( 1, 1 ) = 0.0;
        rResult[2]( 1, 2 ) = -0.125 * ( 1.0 + rPoint[0] );
        rResult[2]( 2, 0 ) = -0.125 * ( 1.0 + rPoint[1] );
        rResult[2]( 2, 1 ) = -0.125 * ( 1.0 + rPoint[0] );
        rResult[2]( 2, 2 ) = 0.0;

        rResult[3]( 0, 0 ) = 0.0;
        rResult[3]( 0, 1 ) = -0.125 * ( 1.0 - rPoint[2] );
        rResult[3]( 0, 2 ) = 0.125 * ( 1.0 + rPoint[1] );
        rResult[3]( 1, 0 ) = -0.125 * ( 1.0 - rPoint[2] );
        rResult[3]( 1, 1 ) = 0.0;
        rResult[3]( 1, 2 ) = -0.125 * ( 1.0 - rPoint[0] );
        rResult[3]( 2, 0 ) = 0.125 * ( 1.0 + rPoint[1] );
        rResult[3]( 2, 1 ) = -0.125 * ( 1.0 - rPoint[0] );
        rResult[3]( 2, 2 ) = 0.0;

        rResult[4]( 0, 0 ) = 0.0;
        rResult[4]( 0, 1 ) = 0.125 * ( 1.0 + rPoint[2] );
        rResult[4]( 0, 2 ) = -0.125 * ( 1.0 - rPoint[1] );
        rResult[4]( 1, 0 ) = 0.125 * ( 1.0 + rPoint[2] );
        rResult[4]( 1, 1 ) = 0.0;
        rResult[4]( 1, 2 ) = -0.125 * ( 1.0 - rPoint[0] );
        rResult[4]( 2, 0 ) = -0.125 * ( 1.0 - rPoint[1] );
        rResult[4]( 2, 1 ) = -0.125 * ( 1.0 - rPoint[0] );
        rResult[4]( 2, 2 ) = 0.0;

        rResult[5]( 0, 0 ) = 0.0;
        rResult[5]( 0, 1 ) = -0.125 * ( 1.0 + rPoint[2] );
        rResult[5]( 0, 2 ) = 0.125 * ( 1.0 - rPoint[1] );
        rResult[5]( 1, 0 ) = -0.125 * ( 1.0 + rPoint[2] );
        rResult[5]( 1, 1 ) = 0.0;
        rResult[5]( 1, 2 ) = -0.125 * ( 1.0 + rPoint[0] );
        rResult[5]( 2, 0 ) = 0.125 * ( 1.0 - rPoint[1] );
        rResult[5]( 2, 1 ) = -0.125 * ( 1.0 + rPoint[0] );
        rResult[5]( 2, 2 ) = 0.0;

        rResult[6]( 0, 0 ) = 0.0;
        rResult[6]( 0, 1 ) = 0.125 * ( 1.0 + rPoint[2] );
        rResult[6]( 0, 2 ) = 0.125 * ( 1.0 + rPoint[1] );
        rResult[6]( 1, 0 ) = 0.125 * ( 1.0 + rPoint[2] );
        rResult[6]( 1, 1 ) = 0.0;
        rResult[6]( 1, 2 ) = 0.125 * ( 1.0 + rPoint[0] );
        rResult[6]( 2, 0 ) = 0.125 * ( 1.0 + rPoint[1] );
        rResult[6]( 2, 1 ) = 0.125 * ( 1.0 + rPoint[0] );
        rResult[6]( 2, 2 ) = 0.0;

        rResult[7]( 0, 0 ) = 0.0;
        rResult[7]( 0, 1 ) = -0.125 * ( 1.0 + rPoint[2] );
        rResult[7]( 0, 2 ) = -0.125 * ( 1.0 + rPoint[1] );
        rResult[7]( 1, 0 ) = -0.125 * ( 1.0 + rPoint[2] );
        rResult[7]( 1, 1 ) = 0.0;
        rResult[7]( 1, 2 ) = 0.125 * ( 1.0 - rPoint[0] );
        rResult[7]( 2, 0 ) = -0.125 * ( 1.0 + rPoint[1] );
        rResult[7]( 2, 1 ) = 0.125 * ( 1.0 - rPoint[0] );
        rResult[7]( 2, 2 ) = 0.0;

        return rResult;
    }


    /**
     * Input and output
     */

    /**
     * Turn back information as a string.
     *
     * @return String contains information about this geometry.
     * @see PrintData()
     * @see PrintInfo()
     */
    std::string Info() const override
    {
        return "3 dimensional hexahedra with eight nodes in 3D space";
    }

    /**
     * Print information about this object.
     *
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    void PrintInfo( std::ostream& rOStream ) const override
    {
        rOStream << "3 dimensional hexahedra with eight nodes in 3D space";
    }

    /**
     * Print geometry's data into given stream.
     * Prints it's points by the order they stored in the geometry
     * and then center point of geometry.
     *
     * @param rOStream Stream to print into it.
     * @see PrintInfo()
     * @see Info()
     */
    void PrintData( std::ostream& rOStream ) const override
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        Matrix jacobian;
        this->Jacobian( jacobian, PointType() );
        rOStream << "    Jacobian in the origin\t : " << jacobian;
    }

protected:

    /**
     * there are no protected class members
     */

private:

    /**
     * Static Member Variables
     */
    static const GeometryData msGeometryData;


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
    }

    Hexahedra3D8(): BaseType( PointsArrayType(), &msGeometryData ) {}

    /**
     * Private Operations
     */

    /**
     * TODO: TO BE VERIFIED
     */
    /**
     * Calculates the gradients in terms of local coordinateds
     * of all shape functions in a given point.
     *
     * @param rPoint the current point at which the gradients are calculated
     * @return the gradients of all shape functions
     * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
     */
    static Matrix ShapeFunctionsLocalGradients( const CoordinatesArrayType& rPoint )
    {
        Matrix result = ZeroMatrix( 8, 3 );
        result( 0, 0 ) = -0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 - rPoint[2] );
        result( 0, 1 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[2] );
        result( 0, 2 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[1] );
        result( 1, 0 ) =  0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 - rPoint[2] );
        result( 1, 1 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[2] );
        result( 1, 2 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[1] );
        result( 2, 0 ) =  0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 - rPoint[2] );
        result( 2, 1 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[2] );
        result( 2, 2 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[1] );
        result( 3, 0 ) = -0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 - rPoint[2] );
        result( 3, 1 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[2] );
        result( 3, 2 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[1] );
        result( 4, 0 ) = -0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 + rPoint[2] );
        result( 4, 1 ) = -0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[2] );
        result( 4, 2 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 - rPoint[1] );
        result( 5, 0 ) =  0.125 * ( 1.0 - rPoint[1] ) * ( 1.0 + rPoint[2] );
        result( 5, 1 ) = -0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[2] );
        result( 5, 2 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 - rPoint[1] );
        result( 6, 0 ) =  0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 + rPoint[2] );
        result( 6, 1 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[2] );
        result( 6, 2 ) =  0.125 * ( 1.0 + rPoint[0] ) * ( 1.0 + rPoint[1] );
        result( 7, 0 ) = -0.125 * ( 1.0 + rPoint[1] ) * ( 1.0 + rPoint[2] );
        result( 7, 1 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[2] );
        result( 7, 2 ) =  0.125 * ( 1.0 - rPoint[0] ) * ( 1.0 + rPoint[1] );
        return result;
    }


    /**
     * TODO: TO BE VERIFIED
     */
    /**
     * Calculates the values of all shape function in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the matrix of values of every shape function in each integration point
     *
     */
    static Matrix CalculateShapeFunctionsIntegrationPointsValues(
        typename BaseType::IntegrationMethod ThisMethod )
    {
        IntegrationPointsContainerType all_integration_points = AllIntegrationPoints();
        IntegrationPointsArrayType& integration_points = all_integration_points[ThisMethod];

        //number of integration points
        const int integration_points_number = integration_points.size();
        //number of nodes in current geometry
        const int points_number = 8;
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, points_number );
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            shape_function_values( pnt, 0 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            shape_function_values( pnt, 1 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            shape_function_values( pnt, 2 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            shape_function_values( pnt, 3 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            shape_function_values( pnt, 4 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
            shape_function_values( pnt, 5 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
            shape_function_values( pnt, 6 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
            shape_function_values( pnt, 7 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
        }

        return shape_function_values;
    }

    /**
     * TODO: TO BE VERIFIED
     */
    /**
     * Calculates the local gradients of all shape functions in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the vector of the gradients of all shape functions
     * in each integration point
     *
     */
    static ShapeFunctionsGradientsType
    CalculateShapeFunctionsIntegrationPointsLocalGradients(
        typename BaseType::IntegrationMethod ThisMethod )
    {
        IntegrationPointsContainerType all_integration_points =
            AllIntegrationPoints();
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values( integration_points_number );
        //initialising container
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Matrix& result = d_shape_f_values[pnt];
            result = ZeroMatrix( 8, 3 );
            result( 0, 0 ) =
                -0.125 * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 0, 1 ) =
                -0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 0, 2 ) =
                -0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() );
            result( 1, 0 ) =
                0.125 * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 1, 1 ) =
                -0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 1, 2 ) =
                -0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() );
            result( 2, 0 ) =
                0.125 * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 2, 1 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 2, 2 ) =
                -0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() );
            result( 3, 0 ) =
                -0.125 * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 3, 1 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Z() );
            result( 3, 2 ) =
                -0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() );
            result( 4, 0 ) =
                -0.125 * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 4, 1 ) =
                -0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 4, 2 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() );
            result( 5, 0 ) =
                0.125 * ( 1.0 - integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 5, 1 ) =
                -0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 5, 2 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 - integration_points[pnt].Y() );
            result( 6, 0 ) =
                0.125 * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 6, 1 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 6, 2 ) =
                0.125 * ( 1.0 + integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() );
            result( 7, 0 ) =
                -0.125 * ( 1.0 + integration_points[pnt].Y() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 7, 1 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Z() );
            result( 7, 2 ) =
                0.125 * ( 1.0 - integration_points[pnt].X() )
                * ( 1.0 + integration_points[pnt].Y() );
        }

        return d_shape_f_values;
    }

    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        IntegrationPointsContainerType integration_points =
        {
            {
                Quadrature < HexahedronGaussLegendreIntegrationPoints1,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < HexahedronGaussLegendreIntegrationPoints2,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < HexahedronGaussLegendreIntegrationPoints3,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < HexahedronGaussLegendreIntegrationPoints4,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature < HexahedronGaussLegendreIntegrationPoints5,
                3, IntegrationPoint<3> >::GenerateIntegrationPoints()
            }
        };
        return integration_points;
    }

    static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        ShapeFunctionsValuesContainerType shape_functions_values =
        {
            {
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_1 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_2 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_3 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_4 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsValues(
                    GeometryData::GI_GAUSS_5 )
            }
        };
        return shape_functions_values;
    }

    /**
     * TODO: TO BE VERIFIED
     */
    static const ShapeFunctionsLocalGradientsContainerType
    AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
        {
            {
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_1 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_2 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_3 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_4 ),
                Hexahedra3D8<TPointType>::CalculateShapeFunctionsIntegrationPointsLocalGradients(
                    GeometryData::GI_GAUSS_5 )
            }
        };
        return shape_functions_local_gradients;
    }


    /**
     * Private Friends
     */

    template<class TOtherPointType> friend class Hexahedra3D8;


    /**
     * Un accessible methods
     */

};// Class Hexahedra3D8


/**
 * Input and output
 */

/**
 * input stream function
 */
template<class TPointType> inline std::istream& operator >> (
    std::istream& rIStream, Hexahedra3D8<TPointType>& rThis );

/**
 * output stream function
 */
template<class TPointType> inline std::ostream& operator << (
    std::ostream& rOStream, const Hexahedra3D8<TPointType>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}


template<class TPointType> const
GeometryData Hexahedra3D8<TPointType>::msGeometryData(
    3, 3, 3, GeometryData::GI_GAUSS_2,
    Hexahedra3D8<TPointType>::AllIntegrationPoints(),
    Hexahedra3D8<TPointType>::AllShapeFunctionsValues(),
    AllShapeFunctionsLocalGradients()
);

}// namespace Kratos.

#endif // KRATOS_HEXAHEDRA_3D_8_H_INCLUDED  defined 
