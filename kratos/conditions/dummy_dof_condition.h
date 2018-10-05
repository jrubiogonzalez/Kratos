//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//			 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_DUMMY_DOF_CONDITION_H_INCLUDED )
#define  KRATOS_DUMMY_DOF_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/condition.h"

namespace Kratos
{
///@name Kratos Globals
///@{
///@}
///@name Type Definitions
///@{
///@}
///@name  Enum's
///@{
///@}
///@name  Functions
///@{
///@}
///@name Kratos Classes
///@{
/**
 * @class DummyDofCondition
 * @ingroup KratosCore
 * @brief This is conditions adds a dummy dof to the node
 * @author Vicente Mataix Ferrandiz
 */
class DummyDofCondition
    : public Condition
{
public:

    ///@name Type Definitions
    ///@{

    /// We define the base class Condition
    typedef Condition BaseType;
    
    /// Dfinition of the index type
    typedef BaseType::IndexType IndexType;

    /// Definition of the size type
    typedef BaseType::SizeType SizeType;
    
    /// Definition of the node type
    typedef BaseType::NodeType NodeType;

    /// Definition of the properties type
    typedef BaseType::PropertiesType PropertiesType;

    /// Definition of the geometry type with given NodeType
    typedef BaseType::GeometryType GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef BaseType::NodesArrayType NodesArrayType;
    
    /// Counted pointer of DummyDofCondition
    KRATOS_CLASS_POINTER_DEFINITION( DummyDofCondition);
    
    ///@}

public:

    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor.
     * @param NewId The Id of the new created condition
     */
    DummyDofCondition(IndexType NewId = 0);

    /**
     * @brief Constructor using an array of nodes
     * @param NewId The Id of the new created condition
     * @param rThisNodes The array of nodes that will define the geometry that will define the condition
     */
    DummyDofCondition(
        IndexType NewId, 
        const NodesArrayType& rThisNodes
        );

    /**
     * @brief Constructor using Geometry
     * @param NewId The Id of the new created condition
     * @param pGeometry The pointer to the geometry that will define the condition
     */
    DummyDofCondition(
        IndexType NewId, 
        GeometryType::Pointer pGeometry
        );

    /**
     * @brief Constructor using Properties
     * @param NewId The Id of the new created condition
     * @param pGeometry The pointer to the geometry that will define the condition
     * @param pProperties The pointer to the properties that will define the behaviour of the condition
     */
    DummyDofCondition(
        IndexType NewId, 
        GeometryType::Pointer pGeometry, 
        PropertiesType::Pointer pProperties
        );

    ///Copy constructor
    DummyDofCondition(DummyDofCondition const& rOther);

    /// Destructor.
    ~DummyDofCondition() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    DummyDofCondition& operator=(DummyDofCondition const& rOther);

    ///@}
    ///@name Operations
    ///@{
   
    /**
     * @brief Creates a new condition pointer
     * @param NewId the ID of the new condition
     * @param ThisNodes the nodes of the new condition
     * @param pProperties the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(
        IndexType NewId, 
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new condition pointer
     * @param NewId the ID of the new condition
     * @param pGeom the geometry to be employed
     * @param pProperties the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new condition pointer and clones the previous condition data
     * @param NewId the ID of the new condition
     * @param ThisNodes the nodes of the new condition
     * @param pProperties the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone (
        IndexType NewId, 
        NodesArrayType const& ThisNodes
        ) const override;

    /**
     * @brief This determines the condition equation ID vector for all condition DOFs
     * @param rResult the condition equation ID vector
     * @param rCurrentProcessInfo The current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief It determines the condition list of DOFs
     * @param ConditionDofList the list of DOFs
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetDofList(
        DofsVectorType& rConditionDofList,
        ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This method provides the place to perform checks on the completeness of the input and the compatibility with the problem options as well as the contitutive laws selected
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rCurrentProcessInfo The current process info instance
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Contact DoF Condition #" << Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Contact DoF Condition #" << Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
    }
    
    ///@}   

private:

    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;
    
    ///@}

}; // Class DummyDofCondition

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}

} // namespace Kratos.
#endif // KRATOS_DUMMY_DOF_CONDITION_H_INCLUDED  defined
