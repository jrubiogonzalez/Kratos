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

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "conditions/dummy_dof_condition.h"

namespace Kratos
{
DummyDofCondition::DummyDofCondition(IndexType NewId)
    : BaseType(NewId)
{
}

/***********************************************************************************/
/***********************************************************************************/

DummyDofCondition::DummyDofCondition(
    IndexType NewId, 
    const NodesArrayType& rThisNodes
    ) : BaseType(NewId, rThisNodes)
{
}

/***********************************************************************************/
/***********************************************************************************/

DummyDofCondition::DummyDofCondition(
    IndexType NewId, 
    GeometryType::Pointer pGeometry
    ) : BaseType(NewId, pGeometry)
{
}

/***********************************************************************************/
/***********************************************************************************/

DummyDofCondition::DummyDofCondition(
    IndexType NewId, 
    GeometryType::Pointer pGeometry, 
    PropertiesType::Pointer pProperties
    ) : BaseType(NewId,pGeometry, pProperties)
{
}

/***********************************************************************************/
/***********************************************************************************/

DummyDofCondition::DummyDofCondition(DummyDofCondition const& rOther)
    : BaseType(rOther)
{
}

/***********************************************************************************/
/***********************************************************************************/

DummyDofCondition::~DummyDofCondition()
{
}

/***********************************************************************************/
/***********************************************************************************/

DummyDofCondition& DummyDofCondition::operator=(DummyDofCondition const& rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Condition::operator=(rOther);

    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer DummyDofCondition::Create(
    IndexType NewId, 
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_shared<DummyDofCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer DummyDofCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    KRATOS_TRY
    return Kratos::make_shared<DummyDofCondition>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer DummyDofCondition::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_shared<DummyDofCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void DummyDofCondition::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo
    )
{
    // NEEDED TO DEFINE GLOBAL IDS FOR THE CORRECT ASSEMBLY
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if ( rResult.size() != dimension )
        rResult.resize( dimension );

    rResult[0] = GetGeometry()[0].GetDof( NODAL_VAUX_X ).EquationId();
    rResult[1] = GetGeometry()[0].GetDof( NODAL_VAUX_Y ).EquationId();
    if( dimension == 3)
        rResult[2] = GetGeometry()[0].GetDof( NODAL_VAUX_Z ).EquationId();
}

/***********************************************************************************/
/***********************************************************************************/

void DummyDofCondition::GetDofList(
    DofsVectorType& rConditionDofList,
    ProcessInfo& rCurrentProcessInfo
    )
{
    // NEEDED TO DEFINE THE DOFS OF THE ELEMENT
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if ( rConditionDofList.size() != dimension )
        rConditionDofList.resize( dimension );

    rConditionDofList[0] = GetGeometry()[0].pGetDof( NODAL_VAUX_X );
    rConditionDofList[1] = GetGeometry()[0].pGetDof( NODAL_VAUX_Y );
    if( dimension == 3 )
        rConditionDofList[2] = GetGeometry()[0].pGetDof( NODAL_VAUX_Z );
}

/***********************************************************************************/
/***********************************************************************************/

int DummyDofCondition::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF( this->Id() < 1 ) << "Condition found with Id " << this->Id() << std::endl;

    return 0;

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void DummyDofCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition )
}

void DummyDofCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition )
}

} // Namespace Kratos


