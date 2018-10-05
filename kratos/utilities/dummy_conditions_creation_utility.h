//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_DUMMY_COND_UTILS )
#define  KRATOS_DUMMY_COND_UTILS

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/checks.h"

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
 * @class DummyConditionsCreationUtility
 * @ingroup KratosCore
 * @brief This do things
 * @author Vicente Mataix Ferrandiz
 */
class DummyConditionsCreationUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// We create the Pointer related to DummyConditionsCreationUtility
    KRATOS_CLASS_POINTER_DEFINITION(DummyConditionsCreationUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */

    /** Destructor.
     */

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    static void CreateConstraints(ModelPart& rModelPart)
    {
        Matrix rel_mat = ZeroMatrix(3, 3);
        Vector const_vector = ZeroVector(3);
        for (auto& node : rModelPart.Nodes()) {
            ModelPart::DofsVectorType dof_master, dof_slave;
            dof_slave.resize(3);
            dof_slave[0] = node.pGetDof(DISPLACEMENT_X);
            dof_slave[1] = node.pGetDof(DISPLACEMENT_Y);
            dof_slave[2] = node.pGetDof(DISPLACEMENT_Z);
            dof_master.resize(3);
            dof_master[0] = node.pGetDof(NODAL_VAUX_X);
            dof_master[1] = node.pGetDof(NODAL_VAUX_Y);
            dof_master[2] = node.pGetDof(NODAL_VAUX_Z);
            const array_1d<double, 3>& normal = node.FastGetSolutionStepValue(NORMAL);
            const array_1d<double, 3>& tangent_xi = node.FastGetSolutionStepValue(TANGENT_XI);
            const array_1d<double, 3>& tangent_eta = node.FastGetSolutionStepValue(TANGENT_ETA);

            for (int i = 0; i < 3; ++i) {
                rel_mat(i, 0) = normal[i];
                rel_mat(i, 1) = tangent_xi[i];
                rel_mat(i, 2) = tangent_eta[i];
            }

            rModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", node.Id(), dof_master, dof_slave, rel_mat, const_vector);
        }
    }

    ///@}
    ///@name Acces
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Acces
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; /* Class DummyConditionsCreationUtility */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_DUMMY_COND_UTILS  defined */
