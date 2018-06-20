//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#if !defined(KRATOS_SHALLOW_WATER_VARIABLES_UTILITY_H_INCLUDED)
#define  KRATOS_SHALLOW_WATER_VARIABLES_UTILITY_H_INCLUDED

// System includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "processes/find_nodal_neighbours_process.h"

// External includes 

// Project includes
#include "shallow_water_application.h"

namespace Kratos
{

typedef std::size_t IndexType;

typedef std::vector<IndexType> IndexVectorType;

typedef std::map<int, Properties::Pointer> PropertiesMapType;

class ShallowWaterVariablesUtility
{
public:

    KRATOS_CLASS_POINTER_DEFINITION(ShallowWaterVariablesUtility);

    ShallowWaterVariablesUtility(ModelPart& rModelPart, const double& rDryHeight = 1e-3)
        : mrModelPart(rModelPart)  
    {
        KRATOS_TRY
        
        std::cout << "Initializing shallow water variables utility" << std::endl; 
        mWaterHeightConvert = mrModelPart.GetProcessInfo()[WATER_HEIGHT_UNIT_CONVERTER];
        mDryHeight = rDryHeight;
        mZeroValue = 1e-8;
        
        KRATOS_CATCH("")
    }

    ~ShallowWaterVariablesUtility()
    {}

    /**
     * This method computes the free surface elevation as HEIGHT + BATHYMETRY
     */
    void ComputeFreeSurfaceElevation()
    {
        KRATOS_TRY

        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_nodes.size()); i++)
        {
            ModelPart::NodesContainerType::iterator inode = r_nodes.begin() + i;
            inode->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) = inode->FastGetSolutionStepValue(HEIGHT) + (inode->FastGetSolutionStepValue(BATHYMETRY) / mWaterHeightConvert);
        }

        KRATOS_CATCH("")
    }

    /**
     * This method computes the water height as FREE_SURFACE_ELEVATION - BATHYMETRY
     */
    void ComputeHeightFromFreeSurface()
    {
        KRATOS_TRY

        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_nodes.size()); i++)
        {
            ModelPart::NodesContainerType::iterator inode = r_nodes.begin() + i;
            inode->FastGetSolutionStepValue(HEIGHT) = inode->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION) - (inode->FastGetSolutionStepValue(BATHYMETRY) / mWaterHeightConvert);
        }

        KRATOS_CATCH("")
    }

    /** 
     * This method computes the velocity as the MOMENTUM / HEIGHT
     */
    void ComputeVelocity()
    {
        KRATOS_TRY

        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_nodes.size()); i++)
        {
            ModelPart::NodesContainerType::iterator inode = r_nodes.begin() + i;
            inode->GetSolutionStepValue(VELOCITY) = inode->FastGetSolutionStepValue(MOMENTUM) / (inode->FastGetSolutionStepValue(HEIGHT) * mWaterHeightConvert);
        }

        KRATOS_CATCH("")
    }

    /** 
     * This method computes the momentum as the VELOCITY * HEIGHT
     */
    void ComputeMomentum()
    {
        KRATOS_TRY

        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_nodes.size()); i++)
        {
            ModelPart::NodesContainerType::iterator inode = r_nodes.begin() + i;
            inode->GetSolutionStepValue(MOMENTUM) = inode->FastGetSolutionStepValue(VELOCITY) * (inode->FastGetSolutionStepValue(HEIGHT) * mWaterHeightConvert);
        }

        KRATOS_CATCH("")
    }

    /**
     * Set the dry nodes to zero
     */
    void CheckDryConservedVariables()
    {
        KRATOS_TRY

        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_nodes.size()); i++)
        {
            ModelPart::NodesContainerType::iterator inode = r_nodes.begin() + i;
            if (inode->FastGetSolutionStepValue(HEIGHT) < mDryHeight &&
                inode->FastGetSolutionStepValue(RAIN)   < mDryHeight )
            {
                inode->FastGetSolutionStepValue(HEIGHT)     = mZeroValue;
                inode->FastGetSolutionStepValue(MOMENTUM_X) = 0;
                inode->FastGetSolutionStepValue(MOMENTUM_Y) = 0;
            }
        }
        KRATOS_CATCH("")
    }

    /**
     * Set the dry nodes to zero
     */
    void CheckDryPrimitiveVariables()
    {
        KRATOS_TRY

        ModelPart::NodesContainerType& r_nodes = mrModelPart.Nodes();
        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_nodes.size()); i++)
        {
            ModelPart::NodesContainerType::iterator inode = r_nodes.begin() + i;
            if (inode->FastGetSolutionStepValue(HEIGHT) < mDryHeight &&
                inode->FastGetSolutionStepValue(RAIN)   < mDryHeight )
            {
                inode->FastGetSolutionStepValue(HEIGHT)     = mZeroValue;
                inode->FastGetSolutionStepValue(VELOCITY_X) = 0;
                inode->FastGetSolutionStepValue(VELOCITY_Y) = 0;
            }
        }
        KRATOS_CATCH("")
    }

    /**
     * This method looks for the dry elements and sets the proper flag
     */
    void SetDryWetState()
    {
        KRATOS_TRY
        
        // Getting the elements from the model
        const int nelem = static_cast<int>(mrModelPart.Elements().size());
        int nodes;
        bool wet_node;
        
        // And now, if an element has all nodes dry, it is not active
        #pragma omp parallel for
        for(int k = 0; k < nelem; k++)
        {
            ModelPart::ElementsContainerType::iterator it = mrModelPart.ElementsBegin() + k;
            nodes = it->GetGeometry().size();
            wet_node = false;
            for(int l = 0; l < nodes; l++)
            {
                if (it->GetGeometry()[l].FastGetSolutionStepValue(HEIGHT) >= mDryHeight ||
                    it->GetGeometry()[l].FastGetSolutionStepValue(RAIN)   >= mDryHeight )
                    wet_node = true;  // It means there is almost a wet node
            }

            if (wet_node)
            {
                it->Set(FLUID, true);
                it->Set(ACTIVE, true);
            }
            else
            {
                it->Set(FLUID, false);
                it->Set(ACTIVE, false);
            }
        }
        
        KRATOS_CATCH("")
    }

    /**
     * This method creates the dry properties as a copy of the wet properties
     * The only difference between them is for visualization purpose
     */
    void DefineDryProperties()
    {
        // Create a copy for each property
        const int nprop = static_cast<int>(mrModelPart.NumberOfProperties());
        ModelPart::PropertiesContainerType::iterator prop_begin = mrModelPart.PropertiesBegin();
        IndexType last_id = 0;
        IndexVectorType prop_id;

        for (int i = 0; i < nprop; i++)
        {
            ModelPart::PropertiesContainerType::iterator prop = prop_begin + i;

            if (prop->Id() > last_id)
                last_id = prop->Id();
            prop_id.push_back(prop->Id());
        }

        // for (int i = 0; i < nprop; i++)
        // {
        //     ModelPart::PropertiesContainerType::iterator prop = prop_begin + i;
        for (auto id : prop_id)
        {
            // Get pointers to the properties and create the dry property
            Properties::Pointer wet_prop = mrModelPart.pGetProperties(id); // This work around is inefficient. TODO: find another way
            Properties::Pointer dry_prop(new Properties(*wet_prop));
            dry_prop->SetId(++last_id);

            // Add the new property and add them to the maps
            mrModelPart.AddProperties(dry_prop, dry_prop->Id());
            mWetToDryPropertiesMap[wet_prop->Id()] = dry_prop;
            mDryToWetPropertiesMap[dry_prop->Id()] = wet_prop;
        }
    }

    /**
     * This method assign the wet and dry properties
     * Wet and dry are tween properties 
     * The only difference between them is for visualization purpose
     * ExecuteBeforOutputStep
     * @see DefineDryProperties
     */
    void AssignDryWetProperties()
    {
        const int nelem = static_cast<int>(mrModelPart.Elements().size());
        ModelPart::ElementsContainerType::iterator elem_begin = mrModelPart.ElementsBegin();

        #pragma omp parallel for
        for (int i = 0; i < nelem; i++)
        {
            ModelPart::ElementsContainerType::iterator elem = elem_begin + i;

            if (elem->Is(FLUID))
            {
                auto search = mDryToWetPropertiesMap.find(elem->GetProperties().Id());
                if (search != mDryToWetPropertiesMap.end()) // The element was dry
                    elem->SetProperties(search->second);
            }
            else
            {
                auto search = mWetToDryPropertiesMap.find(elem->GetProperties().Id());
                if (search != mWetToDryPropertiesMap.end()) // The element was wet
                {
                    elem->SetProperties(search->second);
                }
            }
        }
    }

    /**
     * This method sets the mesh position for visualization purpose
     * ExecuteBeforeOutputStep
     * @see ResetMeshPosition
     */
    void SetMeshPosition()
    {
        // Move mesh to the current position
        const int nodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = mrModelPart.NodesBegin();

        #pragma omp parallel for
        for(int i = 0; i < nodes; i++)
        {
            ModelPart::NodesContainerType::iterator node = node_begin + i;

            if (node->FastGetSolutionStepValue(HEIGHT) <= mDryHeight)
            {
                double value = node->FastGetSolutionStepValue(BATHYMETRY);
                node->Z() = value;
                node->FastGetSolutionStepValue(DISPLACEMENT_Z) = value;
            }
            else
            {
                double value = node->FastGetSolutionStepValue(FREE_SURFACE_ELEVATION);
                node->Z() = value;
                node->FastGetSolutionStepValue(DISPLACEMENT_Z) = value;
            }
        }
    }

    /**
     * This method resets the mesh to the original position (Z0 = 0)
     * ExecuteAfterOutputStep
     * @see SetMeshPosition
     */
    void ResetMeshPosition()
    {
        // Move mesh to the original position
        const int nodes = static_cast<int>(mrModelPart.Nodes().size());
        ModelPart::NodesContainerType::iterator node_begin = mrModelPart.NodesBegin();

        #pragma omp parallel for
        for(int i = 0; i < nodes; i++)
        {
            ModelPart::NodesContainerType::iterator node = node_begin + i;

            node->Z() = node->Z0();
        }
    }

    /**
     * This method sets the all the elements active for visualization purpose
     * ExecuteBeforeOutputStep
     * @see AssignDryWetProperties
     * @see SetDryWetState
     */
    void SetElementsActive()
    {
        const int nelem = static_cast<int>(mrModelPart.Elements().size());
        ModelPart::ElementsContainerType::iterator elem_begin = mrModelPart.ElementsBegin();

        #pragma omp parallel for
        for (int i = 0; i < nelem; i++)
        {
            ModelPart::ElementsContainerType::iterator elem = elem_begin + i;

            elem->Set(ACTIVE, true);
        }
    }

protected:

private:

    ModelPart& mrModelPart;
    double mWaterHeightConvert;
    double mDryHeight;
    double mZeroValue;
    PropertiesMapType mWetToDryPropertiesMap;
    PropertiesMapType mDryToWetPropertiesMap;

}; // class ShallowWaterVariablesUtility

} // namespace Kratos

# endif // KRATOS_SHALLOW_WATER_VARIABLES_UTILITY_H_INCLUDED
