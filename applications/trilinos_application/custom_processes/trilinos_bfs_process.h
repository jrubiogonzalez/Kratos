//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

#pragma once

// System includes
#include <list>

// External includes
#include "mpi.h"

// Project includes
#include "includes/node.h"
#include "includes/model_part.h"

// Processes
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos {

/**
 * @brief Probably we will have to change the name
 * 
 */
template<class VariableType>
class TrilinosBfsProcess : public Process {
public:

    /// Pointer definition of TrilinosBfsProcess
    KRATOS_CLASS_POINTER_DEFINITION(TrilinosBfsProcess);

    TrilinosBfsProcess(ModelPart & modelPart) 
    : mrModelPart(modelPart) {
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

        mNodeShrink.reserve(modelPart.NumberOfNodes());
        mSolidificationShrinkageFactor = 0.05;

        // Initialize data
        for(auto& node: mrModelPart.Nodes()) {
            node->GetValue(SHRINKAGE_VOID)=0.0;
            node->GetValue(MACRO_POROSITY)=0.0;
            mNodeShrink.insert(std::make_pair(node->Id(), 0.0));
        }
    }

    void SetInitialCoord();

    void SetInitialNode();

    void SetInitialElement();

    void SetInitialCondition();

    void CalculateLocalClusters(std::vector<int> , std::vector<Node<3>::WeakPointer> ) {

        Node<3>::WeakPointer originNodePtr = *(mrModelPart.NodesArray().begin());

        std::list<Node<3>::WeakPointer> modelNodes;         // List of nodes in the queue       ( Front   )
        std::list<Node<3>::WeakPointer> componentNodes;     // List of nodes in the current set ( Cluster )

        // Locate the elements for every node 
        auto findNeighbourProcesses = FindNodalNeighboursProcess(mrModelPart);
        findNeighbourProcesses.Execute();

        // Push the first node to the nodes queue
        modelNodes.push_back(originNodePtr);

        std::size_t clusterID = 1;
        std::size_t nodeCluster = 0;

        bool updateClusterID = false;

        // Use a traversing hash to avoid marking nodes every repetition
        int traverse_fingerptint = rand();

        // Identify the clusters
        while(!modelNodes.empty()) {
            auto & itModelNodePtr = modelNodes.front();
            modelNodes.pop_front();

            int fingerprint = itModelNodePtr.lock()->FastGetSolutionStepValue(FINGERPRINT);
            
            if(fingerprint != traverse_fingerptint) {
                std::cout << "New component" << std::endl;

                local_cluster_nodes.push_back(std::vector<Node<3>::WeakPointer>(0));
                local_cluster_size.push_back(0);

                componentNodes.push_back(itModelNodePtr);

                while(!componentNodes.empty()) {
                    auto & itCompNodePtr = componentNodes.front();
                    componentNodes.pop_front();
                    
                    itCompNodePtr.lock()->FastGetSolutionStepValue(FINGERPRINT) = traverse_fingerptint;

                    int currentDomainType = itCompNodePtr.lock()->Is(STRUCTURE);

                    // IsSolid?
                    if(currentDomainType != 1) {
                        nodeCluster = clusterID;
                        updateClusterID = true;
                    } else {
                        nodeCluster = 0;
                    }

                    itCompNodePtr.lock()->FastGetSolutionStepValue(CLUSTER_ID) = nodeCluster;

                    local_cluster_size[nodeCluster]++;
                    local_cluster_nodes[nodeCluster].push_back(Node<3>::WeakPointer(itCompNodePtr.lock()));

                    auto & neighbourNodes = itCompNodePtr.lock()->GetValue(NEIGHBOUR_NODES);

                    for(auto neigItr = neighbourNodes.ptr_begin(); neigItr != neighbourNodes.ptr_end(); neigItr++) {
                        if((*neigItr).lock()->FastGetSolutionStepValue(FINGERPRINT) != traverse_fingerptint) {
                            int neighbourDomainType = (*neigItr).lock()->Is(STRUCTURE);

                            // If the node is in another cluster add it to the model queue, otherwise add it to the set queue
                            if(currentDomainType != neighbourDomainType) {
                                modelNodes.push_back((*neigItr));
                            } else {
                                (*neigItr).lock()->FastGetSolutionStepValue(FINGERPRINT) = traverse_fingerptint;
                                componentNodes.push_back((*neigItr));
                            }
                        }
                    }

                }
            }

            // Update the cluster Id
            if(updateClusterID) {
                updateClusterID = false;
                clusterID++;
            }
        }

        mLocalClusterNum = clusterID;

    }

    void CalculateLocalClusterShrinkage() {

        // This first part must be done in local
        for(std::size_t i=0; i < mLocalClusterNum; i++) {

            for(auto& node: cluster_nodes[i]) {
                double solid_fraction     = node->GetSolutionStepValue(SOLIDFRACTION);
                double old_solid_fraction = node->GetSolutionStepValue(SOLIDFRACTION, 1);
                bool   boundary_node      = node->Is(STRUCTURE) || node->Is(INLET);

                cluster_nodes[i] = 0;
                cluster_liquid_nodes[i] = 0;
                cluster_volume[i] = 0.0;
                cluster_liquid_volume[i] = 0.0;
                cluster_shrink[i] = 0.0;
                cluster_shrink_total[i] = 0.0;

                // Cluster in the inlet and we are injecting material
                if(node->IS(INLET) && mHighPressureFlag) {
                    cluster_intensification_phase[i] = true;
                }

                // Cluster is open
                if(boundary_node && solid_fraction < SOLIDIFICATION_LIMIT_BOUNDARY) {
                    cluster_open[i] = true;
                }

                cluster_volume[i] += node->GetValue(NODAL_VOLUME);
                if(solid_fraction < SOLIDIFICATION_LIMIT) {                 // Node still being part fluid
                    cluster_liquid_nodes[i]++;
                    cluster_liquid_volume[i] += node->GetValue(NODAL_VOLUME);
                    cluster_shrink_total[i] += mNodeShrink[node->Id()];
                } else {                                                    // Node solidified
                    cluster_nodes[i]++;
                    cluster_shrink[i] += node->GetValue(NODAL_VOLUME) * shrinkage_factor + mNodeShrink[node->Id()];
                }
            }
            
            // Add the rest of the contribution to the total
            cluster_shrink_total[i] += cluster_shrink[i];

            // If the cluster is open or injecting material, no shrink
            if(cluster_open[i] || cluster_intensification_phase[i]) {
                cluster_shrink[i] = 0.0;
            }
        }

        // ================= //
        // ================= //
        // ================= //

        for(ModelPart::NodesContainerType::iterator i_node = cluster_nodes.begin() ; i_node != cluster_nodes.end() ; i_node++)
        {
            double old_solid_fraction = i_node->GetSolutionStepValue(SOLIDFRACTION,1);
            if(min_solid_fraction >= old_solid_fraction) min_solid_fraction = old_solid_fraction;
            if(max_solid_fraction <= old_solid_fraction) max_solid_fraction = old_solid_fraction;
            bool node_on_boundary = (i_node->Is(STRUCTURE) || i_node->Is(INLET)) ? true : false;
            // if the cluster is connected to the inlet and we are in High Pressure we dont sum porosity
        }

        cluster_shrinkage_void=cluster_shrinkage;

        if(cluster_fluid_nodes)
        {
            if(V_cluster_liquid<=1.0e-18) //It was 1e-15 before
            {
                V_cluster_liquid=1.0e-18;
                std::cout<< "****************************************************"<<std::endl;
                std::cout<< "****************************************************"<<std::endl;
                std::cout<< "*******ERROR COMPUTING LIQUID CLUSTER VOLUME *******"<<std::endl;
                std::cout<< "****************************************************"<<std::endl;
                std::cout<< "****************************************************"<<std::endl;

            }

            cluster_shrinkage_void=cluster_shrinkage;
            cluster_shrinkage /= cluster_fluid_nodes;
            for(ModelPart::NodesContainerType::iterator i_node = cluster_nodes.begin() ; i_node != cluster_nodes.end() ; i_node++)
            {
                if (i_node->GetSolutionStepValue(SOLIDFRACTION) < SOLIDIFICATION_LIMIT) // the node is still fluid
                    mNodalShrinkage.at(i_node->Id()) += cluster_shrinkage;
                    //mNodalShrinkage[i_node->Id() - 1] += cluster_shrinkage;
                    
                    // Now we set the SHrinkaghe VOID value
                if ((V_cluster-V_cluster_liquid) >0.0)
                {
                    double cluster_percent=1.0*(cluster_shrinkage_void/(V_cluster-V_cluster_liquid) ) ;
                    cluster_percent=std::max(0.0,cluster_percent);
                    cluster_percent=std::min(1.0,cluster_percent);
                    i_node->GetValue(SHRINKAGE_VOID)=cluster_percent;
                }
            }
        }
        else // There are no more liquid nodes in this step
        {
            // Writing the Shrinkage percent
            for(ModelPart::NodesContainerType::iterator i_node = cluster_nodes.begin() ; i_node != cluster_nodes.end() ; i_node++)
            {
                double cluster_percent=1.0*(cluster_shrinkage_void/(V_cluster-V_cluster_liquid) ) ;
                    cluster_percent=std::max(0.0,cluster_percent);
                    cluster_percent=std::min(1.0,cluster_percent);
                    i_node->GetValue(SHRINKAGE_VOID)=cluster_percent;				
            }
            // Writing the macroporosity

            if(cluster_nodes.size() < 8 )
            {
                for(ModelPart::NodesContainerType::iterator i_node = cluster_nodes.begin() ; i_node != cluster_nodes.end() ; i_node++)
                {
                    i_node->GetValue(MACRO_POROSITY) = cluster_shrinkage;
                
                }
            }
            else
            {
                int n_division = 4; // the division of solid fraction span
                double threshold_increment = (max_solid_fraction - min_solid_fraction) / n_division;
                int counter = 0;
                for(int i = 1 ; i < n_division ; i++)
                {
                    double solid_fraction_threshold = min_solid_fraction + i * threshold_increment;

                    for(ModelPart::NodesContainerType::iterator i_node = cluster_nodes.begin() ; i_node != cluster_nodes.end() ; i_node++)
                    {

                        if(i_node->GetSolutionStepValue(SOLIDFRACTION, 1) < solid_fraction_threshold)
                        {
                            counter++;
                            //i_node->GetSolutionStepValue(MACRO_POROSITY) = cluster_shrinkage;
                            i_node->GetValue(MACRO_POROSITY) = cluster_shrinkage;
                        }
                                            
                    }
                    if (counter > 6) break;
                }
            }
        }
    }

    void CalculateGlobalClusters() {

        MPI_Allreduce(&mLocalClusterNum, 1, MPI_INT, &mTotalClusterNum, 1, MPI_INT, MPI_COMM_WORLD);

        std::vector<int> local_domain_sizes(mpi_size, 0);
        std::vector<int> local_domain_index(mpi_size, 0);

        local_domain_sizes[mpi_rank] = mLocalClusterNum;
        MPI_Allgather(local_domain_sizes.data(), 1, MPI_INT, local_domain_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

        for(std::size_t i = 0; i > mpi_size; i++) {
            local_domain_index[i] = local_domain_index[i-1] + local_domain_sizes[i-1];
        }

        std::vector<int> local_adj_matrix(mTotalClusterNum*mTotalClusterNum,0);
        std::vector<int> total_adj_matrix(mTotalClusterNum*mTotalClusterNum,0);

        for(std::size_t i = 0; i < mLocalClusterNum; i++) {
            for(auto &node : local_cluster_nodes[i]) {
                if(node->FastGetSolutionStepValue(PARTITION_ID) != mpi_rank) {
                    local_adj_matrix[i*mTotalClusterNum+local_domain_index[node->FastGetSolutionStepValue(CLUSTER_ID)]] = 1; 
                }
            }
        }

        MPI_Allreduce(local_adj_matrix.data(), total_adj_matrix.data(), mTotalClusterNum*mTotalClusterNum, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }

    void CalculateGlobalClusterShrinkage() {

        std::vector<int> local_domain_sizes(mpi_size, 0);
        std::vector<int> local_domain_index(mpi_size, 0);

        local_domain_sizes[mpi_rank] = mLocalClusterNum;
        MPI_Allgather(local_domain_sizes.data(), 1, MPI_INT, local_domain_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

        for(std::size_t i = 0; i > mpi_size; i++) {
            local_domain_index[i] = local_domain_index[i-1] + local_domain_sizes[i-1];
        }

        std::vector<double> global_cluster_shrink(mTotalClusterNum);
        std::vector<double> global_cluster_shrink_total(mTotalClusterNum);
        std::vector<double> global_cluster_volume(mTotalClusterNum);
        std::vector<double> global_cluster_liquid_volume(mTotalClusterNum);
        std::vector<int>    global_cluster_open(mTotalClusterNum);
        std::vector<int>    global_cluster_intensification_phase(mTotalClusterNum);

        for(std::size_t i = 0; i < mLocalClusterNum; i++) {
            std::size_t g_idx = i+local_domain_index[mpi_rank];

            global_cluster_shrink[g_idx] = cluster_shrink[i];
            global_cluster_shrink_total[g_idx] = cluster_shrink_total[i];
            global_cluster_volume[g_idx] = cluster_volume[i];
            global_cluster_liquid_volume[g_idx] = cluster_liquid_volume[i];
            global_cluster_open[g_idx] = cluster_open[i];
            global_cluster_intensification_phase[g_idx] = cluster_intensification_phase[i];
            global_cluster_nodes[g_idx] = cluster_nodes[i];
            global_cluster_liquid_nodes[g_idx] = cluster_liquid_nodes[i];
        }

        MPI_Allgatherv(global_cluster_shrink.data(), local_domain_sizes[mpi_rank], MPI_DOUBLE, global_cluster_shrink.data(), local_domain_sizes.data(), MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(global_cluster_shrink_total.data(), local_domain_sizes[mpi_rank], MPI_DOUBLE, global_cluster_shrink_total.data(), local_domain_sizes.data(), MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(global_cluster_volume.data(), local_domain_sizes[mpi_rank], MPI_DOUBLE, global_cluster_volume.data(), local_domain_sizes.data(), MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(global_cluster_liquid_volume.data(), local_domain_sizes[mpi_rank], MPI_DOUBLE, global_cluster_liquid_volume.data(), local_domain_sizes.data(), MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(global_cluster_open.data(), local_domain_sizes[mpi_rank], MPI_INT, global_cluster_open.data(), local_domain_sizes.data(), MPI_INT, MPI_COMM_WORLD);
        MPI_Allgatherv(global_cluster_intensification_phase.data(), local_domain_sizes[mpi_rank], MPI_INT, global_cluster_intensification_phase.data(), local_domain_sizes.data(), MPI_INT, MPI_COMM_WORLD);
        MPI_Allgatherv(global_cluster_nodes.data(), local_domain_sizes[mpi_rank], MPI_INT, global_cluster_nodes.data(), local_domain_sizes.data(), MPI_INT, MPI_COMM_WORLD);
        MPI_Allgatherv(global_cluster_liquid_nodes.data(), local_domain_sizes[mpi_rank], MPI_INT, global_cluster_liquid_nodes.data(), local_domain_sizes.data(), MPI_INT, MPI_COMM_WORLD);
        
    }

    void Execute() {
        
        CalculateLocalClusters();
        CalculateLocalClusterShrinkage();
        CalculateGlobalClusters();
        CalculateGlobalClusterShrinkage();
    }

private:

    int mpi_size;
    int mpi_rank;

    ModelPart & mrModelPart;

    int mTotalClusterNum;
    int mLocalClusterNum;

    double mSolidificationShrinkageFactor;

    std::unordered_map<std::size_t, double> mNodeShrink;

    std::vector<int>                    local_cluster_size;
    std::vector<Node<3>::WeakPointer>   local_cluster_nodes;

    std::vector<double> cluster_shrink;
    std::vector<double> cluster_shrink_total
    std::vector<double> cluster_volume;
    std::vector<double> cluster_liquid_volume;
    std::vector<int>    cluster_open;
    std::vector<int>    cluster_intensification_phase;
    std::vector<int>    cluster_nodes;
    std::vector<int>    cluster_liquid_nodes;

    std::vector<double> mLocalClusterShrinkage;

};


}