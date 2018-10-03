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
#include <vector>

// External includes
#include "mpi.h"

// Project includes
#include "includes/node.h"
#include "includes/model_part.h"
#include "includes/c2c_variables.h"

// Processes
#include "processes/find_nodal_neighbours_process.h"

namespace Kratos {


/**
 * @brief auxiliary class to represent the global graph of
 * connections between partitions
 * 
 */
struct GraphNode {
    std::size_t mId;
    std::size_t mClusterId;
    bool mVisited;
    std::vector<std::shared_ptr<GraphNode>> mNeighbours;

    GraphNode(std::size_t Id) : mId(Id), mVisited(false) {}
};

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
            // node->GetValue(SHRINKAGE_VOID)=0.0;
            node.GetValue(MACRO_POROSITY)=0.0;
            mNodeShrink.insert(std::make_pair(node.Id(), 0.0));
        }
    }

    void SetInitialCoord();

    void SetInitialNode();

    void SetInitialElement();

    void SetInitialCondition();

    void CalculateLocalClusters(std::vector<int> & cluster_size, std::vector<std::vector<Node<3>::WeakPointer>> & cluster_nodes) {

        std::list<Node<3>::WeakPointer> modelNodes;         // List of nodes in the queue       ( Front   )
        std::list<Node<3>::WeakPointer> componentNodes;     // List of nodes in the current set ( Cluster )

        // Locate the elements for every node 
        auto findNeighbourProcesses = FindNodalNeighboursProcess(mrModelPart);
        findNeighbourProcesses.Execute();

        // Push the first node to the nodes queue

        std::size_t clusterID = 1;
        std::size_t nodeCluster = 0;

        bool updateClusterID = false;

        // Use a traversing hash to avoid marking nodes every repetition
        int traverse_fingerptint = rand();

        // Initialize the local clusters
        cluster_nodes.push_back(std::vector<Node<3>::WeakPointer>(0));
        cluster_size.push_back(0);

        Node<3>::WeakPointer originNodePtr = *(mrModelPart.NodesArray().begin());
        modelNodes.push_back(originNodePtr);

        // Identify the clusters
        while(!modelNodes.empty()) {
            auto & itModelNodePtr = modelNodes.front();
            modelNodes.pop_front();

            int fingerprint = itModelNodePtr.lock()->FastGetSolutionStepValue(FINGERPRINT);
            
            if(fingerprint != traverse_fingerptint) {
                std::cout << "New component" << std::endl;

                cluster_nodes.push_back(std::vector<Node<3>::WeakPointer>(0));
                cluster_size.push_back(0);

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

                    cluster_size[nodeCluster]++;
                    cluster_nodes[nodeCluster].push_back(Node<3>::WeakPointer(itCompNodePtr.lock()));

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

    /**
     * @brief Calculate cluster shringkage. The first part can be done
     * distributetly without communication
     * 
     * @param cluster_nodes Vector of Vector of nodes containing the list of nodes of eavery
     * local cluster identified.
     */
    void CalculateLocalClusterShrinkage(std::vector<std::vector<Node<3>::WeakPointer>> & cluster_nodes) {

        double SOLIDIFICATION_LIMIT_BOUNDARY = 0.0;
        double SOLIDIFICATION_LIMIT = 0.0;
        double shrinkage_factor = 1.0;

        // Store the statistics of every cluster. We cannot do cluster by cluster
        // as in a later stage some clusters will need to merge.
        cluster_shrink.resize(mLocalClusterNum);
        cluster_shrink_total.resize(mLocalClusterNum);
        cluster_shrink_void.resize(mLocalClusterNum);
        cluster_volume.resize(mLocalClusterNum);
        cluster_liquid_volume.resize(mLocalClusterNum);
        cluster_min_solid_fraction.resize(mLocalClusterNum);
        cluster_max_solid_fraction.resize(mLocalClusterNum);
        cluster_open.resize(mLocalClusterNum);
        cluster_intensification_phase.resize(mLocalClusterNum);
        cluster_nodes_count.resize(mLocalClusterNum);
        cluster_liquid_nodes_count.resize(mLocalClusterNum);

        // This first part must be done in local
        for(std::size_t i=0; i < mLocalClusterNum; i++) {

            for(auto& node_wk_ptr: cluster_nodes[i]) {
                auto node = node_wk_ptr.lock();

                double solid_fraction     = node->GetSolutionStepValue(SOLIDFRACTION);
                double old_solid_fraction = node->GetSolutionStepValue(SOLIDFRACTION, 1);
                bool   boundary_node      = node->Is(STRUCTURE) || node->Is(INLET);

                cluster_volume[i] = 0.0;
                cluster_liquid_volume[i] = 0.0;
                cluster_min_solid_fraction[i] = std::numeric_limits<double>::max();
                cluster_max_solid_fraction[i] = std::numeric_limits<double>::min();
                cluster_shrink[i] = 0.0;
                cluster_shrink_total[i] = 0.0;
                cluster_shrink_void[i] = 0.0;
                cluster_nodes_count[i] = 0;
                cluster_liquid_nodes_count[i] = 0;

                // Cluster connected to the inlet and we are injecting material
                if(node->Is(INLET) && mHighPressureFlag) {
                    cluster_intensification_phase[i] = true;
                }

                // Is cluster open?
                if(boundary_node && solid_fraction < SOLIDIFICATION_LIMIT_BOUNDARY) {
                    cluster_open[i] = true;
                }

                // Add the volume contribution
                cluster_volume[i] += node->GetValue(NODAL_VOLUME);

                if(solid_fraction < SOLIDIFICATION_LIMIT) {                 // Node still liquid
                    cluster_liquid_nodes_count[i]++;
                    cluster_liquid_volume[i] += node->GetValue(NODAL_VOLUME);
                    cluster_shrink_total[i] += mNodeShrink[node->Id()];
                } else {                                                    // Node solidified
                    cluster_nodes_count[i]++;
                    cluster_shrink[i] += node->GetValue(NODAL_VOLUME) * shrinkage_factor + mNodeShrink[node->Id()];
                }
            }
            
            // Add the rest of the contribution to the total
            cluster_shrink_total[i] += cluster_shrink[i];

            // If the cluster is open or in intensification phase, no shrinkage
            if(cluster_open[i] || cluster_intensification_phase[i]) {
                cluster_shrink[i] = 0.0;
            }

            cluster_shrink_void[i] = cluster_shrink[i];

            // End of the local part. Next stages need the real clusters.
        }
    }

    void CalculateGlobalClusters(std::vector<int> & local_domain_sizes, std::vector<int> & local_domain_index, std::vector<int> & local_adj_matrix, std::vector<int> & total_adj_matrix, std::vector<std::shared_ptr<GraphNode>> &global_adj_graph) {

        mrModelPart.GetCommunicator().SynchronizeNodalSolutionStepsData();

        for(std::size_t i = 1; i < mLocalClusterNum; i++) {
            for(auto &node : local_cluster_nodes[i]) {
                if(node.lock()->FastGetSolutionStepValue(PARTITION_INDEX) != mpi_rank && node.lock()->FastGetSolutionStepValue(CLUSTER_ID) != 0) {

                    std::size_t idx = i + local_domain_index[mpi_rank];
                    std::size_t idr = local_domain_index[node.lock()->FastGetSolutionStepValue(PARTITION_INDEX)]+node.lock()->FastGetSolutionStepValue(CLUSTER_ID);

                    std::cout << "Cluster " << i << " from partition " << mpi_rank << " " << local_domain_index[node.lock()->FastGetSolutionStepValue(PARTITION_INDEX)] << " is neighbouring " <<  node.lock()->FastGetSolutionStepValue(CLUSTER_ID) << " from partition " << node.lock()->FastGetSolutionStepValue(PARTITION_INDEX) << std::endl;
                    local_adj_matrix[idx*mTotalClusterNum+idr] = 1;
                    local_adj_matrix[idr*mTotalClusterNum+idx] = 1; // if neighb(a,b) = 1, then neighb(b,a) = 1 as well
                }
            }
        }

        // // Debug: Print the global adj matrix:
        // for(std::size_t p = 0; p < mpi_size; p++) {
        //     MPI_Barrier(MPI_COMM_WORLD);
        //     if(mpi_rank == p) {
        //         for(std::size_t i = 0; i < mTotalClusterNum; i++) {
        //             for(std::size_t j = 0; j < mTotalClusterNum; j++) {
        //                 std::cout << local_adj_matrix[i*mTotalClusterNum+j] << " ";
        //             }
        //             std::cout << std::endl;
        //         }
        //         std::cout << "=====" << std::endl;
        //     }
        // }

        MPI_Allreduce(local_adj_matrix.data(), total_adj_matrix.data(), mTotalClusterNum*mTotalClusterNum, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        // Debug: Print the global adj matrix:
        if(!mpi_rank) {
            for(std::size_t i = 0; i < mTotalClusterNum; i++) {
                for(std::size_t j = 0; j < mTotalClusterNum; j++) {
                    std::cout << total_adj_matrix[mTotalClusterNum*i+j] << " ";
                }
                std::cout << std::endl;
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // Create the graph
        for(std::size_t i = 0; i < mTotalClusterNum; i++) {
            global_adj_graph.push_back(std::shared_ptr<GraphNode>(new GraphNode(i)));
        }

        // Populate the connections
        for(std::size_t i = 0; i < mTotalClusterNum; i++) {
            for(std::size_t j = 0; j < mTotalClusterNum; j++) {
                if(total_adj_matrix[mTotalClusterNum*i+j] != 0 && i != j) {
                    global_adj_graph[i]->mNeighbours.push_back(global_adj_graph[j]);
                }
            }
        }
    }

    void CalculateGlobalClusterShrinkage(std::vector<int> & local_domain_sizes, std::vector<int> & local_domain_index, std::vector<int> & total_adj_matrix, std::vector<std::shared_ptr<GraphNode>> &global_adj_graph) {

        double SOLIDIFICATION_LIMIT_BOUNDARY = 0.0;
        double SOLIDIFICATION_LIMIT = 0.0;
        double shrinkage_factor = 1.0;
        
        std::vector<double> global_cluster_shrink(mTotalClusterNum);
        std::vector<double> global_cluster_shrink_total(mTotalClusterNum);
        std::vector<double> global_cluster_shrink_void(mTotalClusterNum);
        std::vector<double> global_cluster_volume(mTotalClusterNum);
        std::vector<double> global_cluster_liquid_volume(mTotalClusterNum);
        std::vector<int>    global_cluster_open(mTotalClusterNum);
        std::vector<int>    global_cluster_intensification_phase(mTotalClusterNum);
        std::vector<int>    global_cluster_nodes_count(mTotalClusterNum);
        std::vector<int>    global_cluster_liquid_nodes_count(mTotalClusterNum);

        for(std::size_t i = 0; i < mLocalClusterNum; i++) {
            std::size_t g_idx = i+local_domain_index[mpi_rank];

            global_cluster_shrink[g_idx] = cluster_shrink[i];
            global_cluster_shrink_total[g_idx] = cluster_shrink_total[i];
            global_cluster_shrink_void[g_idx] = cluster_shrink_void[i];
            global_cluster_volume[g_idx] = cluster_volume[i];
            global_cluster_liquid_volume[g_idx] = cluster_liquid_volume[i];
            global_cluster_open[g_idx] = cluster_open[i];
            global_cluster_intensification_phase[g_idx] = cluster_intensification_phase[i];
            global_cluster_nodes_count[g_idx] = cluster_nodes_count[i];
            global_cluster_liquid_nodes_count[g_idx] = cluster_liquid_nodes_count[i];
        }

        MPI_Allgatherv(global_cluster_shrink.data(), local_domain_sizes[mpi_rank], MPI_DOUBLE, global_cluster_shrink.data(), local_domain_sizes.data(), local_domain_index.data(), MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(global_cluster_shrink_total.data(), local_domain_sizes[mpi_rank], MPI_DOUBLE, global_cluster_shrink_total.data(), local_domain_sizes.data(), local_domain_index.data(), MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(global_cluster_volume.data(), local_domain_sizes[mpi_rank], MPI_DOUBLE, global_cluster_volume.data(), local_domain_sizes.data(), local_domain_index.data(), MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(global_cluster_liquid_volume.data(), local_domain_sizes[mpi_rank], MPI_DOUBLE, global_cluster_liquid_volume.data(), local_domain_sizes.data(), local_domain_index.data(), MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgatherv(global_cluster_open.data(), local_domain_sizes[mpi_rank], MPI_INT, global_cluster_open.data(), local_domain_sizes.data(), local_domain_index.data(), MPI_INT, MPI_COMM_WORLD);
        MPI_Allgatherv(global_cluster_intensification_phase.data(), local_domain_sizes[mpi_rank], MPI_INT, global_cluster_intensification_phase.data(), local_domain_sizes.data(), local_domain_index.data(), MPI_INT, MPI_COMM_WORLD);
        MPI_Allgatherv(global_cluster_nodes_count.data(), local_domain_sizes[mpi_rank], MPI_INT, global_cluster_nodes_count.data(), local_domain_sizes.data(), local_domain_index.data(), MPI_INT, MPI_COMM_WORLD);
        MPI_Allgatherv(global_cluster_liquid_nodes_count.data(), local_domain_sizes[mpi_rank], MPI_INT, global_cluster_liquid_nodes_count.data(), local_domain_sizes.data(), local_domain_index.data(), MPI_INT, MPI_COMM_WORLD);

        // At this point:
        //  - All partitions have the information of all clusters.
        //  - All partitions have the global cluster graph

        // Calculate the True clusters
        std::list<std::shared_ptr<GraphNode>> frontNodes;  // List of nodes in the queue ( Front )
        std::size_t true_cluster_id = 0;
        std::vector<std::list<std::size_t>> true_cluster_comps_id;
        std::vector<int> true_cluster_comps_rev(mTotalClusterNum);

        for(std::size_t i = 0; i < mTotalClusterNum; i++) {
            if(!global_adj_graph[i]->mVisited) {
                true_cluster_id++;
                if(mpi_rank == 0)
                    std::cout << "New true cluster " << true_cluster_id << std::endl;
                frontNodes.push_back(global_adj_graph[i]);
                true_cluster_comps_id.push_back(std::list<std::size_t>(0));
            }

            while(!frontNodes.empty()) {
                auto thisNode = frontNodes.front();
                frontNodes.pop_front();

                thisNode->mClusterId = true_cluster_id-1;
                thisNode->mVisited = true;
                true_cluster_comps_rev[thisNode->mId] = thisNode->mClusterId;
                true_cluster_comps_id[true_cluster_id-1].push_back(thisNode->mId);

                if(mpi_rank == 0)
                    std::cout << "\t- Adding global cluster " << thisNode->mId << std::endl;

                for(auto & neighbour : thisNode->mNeighbours) {
                    if(!neighbour->mVisited) {
                        frontNodes.push_back(neighbour);
                    }
                }
            }
        }

        std::vector<double> true_cluster_shrink(true_cluster_id, 0.0);
        std::vector<double> true_cluster_shrink_total(true_cluster_id, 0.0);
        std::vector<double> true_cluster_shrink_void(true_cluster_id, 0.0);
        std::vector<double> true_cluster_volume(true_cluster_id, 0.0);
        std::vector<double> true_cluster_liquid_volume(true_cluster_id, 0.0);
        std::vector<int>    true_cluster_open(true_cluster_id, 0);
        std::vector<int>    true_cluster_intensification_phase(true_cluster_id, 0);
        std::vector<int>    true_cluster_nodes_count(true_cluster_id, 0);
        std::vector<int>    true_cluster_liquid_nodes_count(true_cluster_id, 0);

        std::vector<std::vector<Node<3>::WeakPointer>> true_cluster_local_nodes(true_cluster_id, std::vector<Node<3>::WeakPointer>(0));

        for(std::size_t i = 0; i < true_cluster_comps_id.size(); i++) {
            for(auto & subcluster_id : true_cluster_comps_id[i]) {
                true_cluster_shrink[i] += global_cluster_shrink[subcluster_id];
                true_cluster_shrink_total[i] += global_cluster_shrink_total[subcluster_id];
                true_cluster_shrink_void[i] += global_cluster_shrink_void[subcluster_id];
                true_cluster_volume[i] += global_cluster_volume[subcluster_id];
                true_cluster_liquid_volume[i] += global_cluster_liquid_volume[subcluster_id];
                true_cluster_open[i] += global_cluster_open[subcluster_id];
                true_cluster_intensification_phase[i] += global_cluster_intensification_phase[subcluster_id];
                true_cluster_nodes_count[i] += global_cluster_nodes_count[subcluster_id];
                true_cluster_liquid_nodes_count[i] += global_cluster_liquid_nodes_count[subcluster_id];
            }
        }

        std::cout << true_cluster_comps_rev.size() << " " << std::endl;
        std::cout << true_cluster_local_nodes.size() << " " << std::endl;

        // // Print what we have found
        // if(mpi_rank == 0) {
        //     for(std::size_t i = 0; i < true_cluster_comps_id.size(); i++) {
        //         std::cout << "Found a true cluster (" << i <<  ") with id's: ";
        //         for(auto & v : true_cluster_comps_id[i]) {
        //             std::cout << v << "\t";
        //         }
        //         std::cout << std::endl;
        //     }
        // } 

        for(auto & node : mrModelPart.NodesArray()) {
            std::size_t node_cluster_id = node->FastGetSolutionStepValue(CLUSTER_ID);
            std::size_t node_true_cluster_id = node_cluster_id + local_domain_index[mpi_rank];
            true_cluster_local_nodes[true_cluster_comps_rev[node_true_cluster_id]].push_back(node);
            if(node->FastGetSolutionStepValue(CLUSTER_ID) != 0) {
                node->FastGetSolutionStepValue(CLUSTER_ID) = true_cluster_comps_rev[node_true_cluster_id];
            }
        } 

        // Calculate the second part of the shrinkage

        for(std::size_t i = 0; i < true_cluster_comps_id.size(); i++) {
            // If cluster has liquid nodes
            if(true_cluster_liquid_nodes_count[i]) {
                
                // If the volume is less than a treshhold, error.
                if(true_cluster_liquid_volume[i] <= 1.0e-18) {
                    true_cluster_liquid_volume[i] = 1.0e-18;

                    std::cout<< "****************************************************"<<std::endl;
                    std::cout<< "****************************************************"<<std::endl;
                    std::cout<< "*******ERROR COMPUTING LIQUID CLUSTER VOLUME *******"<<std::endl;
                    std::cout<< "****************************************************"<<std::endl;
                    std::cout<< "****************************************************"<<std::endl;
                }

                true_cluster_shrink_void[i] = true_cluster_shrink[i]; // Charlie: this line is repeated?
                true_cluster_shrink[i] /= true_cluster_liquid_nodes_count[i];

                for(auto & cluster_node_weak : true_cluster_local_nodes[i]) {
                    auto cluster_node = cluster_node_weak.lock();
                    if(cluster_node->GetSolutionStepValue(SOLIDFRACTION) < SOLIDIFICATION_LIMIT) {
                        mNodeShrink[cluster_node->Id()] += cluster_shrink[i];
                    }

                    if((cluster_volume[i] - cluster_liquid_volume[i]) > 0.0) {
                        double cluster_percent = 1.0 * (cluster_shrink_void[i] / (cluster_volume[i]-cluster_liquid_volume[i]));
                        cluster_percent = std::max(0.0, cluster_percent);
                        cluster_percent = std::min(1.0, cluster_percent);
                        // cluster_node->GetValue(SHRINKAGE_VOID) = cluster_percent;
                    }
                }

            } else {
                // The cluster has no liquid nodes

                for(auto & cluster_node_weak : true_cluster_local_nodes[i]) {
                    auto cluster_node = cluster_node_weak.lock();

                    double cluster_percent = 1.0 * (cluster_shrink_void[i] / (cluster_volume[i]-cluster_liquid_volume[i]));
                    cluster_percent = std::max(0.0, cluster_percent);
                    cluster_percent = std::min(1.0, cluster_percent);
                    // cluster_node->GetValue(SHRINKAGE_VOID) = cluster_percent;
                }

                // Writing Macroporosity
                if(true_cluster_nodes_count[i] < 8) {
                    for(auto & cluster_node_weak : true_cluster_local_nodes[i]) {
                        auto cluster_node = cluster_node_weak.lock();
                        cluster_node->GetValue(MACRO_POROSITY) = cluster_shrink[i];
                    }
                } else {
                    int n_division = 4; // The division of solid fraction
                    int counter = 0;    // 
                    double treshold_increment = (max_solid_fraction - min_solid_fraction) / n_division;

                    for(std::size_t j = 0; j < n_division; j++) {
                        double solid_fraction_threshold = min_solid_fraction + j * treshold_increment;

                        for(auto & cluster_node_weak : true_cluster_local_nodes[i]) {
                            auto cluster_node = cluster_node_weak.lock();

                            if(cluster_node->GetSolutionStepValue(SOLIDFRACTION, 1) < solid_fraction_threshold) {
                                counter++;
                                cluster_node->GetValue(MACRO_POROSITY) = cluster_shrink[i];
                            }
                        }
                    }
                    if (counter > 6) break;
                }
            }
        }

    }

    void Execute() {
        
        std::cout << "Calculating local clusters..." << std::endl;
        CalculateLocalClusters(local_cluster_size, local_cluster_nodes);
        std::cout << "Local clusters calculated \t OK!" << std::endl;

        std::cout << "Calculating local shrinkage..." << std::endl;
        CalculateLocalClusterShrinkage(local_cluster_nodes);
        std::cout << "Local shrinkage calculated \t OK!" << std::endl;

        MPI_Allreduce(&mLocalClusterNum, &mTotalClusterNum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        std::vector<int> local_adj_matrix(mTotalClusterNum*mTotalClusterNum,0);
        std::vector<int> total_adj_matrix(mTotalClusterNum*mTotalClusterNum,0);

        std::vector<int> local_domain_sizes(mpi_size, 0);
        std::vector<int> local_domain_index(mpi_size, 0);

        std::vector<std::shared_ptr<GraphNode>> global_adj_graph;

        local_domain_sizes[mpi_rank] = mLocalClusterNum;
        MPI_Allgather(&mLocalClusterNum, 1, MPI_INT, local_domain_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

        for(std::size_t i = 0; i < mpi_size; i++) {
            local_domain_index[i] = local_domain_index[i-1] + local_domain_sizes[i-1];
            std::cout << "XXXX= " << local_domain_index[i] << std::endl;
        }

        std::cout << "Calculating global clusters..." << std::endl;
        CalculateGlobalClusters(local_domain_sizes, local_domain_index, local_adj_matrix, total_adj_matrix, global_adj_graph);
        std::cout << "Global clusters calculated \t OK!" << std::endl;

        std::cout << "Calculating global shrinkage..." << std::endl;
        CalculateGlobalClusterShrinkage(local_domain_sizes, local_domain_index, total_adj_matrix, global_adj_graph);
        std::cout << "Global shrinkage clauclated \t OK!" << std::endl;
    }

private:

    int mpi_size;
    int mpi_rank;

    ModelPart & mrModelPart;

    int mTotalClusterNum;
    int mLocalClusterNum;

    double mSolidificationShrinkageFactor;
    bool mHighPressureFlag;

    double min_solid_fraction;
    double max_solid_fraction;

    std::unordered_map<std::size_t, double> mNodeShrink;

    std::vector<int>                                local_cluster_size;
    std::vector<std::vector<Node<3>::WeakPointer>>  local_cluster_nodes;

    std::vector<double> cluster_shrink;
    std::vector<double> cluster_shrink_total;
    std::vector<double> cluster_shrink_void;
    std::vector<double> cluster_volume;
    std::vector<double> cluster_liquid_volume;
    std::vector<double> cluster_min_solid_fraction;
    std::vector<double> cluster_max_solid_fraction;
    std::vector<int>    cluster_open;
    std::vector<int>    cluster_intensification_phase;
    std::vector<int>    cluster_nodes_count;
    std::vector<int>    cluster_liquid_nodes_count;

    std::vector<double> mLocalClusterShrinkage;

};


}