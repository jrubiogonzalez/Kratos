# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Import Kratos core and apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Additional imports
from analyzer_base import AnalyzerBaseClass

# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters = Parameters(parameter_file.read())

# Definition of external analyzer
class CustomAnalyzer(AnalyzerBaseClass):
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):

        # Constraint 1
        constraint_node_id = 893
        if communicator.isRequestingValueOf("point_distance_893"):
            value = current_design.Nodes[constraint_node_id].Y
            communicator.reportValue("point_distance_893", value)

        if communicator.isRequestingGradientOf("point_distance_893"):
            gradient = {}
            for node in current_design.Nodes:
                if node.Id == constraint_node_id:
                    gradient[node.Id] = [0.0,1.0,0.0]
                else:
                    gradient[node.Id] = [0.0,0.0,0.0]
            communicator.reportGradient("point_distance_893", gradient)

        # Constraint 2
        constraint_node_id = 1200
        if communicator.isRequestingValueOf("point_distance_1200"):
            value = current_design.Nodes[constraint_node_id].Y
            communicator.reportValue("point_distance_1200", value)

        if communicator.isRequestingGradientOf("point_distance_1200"):
            gradient = {}
            for node in current_design.Nodes:
                if node.Id == constraint_node_id:
                    gradient[node.Id] = [0.0,1.0,0.0]
                else:
                    gradient[node.Id] = [0.0,0.0,0.0]
            communicator.reportGradient("point_distance_1200", gradient)            


# Defining the model_part
optimization_model_part = ModelPart(parameters["optimization_settings"]["design_variables"]["optimization_model_part_name"].GetString())
optimization_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, parameters["optimization_settings"]["design_variables"]["domain_size"].GetInt())

# Create optimizer and perform optimization
import optimizer_factory
optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], optimization_model_part, CustomAnalyzer())
optimizer.Optimize()