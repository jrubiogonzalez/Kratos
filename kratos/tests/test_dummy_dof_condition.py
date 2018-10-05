from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest
import os
import math

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestDummyDofCondition(KratosUnittest.TestCase):

    def test_dummy_dof_condition(self):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        model = KratosMultiphysics.Model()
        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_VAUX)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TANGENT_XI)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TANGENT_ETA)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_dummy_dof_condition"))
        model_part_io.ReadModelPart(model_part)
        model.AddModelPart(model_part)

        # Set CL
        cl = StructuralMechanicsApplication.LinearElastic3DLaw()
        model_part.GetProperties()[1].SetValue(KratosMultiphysics.CONSTITUTIVE_LAW,cl)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.NODAL_VAUX_X, KratosMultiphysics.FORCE_RESIDUAL_X,model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.NODAL_VAUX_Y, KratosMultiphysics.FORCE_RESIDUAL_Y,model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.NODAL_VAUX_Z, KratosMultiphysics.FORCE_RESIDUAL_Z,model_part)

        # Set vectors
        normal = KratosMultiphysics.Vector(3)
        normal[0] = 0.0
        normal[1] = 0.0
        normal[2] = 1.0
        tangent_xi = KratosMultiphysics.Vector(3)
        tangent_xi[0] = 1.0
        tangent_xi[1] = 0.0
        tangent_xi[2] = 0.0
        tangent_eta = KratosMultiphysics.Vector(3)
        tangent_eta[0] = 0.0
        tangent_eta[1] = 1.0
        tangent_eta[2] = 0.0
        for node in model_part.GetSubModelPart("DISPLACEMENT_Displacement_Auto2").Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.NORMAL, normal)
            node.SetSolutionStepValue(KratosMultiphysics.TANGENT_XI, tangent_xi)
            node.SetSolutionStepValue(KratosMultiphysics.TANGENT_ETA, tangent_eta)

        # Set delta disp
        delta_disp = KratosMultiphysics.Vector(3)
        delta_disp[0] = 0.0
        delta_disp[1] = 0.0
        delta_disp[2] = 1.0e-5
        for node in model_part.GetSubModelPart("DISPLACEMENT_Displacement_Auto1").Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, delta_disp)
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)

        # Add constraints
        KratosMultiphysics.DummyConditionsCreationUtility.CreateConstraints(model_part.GetSubModelPart("DISPLACEMENT_Displacement_Auto2"))

        linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolverWithConstraints(linear_solver)
        scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-4,1e-9)
        max_iters = 20
        compute_reactions = True
        reform_step_dofs = True
        calculate_norm_dx = False
        move_mesh_flag = True
        strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(model_part,
                                                                        scheme,
                                                                        linear_solver,
                                                                        convergence_criterion,
                                                                        builder_and_solver,
                                                                        max_iters,
                                                                        compute_reactions,
                                                                        reform_step_dofs,
                                                                        move_mesh_flag)
        convergence_criterion.SetEchoLevel(1)
        strategy.SetEchoLevel(1)

        strategy.Check()
        strategy.Solve()

        # DEBUG
        self._post_process(model_part)

        for node in model_part.GetSubModelPart("DISPLACEMENT_Displacement_Auto2").Nodes:
            disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
            self.assertGreater(disp[0], 0.0)
            self.assertGreater(disp[1], 0.0)

    def _post_process(self, model_part):
        from gid_output_process import GiDOutputProcess
        gid_output = GiDOutputProcess(model_part,
                                    "gid_output",
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditions",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                "nodal_results" : ["DISPLACEMENT","REACTION","NODAL_VAUX","NORMAL","TANGENT_XI", "TANGENT_ETA"]
                                            }
                                        }
                                        """)
                                    )

        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()


if __name__ == '__main__':
    KratosUnittest.main()

