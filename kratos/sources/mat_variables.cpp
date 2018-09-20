//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//			 Kratos default license: kratos/license.txt
//
//  Main authors:    Josep Maria Carbonell
//

// System includes
#include <string>
#include <iostream>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/mat_variables.h"
#include "includes/kernel.h"
#include "includes/node.h"

#include "includes/kratos_flags.h"

//commented variables are defined in the variables.h and variables.cpp

namespace Kratos
{
  //solution
  KRATOS_CREATE_VARIABLE( std::string, CONSTITUTIVE_LAW_NAME )
  KRATOS_CREATE_VARIABLE( bool, IMPLEX ) 
  KRATOS_CREATE_VARIABLE( bool, IMPLEX_CONTACT ) 

  //elasticity
  //KRATOS_CREATE_VARIABLE( double, YOUNG_MODULUS )
  //KRATOS_CREATE_VARIABLE( double, POISSON_RATIO )
  //KRATOS_CREATE_VARIABLE( double, DENSITY )
  //KRATOS_CREATE_VARIABLE( double, THICKNESS )
  //KRATOS_CREATE_VARIABLE( double, EQUIVALENT_YOUNG_MODULUS )
  //KRATOS_CREATE_VARIABLE( double, BULK_MODULUS )

  KRATOS_CREATE_VARIABLE( double, SHEAR_MODULUS )
  KRATOS_CREATE_VARIABLE( double, LAME_MU )
  KRATOS_CREATE_VARIABLE( double, LAME_LAMBDA )
  KRATOS_CREATE_VARIABLE( double, C10 )
  KRATOS_CREATE_VARIABLE( double, C20 )
  KRATOS_CREATE_VARIABLE( double, C30 )

  //viscosity
  //KRATOS_CREATE_VARIABLE( double, DYNAMIC_VISCOSITY )
  //KRATOS_CREATE_VARIABLE( double, VISCOSITY )
  
  //damage
  KRATOS_CREATE_VARIABLE( double, DAMAGE_VARIABLE )
  KRATOS_CREATE_VARIABLE( double, DAMAGE_THRESHOLD )
  KRATOS_CREATE_VARIABLE( double, STRENGTH_RATIO )
  KRATOS_CREATE_VARIABLE( double, FRACTURE_ENERGY )
  KRATOS_CREATE_VARIABLE( double, RESIDUAL_STRENGTH )

  //plasticity
  KRATOS_CREATE_VARIABLE( double, PLASTIC_STRAIN )
  KRATOS_CREATE_VARIABLE( double, DELTA_PLASTIC_STRAIN )
  KRATOS_CREATE_VARIABLE( double, NORM_ISOCHORIC_STRESS )
  KRATOS_CREATE_VARIABLE( double, PLASTIC_STRAIN_RATE )
  
  //hardening
  KRATOS_CREATE_VARIABLE( double, ISOTROPIC_HARDENING_MODULUS )
  KRATOS_CREATE_VARIABLE( double, KINEMATIC_HARDENING_MODULUS )
  KRATOS_CREATE_VARIABLE( double, HARDENING_EXPONENT )
  KRATOS_CREATE_VARIABLE( double, REFERENCE_HARDENING_MODULUS )
  KRATOS_CREATE_VARIABLE( double, INFINITY_HARDENING_MODULUS )
  KRATOS_CREATE_VARIABLE( double, SOFTENING_SLOPE )

  //baker-johnson-cook parameters
  KRATOS_CREATE_VARIABLE( double, JC_PARAMETER_A )
  KRATOS_CREATE_VARIABLE( double, JC_PARAMETER_B )
  KRATOS_CREATE_VARIABLE( double, JC_PARAMETER_C )
  KRATOS_CREATE_VARIABLE( double, JC_PARAMETER_m )
  KRATOS_CREATE_VARIABLE( double, JC_PARAMETER_n )
  KRATOS_CREATE_VARIABLE( double, JC_PARAMETER_K ) 
  
  //thermal
  //KRATOS_CREATE_VARIABLE( double, THERMAL_EXPANSION_COEFFICIENT )
  KRATOS_CREATE_VARIABLE( double, REFERENCE_CONDUCTIVITY )
  KRATOS_CREATE_VARIABLE( double, HARDNESS_CONDUCTIVITY )  
  KRATOS_CREATE_VARIABLE( double, REFERENCE_TEMPERATURE )
  KRATOS_CREATE_VARIABLE( double, MELD_TEMPERATURE )
  KRATOS_CREATE_VARIABLE( double, PLASTIC_DISSIPATION )
  KRATOS_CREATE_VARIABLE( double, DELTA_PLASTIC_DISSIPATION )
  
  //anisotropy
  KRATOS_CREATE_VARIABLE( double, YOUNG_MODULUS_X )
  KRATOS_CREATE_VARIABLE( double, YOUNG_MODULUS_Y )
  KRATOS_CREATE_VARIABLE( double, YOUNG_MODULUS_Z )
  KRATOS_CREATE_VARIABLE( double, SHEAR_MODULUS_XY )
  KRATOS_CREATE_VARIABLE( double, SHEAR_MODULUS_YZ )
  KRATOS_CREATE_VARIABLE( double, SHEAR_MODULUS_XZ )
  KRATOS_CREATE_VARIABLE( double, POISSON_RATIO_XY )
  KRATOS_CREATE_VARIABLE( double, POISSON_RATIO_YZ )
  KRATOS_CREATE_VARIABLE( double, POISSON_RATIO_XZ )
  
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MATERIAL_ORIENTATION_DX )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MATERIAL_ORIENTATION_DY )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MATERIAL_ORIENTATION_DZ )

  //critical state
  KRATOS_CREATE_VARIABLE( double, CRITICAL_STATE_LINE )
  KRATOS_CREATE_VARIABLE( double, PRE_CONSOLIDATION_STRESS )
  KRATOS_CREATE_VARIABLE( double, OVER_CONSOLIDATION_RATIO )
  KRATOS_CREATE_VARIABLE( double, INITIAL_SHEAR_MODULUS )
  KRATOS_CREATE_VARIABLE( double, NORMAL_COMPRESSION_SLOPE )
  KRATOS_CREATE_VARIABLE( double, SWELLING_SLOPE )
  KRATOS_CREATE_VARIABLE( double, ALPHA_SHEAR )

  //strain
  //KRATOS_CREATE_VARIABLE( double, DETERMINANT_F )
  //KRATOS_CREATE_VARIABLE( Vector, DEFORMATION_GRADIENT);
  KRATOS_CREATE_VARIABLE( Vector, INITIAL_STRAIN_VECTOR )
    
  KRATOS_CREATE_VARIABLE( Vector, GREEN_LAGRANGE_STRAIN_VECTOR )
  //KRATOS_CREATE_VARIABLE( Vector, GREEN_LAGRANGE_STRAIN_TENSOR )
    
  //KRATOS_CREATE_VARIABLE( Vector, HENCKY_STRAIN_VECTOR);
  //KRATOS_CREATE_VARIABLE( Matrix, HENCKY_STRAIN_TENSOR);
  
  KRATOS_CREATE_VARIABLE( Vector, ALMANSI_STRAIN_VECTOR )
  KRATOS_CREATE_VARIABLE( Matrix, ALMANSI_STRAIN_TENSOR )
  
  //stress
  KRATOS_CREATE_VARIABLE( Vector, KIRCHHOFF_STRESS_VECTOR )
  KRATOS_CREATE_VARIABLE( Matrix, KIRCHHOFF_STRESS_TENSOR )
  //KRATOS_CREATE_VARIABLE( Vector, CAUCHY_STRESS_VECTOR )
  //KRATOS_CREATE_VARIABLE( Matrix, CAUCHY_STRESS_TENSOR )    
  //KRATOS_CREATE_VARIABLE( Vector, PK2_STRESS_VECTOR )    
  //KRATOS_CREATE_VARIABLE( Matrix, PK2_STRESS_TENSOR )

    // Constitutive matrices
//   KRATOS_CREATE_VARIABLE( Matrix, LOCAL_CONSTITUTIVE_MATRIX )
//   KRATOS_CREATE_VARIABLE( Matrix, CONSTITUTIVE_MATRIX )
  KRATOS_CREATE_VARIABLE( Matrix, CONSTITUTIVE_MATRIX_PK2 )
  KRATOS_CREATE_VARIABLE( Matrix, CONSTITUTIVE_MATRIX_KIRCHHOFF )
  
  // LMV. Temporal place to include Biot's variables
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WATER_DISPLACEMENT )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WATER_VELOCITY )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WATER_ACCELERATION )
  
  void KratosApplication::RegisterMATVariables()
  {

    //solution
    KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_LAW_NAME )
    KRATOS_REGISTER_VARIABLE( IMPLEX ) 
    KRATOS_REGISTER_VARIABLE( IMPLEX_CONTACT ) 

    //elasticity
    //KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS )
    //KRATOS_REGISTER_VARIABLE( POISSON_RATIO )
    //KRATOS_REGISTER_VARIABLE( DENSITY )
    //KRATOS_REGISTER_VARIABLE( THICKNESS )
    //KRATOS_REGISTER_VARIABLE( EQUIVALENT_YOUNG_MODULUS )
    //KRATOS_REGISTER_VARIABLE( BULK_MODULUS )
    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS )
    KRATOS_REGISTER_VARIABLE( LAME_MU )
    KRATOS_REGISTER_VARIABLE( LAME_LAMBDA )
    KRATOS_REGISTER_VARIABLE( C10 )
    KRATOS_REGISTER_VARIABLE( C20 )
    KRATOS_REGISTER_VARIABLE( C30 )

    //viscosity
    //KRATOS_REGISTER_VARIABLE( DYNAMIC_VISCOSITY )
    //KRATOS_REGISTER_VARIABLE( VISCOSITY )
  
    //damage
    KRATOS_REGISTER_VARIABLE( DAMAGE_VARIABLE )
    KRATOS_REGISTER_VARIABLE( DAMAGE_THRESHOLD )
    KRATOS_REGISTER_VARIABLE( STRENGTH_RATIO )
    KRATOS_REGISTER_VARIABLE( FRACTURE_ENERGY )
    KRATOS_REGISTER_VARIABLE( RESIDUAL_STRENGTH )

    //plasticity
    KRATOS_REGISTER_VARIABLE( PLASTIC_STRAIN )
    KRATOS_REGISTER_VARIABLE( DELTA_PLASTIC_STRAIN )
    KRATOS_REGISTER_VARIABLE( NORM_ISOCHORIC_STRESS )
    KRATOS_REGISTER_VARIABLE( PLASTIC_STRAIN_RATE )
      
    //hardening
    KRATOS_REGISTER_VARIABLE( ISOTROPIC_HARDENING_MODULUS )
    KRATOS_REGISTER_VARIABLE( KINEMATIC_HARDENING_MODULUS )
    KRATOS_REGISTER_VARIABLE( HARDENING_EXPONENT )
    KRATOS_REGISTER_VARIABLE( REFERENCE_HARDENING_MODULUS )
    KRATOS_REGISTER_VARIABLE( INFINITY_HARDENING_MODULUS )
    KRATOS_REGISTER_VARIABLE( SOFTENING_SLOPE )

    //baker-johnson-cook parameters
    KRATOS_REGISTER_VARIABLE( JC_PARAMETER_A )
    KRATOS_REGISTER_VARIABLE( JC_PARAMETER_B )
    KRATOS_REGISTER_VARIABLE( JC_PARAMETER_C )
    KRATOS_REGISTER_VARIABLE( JC_PARAMETER_m )
    KRATOS_REGISTER_VARIABLE( JC_PARAMETER_n )
    KRATOS_REGISTER_VARIABLE( JC_PARAMETER_K ) 
  
    //thermal
    //KRATOS_REGISTER_VARIABLE( THERMAL_EXPANSION_COEFFICIENT )
    KRATOS_REGISTER_VARIABLE( REFERENCE_CONDUCTIVITY )
    KRATOS_REGISTER_VARIABLE( HARDNESS_CONDUCTIVITY )  
    KRATOS_REGISTER_VARIABLE( REFERENCE_TEMPERATURE )
    KRATOS_REGISTER_VARIABLE( MELD_TEMPERATURE )
    KRATOS_REGISTER_VARIABLE( PLASTIC_DISSIPATION )
    KRATOS_REGISTER_VARIABLE( DELTA_PLASTIC_DISSIPATION )
  
    //anisotropy
    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_X )
    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_Y )
    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_Z )
    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_XY )
    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_YZ )
    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_XZ )
    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_XY )
    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_YZ )
    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_XZ )
  
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MATERIAL_ORIENTATION_DX )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MATERIAL_ORIENTATION_DY )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MATERIAL_ORIENTATION_DZ )

    //critical state
    KRATOS_REGISTER_VARIABLE( CRITICAL_STATE_LINE )
    KRATOS_REGISTER_VARIABLE( PRE_CONSOLIDATION_STRESS )
    KRATOS_REGISTER_VARIABLE( OVER_CONSOLIDATION_RATIO )
    KRATOS_REGISTER_VARIABLE( INITIAL_SHEAR_MODULUS )
    KRATOS_REGISTER_VARIABLE( NORMAL_COMPRESSION_SLOPE )
    KRATOS_REGISTER_VARIABLE( SWELLING_SLOPE )
    KRATOS_REGISTER_VARIABLE( ALPHA_SHEAR )

    //strain
    //KRATOS_REGISTER_VARIABLE( DETERMINANT_F )
    //KRATOS_REGISTER_VARIABLE( DEFORMATION_GRADIENT);
    KRATOS_REGISTER_VARIABLE( INITIAL_STRAIN_VECTOR )
    
    KRATOS_REGISTER_VARIABLE( GREEN_LAGRANGE_STRAIN_VECTOR )
    //KRATOS_REGISTER_VARIABLE( GREEN_LAGRANGE_STRAIN_TENSOR )
    
    //KRATOS_REGISTER_VARIABLE( HENCKY_STRAIN_VECTOR);
    //KRATOS_REGISTER_VARIABLE( HENCKY_STRAIN_TENSOR);
  
    KRATOS_REGISTER_VARIABLE( ALMANSI_STRAIN_VECTOR )
    KRATOS_REGISTER_VARIABLE( ALMANSI_STRAIN_TENSOR )

    //stress
    KRATOS_REGISTER_VARIABLE( KIRCHHOFF_STRESS_VECTOR )
    KRATOS_REGISTER_VARIABLE( KIRCHHOFF_STRESS_TENSOR )
    //KRATOS_REGISTER_VARIABLE( CAUCHY_STRESS_VECTOR )
    //KRATOS_REGISTER_VARIABLE( CAUCHY_STRESS_TENSOR )    
    //KRATOS_REGISTER_VARIABLE( PK2_STRESS_VECTOR )    
    //KRATOS_REGISTER_VARIABLE( PK2_STRESS_TENSOR )
    
    // Constitutive matrices
//   KRATOS_REGISTER_VARIABLE( LOCAL_CONSTITUTIVE_MATRIX )
//   KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_MATRIX )
  KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_MATRIX_PK2 )
  KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_MATRIX_KIRCHHOFF )

    // LMV. Temporal place to include Biot's variables
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( WATER_DISPLACEMENT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( WATER_VELOCITY )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( WATER_ACCELERATION )

  }


}  // namespace Kratos.
