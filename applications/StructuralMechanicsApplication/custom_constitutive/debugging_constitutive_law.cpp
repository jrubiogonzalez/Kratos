// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/properties.h"
#include "custom_constitutive/debugging_constitutive_law.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

DebuggingConstitutiveLaw::DebuggingConstitutiveLaw()
    : ConstitutiveLaw()
{
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

DebuggingConstitutiveLaw::DebuggingConstitutiveLaw(Kratos::Parameters NewParameters)
    : ConstitutiveLaw()
{
    Kratos::Parameters default_parameters = Kratos::Parameters(R"(
    {
        "name"                      : "DebuggingConstitutiveLaw",
        "constitutive_law_to_debug" : ""
    })" );

    mThisParameters.ValidateAndAssignDefaults(default_parameters);

    const std::string& constitutive_law_to_debug = NewParameters["constitutive_law_to_debug"].GetString();
    mpRealConstitutiveLaw = KratosComponents<ConstitutiveLaw>::Get(constitutive_law_to_debug).Clone();
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

DebuggingConstitutiveLaw::DebuggingConstitutiveLaw(const DebuggingConstitutiveLaw& rOther)
    : ConstitutiveLaw(rOther),
      mpRealConstitutiveLaw(rOther.mpRealConstitutiveLaw),
      mThisParameters(rOther.mThisParameters)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer DebuggingConstitutiveLaw::Clone() const
{
    return Kratos::make_shared<DebuggingConstitutiveLaw>(*this);
}

//********************************CREATE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer DebuggingConstitutiveLaw::Create(Kratos::Parameters NewParameters) const
{
    return Kratos::make_shared<DebuggingConstitutiveLaw>(NewParameters);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

DebuggingConstitutiveLaw::~DebuggingConstitutiveLaw()
{
}

/***********************************************************************************/
/***********************************************************************************/

SizeType DebuggingConstitutiveLaw::WorkingSpaceDimension()
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->WorkingSpaceDimension();

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

SizeType DebuggingConstitutiveLaw::GetStrainSize()
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->GetStrainSize();

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

bool DebuggingConstitutiveLaw::Has(const Variable<bool>& rThisVariable)
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->Has(rThisVariable);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

bool DebuggingConstitutiveLaw::Has(const Variable<int>& rThisVariable)
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->Has(rThisVariable);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

bool DebuggingConstitutiveLaw::Has(const Variable<double>& rThisVariable)
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->Has(rThisVariable);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

bool DebuggingConstitutiveLaw::Has(const Variable<Vector>& rThisVariable)
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->Has(rThisVariable);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

bool DebuggingConstitutiveLaw::Has(const Variable<Matrix>& rThisVariable)
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->Has(rThisVariable);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

bool DebuggingConstitutiveLaw::Has(const Variable<array_1d<double, 3 > >& rThisVariable)
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->Has(rThisVariable);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

bool DebuggingConstitutiveLaw::Has(const Variable<array_1d<double, 6 > >& rThisVariable)
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->Has(rThisVariable);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

bool& DebuggingConstitutiveLaw::GetValue(
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->GetValue(rThisVariable, rValue);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

int& DebuggingConstitutiveLaw::GetValue(
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->GetValue(rThisVariable, rValue);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

double& DebuggingConstitutiveLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->GetValue(rThisVariable, rValue);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

Vector& DebuggingConstitutiveLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->GetValue(rThisVariable, rValue);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& DebuggingConstitutiveLaw::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->GetValue(rThisVariable, rValue);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3 >& DebuggingConstitutiveLaw::GetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->GetValue(rThisVariable, rValue);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 6 >& DebuggingConstitutiveLaw::GetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->GetValue(rThisVariable, rValue);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::SetValue(
    const Variable<bool>& rThisVariable,
    const bool& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::SetValue(
    const Variable<int>& rThisVariable,
    const int& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::SetValue(
    const Variable<Matrix>& rThisVariable,
    const Matrix& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::SetValue(
    const Variable<array_1d<double, 3 >>& rThisVariable,
    const array_1d<double, 3 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::SetValue(
    const Variable<array_1d<double, 6 >>& rThisVariable,
    const array_1d<double, 6 >& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

bool& DebuggingConstitutiveLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, rValue);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

int& DebuggingConstitutiveLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, rValue);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

double& DebuggingConstitutiveLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, rValue);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

Vector& DebuggingConstitutiveLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, rValue);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& DebuggingConstitutiveLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, rValue);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3 >& DebuggingConstitutiveLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<array_1d<double, 3 >>& rThisVariable,
    array_1d<double, 3 >& rValue
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, rValue);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 6 >& DebuggingConstitutiveLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<array_1d<double, 6 >>& rThisVariable,
    array_1d<double, 6 >& rValue
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, rValue);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

bool DebuggingConstitutiveLaw::ValidateInput(const Properties& rMaterialProperties)
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->ValidateInput(rMaterialProperties);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

ConstitutiveLaw::StrainMeasure DebuggingConstitutiveLaw::GetStrainMeasure()
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->GetStrainMeasure();

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

ConstitutiveLaw::StressMeasure DebuggingConstitutiveLaw::GetStressMeasure()
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->GetStressMeasure();

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

bool DebuggingConstitutiveLaw::IsIncremental()
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->IsIncremental();

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::InitializeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->InitializeSolutionStep(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->FinalizeSolutionStep(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::InitializeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->InitializeNonLinearIteration(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::FinalizeNonLinearIteration(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->FinalizeNonLinearIteration(rMaterialProperties, rElementGeometry, rShapeFunctionsValues, rCurrentProcessInfo);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::CalculateMaterialResponsePK1(Parameters& rValues)
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->CalculateMaterialResponsePK1(rValues);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->CalculateMaterialResponsePK2(rValues);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->CalculateMaterialResponseKirchhoff(rValues);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->CalculateMaterialResponseCauchy(rValues);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::InitializeMaterialResponsePK1(Parameters& rValues)
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->InitializeMaterialResponsePK1(rValues);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::InitializeMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->InitializeMaterialResponsePK2(rValues);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::InitializeMaterialResponseKirchhoff(Parameters& rValues)
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->InitializeMaterialResponseKirchhoff(rValues);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::InitializeMaterialResponseCauchy(Parameters& rValues)
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->InitializeMaterialResponseCauchy(rValues);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->FinalizeMaterialResponsePK1(rValues);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->FinalizeMaterialResponsePK2(rValues);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->FinalizeMaterialResponseKirchhoff(rValues);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->FinalizeMaterialResponseCauchy(rValues);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::ResetMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->ResetMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

    KRATOS_CATCH(GetCatchMessage());
}

// /***********************************************************************************/
// /***********************************************************************************/
//
// Vector& DebuggingConstitutiveLaw::TransformStrains(
//     Vector& rStrainVector,
//     const Matrix &rF,
//     StrainMeasure rStrainInitial,
//     StrainMeasure rStrainFinal
//     )
// {
//     KRATOS_TRY;
//
//     // We call to the real CL
//     return mpRealConstitutiveLaw->TransformStrains(rStrainVector, rF, rStressInitial, rStressFinal);
//
//     KRATOS_CATCH(GetCatchMessage());
// }

/***********************************************************************************/
/***********************************************************************************/

Matrix& DebuggingConstitutiveLaw::TransformStresses(
    Matrix& rStressMatrix,
    const Matrix &rF,
    const double &rdetF,
    ConstitutiveLaw::StressMeasure rStressInitial,
    ConstitutiveLaw::StressMeasure rStressFinal
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->TransformStresses(rStressMatrix, rF, rdetF, rStressInitial, rStressFinal);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

Vector& DebuggingConstitutiveLaw::TransformStresses(
    Vector& rStressVector,
    const Matrix &rF,
    const double &rdetF,
    ConstitutiveLaw::StressMeasure rStressInitial,
    ConstitutiveLaw::StressMeasure rStressFinal
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->TransformStresses(rStressVector, rF, rdetF, rStressInitial, rStressFinal);

    KRATOS_CATCH(GetCatchMessage());
}

/**************************CONSTITUTIVE LAW GENERAL FEATURES ***********************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->GetLawFeatures(rFeatures);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

int DebuggingConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // We call to the real CL
    return mpRealConstitutiveLaw->Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::CalculateMaterialResponse(
    const Vector& StrainVector,
    const Matrix& DeformationGradient,
    Vector& StressVector,
    Matrix& AlgorithmicTangent,
    const ProcessInfo& rCurrentProcessInfo,
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    bool CalculateStresses,
    int CalculateTangent,
    bool SaveInternalVariables
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->CalculateMaterialResponse(StrainVector, DeformationGradient, StressVector, AlgorithmicTangent, rCurrentProcessInfo, rMaterialProperties, rElementGeometry, rShapeFunctionsValues, CalculateStresses, CalculateTangent, SaveInternalVariables);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::CalculateVolumetricResponse(
    const double VolumetricStrain,
    const Matrix& DeformationGradient,
    double& VolumetricStress,
    double& AlgorithmicBulk,
    const ProcessInfo& rCurrentProcessInfo,
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    bool CalculateStresses,
    int CalculateTangent,
    bool SaveInternalVariables
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->CalculateVolumetricResponse(VolumetricStrain, DeformationGradient, VolumetricStress, AlgorithmicBulk, rCurrentProcessInfo, rMaterialProperties, rElementGeometry, rShapeFunctionsValues, CalculateStresses, CalculateTangent, SaveInternalVariables);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::CalculateDeviatoricResponse(
    const Vector& StrainVector,
    const Matrix& DeformationGradient,
    Vector& StressVector,
    Matrix& AlgorithmicTangent,
    const ProcessInfo& rCurrentProcessInfo,
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    bool CalculateStresses,
    int CalculateTangent,
    bool SaveInternalVariables
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->CalculateDeviatoricResponse(StrainVector, DeformationGradient, StressVector, AlgorithmicTangent, rCurrentProcessInfo, rMaterialProperties, rElementGeometry, rShapeFunctionsValues, CalculateStresses, CalculateTangent, SaveInternalVariables);

    KRATOS_CATCH(GetCatchMessage());
}

/***********************************************************************************/
/***********************************************************************************/

void DebuggingConstitutiveLaw::CalculateCauchyStresses(
    Vector& Cauchy_StressVector,
    const Matrix& F,
    const Vector& PK2_StressVector,
    const Vector& GreenLagrangeStrainVector
    )
{
    KRATOS_TRY;

    // We call to the real CL
    mpRealConstitutiveLaw->CalculateCauchyStresses(Cauchy_StressVector, F, PK2_StressVector, GreenLagrangeStrainVector);

    KRATOS_CATCH(GetCatchMessage());
}

} // Namespace Kratos
