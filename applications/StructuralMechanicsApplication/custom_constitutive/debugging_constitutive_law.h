// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Fernando Rastellini
//

#if !defined (KRATOS_DEBUGGING_LAW_H_INCLUDED)
#define  KRATOS_DEBUGGING_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// Size type definition
    typedef std::size_t             SizeType;

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
 * @class DebuggingConstitutiveLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This is a dummy constitutive law for debugging pourposes
 * @details This law stores the pointer to the real constitutive law. It stores all the necessary values into the contiguration parameters. This class overrides all the methods in the base class and calls the real constitutive law
 * @author Vicente Mataix Ferrandiz
 * @author Fernando Rastellini
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DebuggingConstitutiveLaw
    : public ConstitutiveLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Process info definition
    typedef ProcessInfo      ProcessInfoType;

    /// Base class definition
    typedef ConstitutiveLaw         BaseType;

    /// Index type definition
    typedef std::size_t            IndexType;

    /// Pointer definition of DebuggingConstitutiveLaw
    KRATOS_CLASS_POINTER_DEFINITION( DebuggingConstitutiveLaw );

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    DebuggingConstitutiveLaw();

    /**
     * @brief Default constructor. (With parameters)
     * @param NewParameters The configuration parameters of the new constitutive law
     */
    DebuggingConstitutiveLaw(Kratos::Parameters NewParameters);

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * @brief It reates a new constitutive law pointer
     * @param NewParameters The configuration parameters of the new constitutive law
     * @return a Pointer to the new constitutive law
     */
    ConstitutiveLaw::Pointer Create(Kratos::Parameters NewParameters) const override;

    /**
     * @brief Copy constructor.
     */
    DebuggingConstitutiveLaw (const DebuggingConstitutiveLaw& rOther);


    /**
     * @brief Destructor.
     */
    ~DebuggingConstitutiveLaw() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @return The working space dimension of the current constitutive law
     */
    SizeType WorkingSpaceDimension() override;

    /**
     * @return The size of the strain vector of the current constitutive law
     */
    SizeType GetStrainSize() override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (boolean)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<bool>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (integer)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<int>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Matrix)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Matrix>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (array of 3 components)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     * @note Fixed size array of 3 doubles (e.g. for 2D stresses, plastic strains, ...)
     */
    bool Has(const Variable<array_1d<double, 3 > >& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (array of 6 components)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     * @note Fixed size array of 6 doubles (e.g. for stresses, plastic strains, ...)
     */
    bool Has(const Variable<array_1d<double, 6 > >& rThisVariable) override;

    /**
     * @brief Returns the value of a specified variable (boolean)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    bool& GetValue(
        const Variable<bool>& rThisVariable,
        bool& rValue
        ) override;

    /**
     * @briefReturns the value of a specified variable (integer)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    int& GetValue(
        const Variable<int>& rThisVariable,
        int& rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& GetValue(
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& GetValue(
        const Variable<Vector>& rThisVariable,
        Vector& rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (Matrix)
     * @param rThisVariable the variable to be returned
     * @return rValue output: the value of the specified variable
     */
    Matrix& GetValue(
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (array of 3 components)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    array_1d<double, 3 > & GetValue(
        const Variable<array_1d<double, 3 > >& rThisVariable,
        array_1d<double, 3 > & rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (array of 6 components)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return the value of the specified variable
     */
    array_1d<double, 6 > & GetValue(
        const Variable<array_1d<double, 6 > >& rThisVariable,
        array_1d<double, 6 > & rValue
        ) override;

    /**
     * @brief Sets the value of a specified variable (boolean)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<bool>& rThisVariable,
        const bool& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (integer)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<int>& rThisVariable,
        const int& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
     void SetValue(
        const Variable<double>& rThisVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<Vector >& rThisVariable,
        const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (Matrix)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<Matrix >& rThisVariable,
        const Matrix& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (array of 3 components)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<array_1d<double, 3 > >& rThisVariable,
        const array_1d<double, 3 > & rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (array of 6 components)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<array_1d<double, 6 > >& rThisVariable,
        const array_1d<double, 6 > & rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculates the value of a specified variable (bool)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    bool& CalculateValue(
        Parameters& rParameterValues,
        const Variable<bool>& rThisVariable,
        bool& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (int)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    int& CalculateValue(
        Parameters& rParameterValues,
        const Variable<int>& rThisVariable,
        int& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(
        Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (Vector)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Vector& CalculateValue(
        Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (Matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
     Matrix& CalculateValue(
        Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (array of 3 components)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
     array_1d<double, 3 > & CalculateValue(
        Parameters& rParameterValues,
        const Variable<array_1d<double, 3 > >& rVariable,
        array_1d<double, 3 > & rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (array of 6 components)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return The value of the specified variable
     */
    array_1d<double, 6 > & CalculateValue(
        Parameters& rParameterValues,
        const Variable<array_1d<double, 6 > >& rVariable,
        array_1d<double, 6 > & rValue
        ) override;

     /**
      * @brief Is called to check whether the provided material parameters in the Properties match the requirements of current constitutive model.
      * @param rMaterialProperties the current Properties to be validated against.
      * @return true, if parameters are correct; false, if parameters are insufficient / faulty
      * @note  this has to be implemented by each constitutive model. Returns false in base class since no valid implementation is contained here.
      */
    bool ValidateInput(const Properties& rMaterialProperties) override;

    /**
     * @brief Returns the expected strain measure of this constitutive law (by default linear strains)
     * @return the expected strain measure
     */
    StrainMeasure GetStrainMeasure() override;

    /**
     * @brief Returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
     * @return the expected stress measure
     */
    StressMeasure GetStressMeasure() override;

    /**
     * @brief Returns whether this constitutive model is formulated in incremental strains/stresses
     * @note By default, all constitutive models should be formulated in total strains
     * @return true, if formulated in incremental strains/stresses, false otherwise
     */
    bool IsIncremental() override;

    /**
     * @brief This is to be called at the very beginning of the calculation
     * @details (e.g. from InitializeElement) in order to initialize all relevant attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues
        ) override;

    /**
     * @brief To be called at the beginning of each solution step
     * @details (e.g. from Element::InitializeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void InitializeSolutionStep(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief To be called at the end of each solution step
     * @details (e.g. from Element::FinalizeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void FinalizeSolutionStep(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief  To be called at the beginning of each step iteration
     * @details (e.g. from Element::InitializeNonLinearIteration)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void InitializeNonLinearIteration(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief To be called at the end of each step iteration
     * @details (e.g. from Element::FinalizeNonLinearIteration)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void FinalizeNonLinearIteration(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1 (Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2 (Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy (Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK1 (Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK2 (Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponseKirchhoff (Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void InitializeMaterialResponseCauchy (Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1 (Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2 (Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff (Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy (Parameters& rValues) override;

    /**
     * @brief This can be used in order to reset all internal variables of the
     * constitutive law (e.g. if a model should be reset to its reference state)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void ResetMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues
        ) override;

//     /**
//      * @brief Methods to transform strain Vectors:
//      * @param rStrainVector the strain tensor in matrix which its stress measure will be changed
//      * @param rF the DeformationGradientF matrix between the configurations
//      * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
//      * @param rStrainInitial the measure of stress of the given  rStrainVector
//      * @param rStrainFinal the measure of stress of the returned rStrainVector
//      */
//     Vector& TransformStrains(
//         Vector& rStrainVector,
//         const Matrix &rF,
//         StrainMeasure rStrainInitial,
//         StrainMeasure rStrainFinal
//         );

    /**
     * @brief Methods to transform stress Matrices:
     * @param rStressMatrix the stress tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStressInitial the measure of stress of the given  rStressMatrix
     * @param rStressFinal the measure of stress of the returned rStressMatrix
     */
    Matrix& TransformStresses(
        Matrix& rStressMatrix,
        const Matrix &rF,
        const double &rdetF,
        StressMeasure rStressInitial,
        StressMeasure rStressFinal
        ) override;

    /**
     * @brief Methods to transform stress Vectors:
     * @param rStressVector the stress tensor in matrix which its stress measure will be changed
     * @param rF the DeformationGradientF matrix between the configurations
     * @param rdetF the determinant of the DeformationGradientF matrix between the configurations
     * @param rStressInitial the measure of stress of the given  rStressVector
     * @param rStressFinal the measure of stress of the returned rStressVector
     */
    Vector& TransformStresses (
        Vector& rStressVector,
        const Matrix &rF,
        const double &rdetF,
        StressMeasure rStressInitial,
        StressMeasure rStressFinal
        ) override;

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return 0 if OK, 1 otherwise
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    //*** OUTDATED METHODS: ***//

    /**
     * @brief Computes the material response in terms of stresses and algorithmic tangent
     * @param StrainVector the current strains (total strains, input)
     * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear strain measure is used)
     * @param StressVector the computed stresses (output)
     * @param algorithmicTangent the material tangent matrix (output)
     * @param rCurrentProcessInfo current ProcessInfo instance
     * @param rMaterialProperties the material's Properties object
     * @param rElementGeometry the element's geometry
     * @param rShapeFunctionsValues the shape functions values in the current integration pointer
     * @param CalculateStresses flag whether or not to compute the stress response
     * @param CalculateTangent flag to determine if to compute the material tangent
     * @note The CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
     * @param SaveInternalVariables flag whether or not to store internal (history) variables
     */
    void CalculateMaterialResponse(
        const Vector& StrainVector,
        const Matrix& DeformationGradient,
        Vector& StressVector,
        Matrix& AlgorithmicTangent,
        const ProcessInfo& rCurrentProcessInfo,
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        bool CalculateStresses = true,
        int CalculateTangent = true,
        bool SaveInternalVariables = true
        ) override;

    /**
     * @brief Computes the volumetric part of the material response in terms of stresses and algorithmic tangent
     * @param StrainVector the current strains (total strains, input)
     * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear
     * @param StressVector the computed stresses (output)
     * @param algorithmicTangent the material tangent matrix (output)
     * @param rCurrentProcessInfo current ProcessInfo instance
     * @param rMaterialProperties the material's Properties object
     * @param rElementGeometry the element's geometry
     * @param rShapeFunctionsValues the shape functions values in the current integration pointer
     * @param CalculateStresses flag whether or not to compute the stress response
     * @param CalculateTangent flag to determine if to compute the material tangent
     * @note The CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
     * @param SaveInternalVariables flag whether or not to store internal (history) variables
     */
    void CalculateVolumetricResponse(
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
        ) override;

    /**
     * @brief Computes the deviatoric part of the material response in terms of stresses and algorithmic tangent
     * @param StrainVector the current strains (total strains, input)
     * @param DeformationGradient the current deformation gradient (can be an empty matrix if a linear
     * @param StressVector the computed stresses (output)
     * @param algorithmicTangent the material tangent matrix (output)
     * @param rCurrentProcessInfo current ProcessInfo instance
     * @param rMaterialProperties the material's Properties object
     * @param rElementGeometry the element's geometry
     * @param rShapeFunctionsValues the shape functions values in the current integration pointer
     * @param CalculateStresses flag whether or not to compute the stress response
     * @param CalculateTangent flag to determine if to compute the material tangent
     * @note the CalculateTangent flag is defined as int to allow for distinctive variants of the tangent
     * @param SaveInternalVariables flag whether or not to store internal (history) variables
     * @todo Add proper definition for algorithmic tangent
     */
    void CalculateDeviatoricResponse(
        const Vector& StrainVector,
        const Matrix& DeformationGradient,
        Vector& StressVector,
        Matrix& AlgorithmicTangent,
        const ProcessInfo& rCurrentProcessInfo,
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        bool CalculateStresses = true,
        int CalculateTangent = true,
        bool SaveInternalVariables = true
        ) override;

    // VM
    void CalculateCauchyStresses(
        Vector& Cauchy_StressVector,
        const Matrix& F,
        const Vector& PK2_StressVector,
        const Vector& GreenLagrangeStrainVector
        ) override;

    ///@}
    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "DebuggingConstitutiveLaw debugs the following law" << mpRealConstitutiveLaw->Info();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DebuggingConstitutiveLaw debugs the following law" << mpRealConstitutiveLaw->Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "DebuggingConstitutiveLaw debugs the following law" << mpRealConstitutiveLaw->Info();
    }

    ///@}
    ///@name Friends
    ///@{
    ///@}
protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ConstitutiveLaw::Pointer mpRealConstitutiveLaw = nullptr; /// The pointer to the real constitutive law
    Kratos::Parameters mThisParameters = Kratos::Parameters(R"({})"); /// The configuration parameters (will contain other values of interest that we may save for later debugging pourposes, such as element id and integration node id, etcetera)

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    ///@}

    std::string GetCatchMessage()
    {
        if (mThisParameters.Has("element_id") && mThisParameters.Has("gp_id"))
            return "Element id: " + std::to_string(mThisParameters["element_id"].GetInt()) + ". Integration GP id: " + std::to_string(mThisParameters["gp_id"].GetInt());
        else
            return "";

        return "";
    }

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw)
//         rSerializer.save("RealConstitutiveLaw", mpRealConstitutiveLaw);
//         rSerializer.save("ThisParameters", mThisParameters);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
//         rSerializer.load("RealConstitutiveLaw", mpRealConstitutiveLaw);
//         rSerializer.load("ThisParameters", mThisParameters);
    }


}; // Class DebuggingConstitutiveLaw
}  // namespace Kratos.
#endif // KRATOS_DEBUGGING_LAW_H_INCLUDED  defined
