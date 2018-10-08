#include "includes/define.h"
#include "includes/element.h"
#include "containers/variable.h"

#include "statistics_record.h"
#include "statistics_data.h"

#include "fluid_dynamics_application_variables.h"

namespace Kratos {

void StatisticsRecord::AddResult(StatisticsSampler::Pointer pResult)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mInitialized) << "Trying to add statistical data after Initialization of the internal storage." << std::endl;

    std::size_t result_size = pResult->GetSize();
    pResult->SetOffset(mDataBufferSize);

    mDataBufferSize += result_size;
    mAverageData.push_back(pResult);

    KRATOS_CATCH("")
}

void StatisticsRecord::InitializeStorage(ModelPart::ElementsContainerType& rElements)
{
    mUpdateBuffer.resize(mDataBufferSize);

    // Note: this should be done on a serial loop to avoid race conditions.
    for (auto it_element = rElements.begin(); it_element != rElements.end(); ++it_element)
    {
        it_element->GetValue(TURBULENCE_STATISTICS_DATA).InitializeStorage(*it_element,mDataBufferSize);
    }
    mInitialized = true;
}

void StatisticsRecord::SampleIntegrationPointResults(ModelPart& rModelPart)
{
    mRecordedSteps++;

    ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
    std::vector<double> dummy;
    for( auto it_elem = rModelPart.ElementsBegin(); it_elem != rModelPart.ElementsEnd(); ++it_elem)
    {
        it_elem->GetValueOnIntegrationPoints(UPDATE_STATISTICS,dummy,r_process_info);
    }
}

void StatisticsRecord::UpdateStatistics(Element* pElement)
{
    KRATOS_DEBUG_ERROR_IF(!pElement->Has(TURBULENCE_STATISTICS_DATA))
    << "Trying to compute turbulent statistics, but " << pElement->Info()
    << " does not have TURBULENCE_STATISTICS_DATA defined." << std::endl;

    auto &r_elemental_statistics = pElement->GetValue(TURBULENCE_STATISTICS_DATA);
    r_elemental_statistics.UpdateMeasurement(pElement, mAverageData, mUpdateBuffer, mRecordedSteps);
}

std::vector<double> StatisticsRecord::OutputForTest(ModelPart::ElementsContainerType& rElements) const
{
    std::vector<double> result;
    for (auto it_element = rElements.begin(); it_element != rElements.end(); ++it_element )
    {
        auto& r_statistics = it_element->GetValue(TURBULENCE_STATISTICS_DATA);
        for (std::size_t g = 0; g < r_statistics.NumberOfIntegrationPoints(); g++)
        {
            auto data_iterator = r_statistics.DataIterator(g);
            for (auto it = data_iterator.begin(); it != data_iterator.end(); ++it)
            {
                result.push_back(*it / mRecordedSteps);
            }
        }
    }
    return result;
}

void StatisticsRecord::PrintToFile(const ModelPart& rModelPart) const
{
    // Open output file
    std::stringstream file_name;
    file_name << "gp_statistics_" << rModelPart.GetCommunicator().MyPID() << ".csv";
    std::ofstream stats_file;
    stats_file.open(file_name.str().c_str(), std::ios::out | std::ios::trunc);

    for (ModelPart::ElementsContainerType::const_iterator it = rModelPart.GetCommunicator().LocalMesh().ElementsBegin();
         it != rModelPart.GetCommunicator().LocalMesh().ElementsEnd(); it++)
    {
        auto &r_elemental_statistics = it->GetValue(TURBULENCE_STATISTICS_DATA);
        r_elemental_statistics.WriteToCSVOutput(stats_file, *it, mAverageData, mRecordedSteps);
    }

    stats_file.close();
}

KRATOS_CREATE_VARIABLE( StatisticsRecord::Pointer, STATISTICS_CONTAINER)

//TODO move somewhere else
KRATOS_CREATE_VARIABLE( StatisticsData, TURBULENCE_STATISTICS_DATA)

std::vector<double> StatisticsRecord::mUpdateBuffer = std::vector<double>();

}