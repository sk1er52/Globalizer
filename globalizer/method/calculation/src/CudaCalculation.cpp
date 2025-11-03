/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      CudaCalculation.cpp                                         //
//                                                                         //
//  Purpose:   Source file for CUDA calculation class                      //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "CudaCalculation.h"

// ------------------------------------------------------------------------------------------------

void CUDACalculation::StartCalculate(InformationForCalculation& inputSet,
    TResultForCalculation& outputSet)
{
    // Перераспределение памяти под буферы, если необходимо
    int newCoordinatesSize = inputSet.trials.size() * parameters.Dimension;
    if (newCoordinatesSize > coordinatesSize)
    {
        coordinatesSize = newCoordinatesSize;
        coordinates = new double[coordinatesSize];
    }

    int newFuncValuesSize = inputSet.trials.size() * pTask->GetNumOfFunc();
    if (newFuncValuesSize > FuncValuesSize)
    {
        FuncValuesSize = newFuncValuesSize;
        FuncValues = new double[FuncValuesSize];
    }

    // Копирование координат из массива испытаний в непрерывный буфер
    for (unsigned i = 0; i < inputSet.trials.size(); i++)
    {
        for (int k = 0; k < parameters.Dimension; k++)
            coordinates[i * parameters.Dimension + k] = inputSet.trials[i]->y[k];
    }

    // Вызов метода задачи, который должен содержать CUDA-вычисления
    pTask->CalculateFuncsInManyPoints(coordinates, 0, int(inputSet.trials.size()), FuncValues);

    // Копирование результатов из буфера обратно в объекты Trial
    for (unsigned i = 0; i < inputSet.trials.size(); i++)
    {
        outputSet.trials[i] = inputSet.trials[i];
        for (int k = 0; k < pTask->GetNumOfFunc(); k++)
            outputSet.trials[i]->FuncValues[k] = FuncValues[i + k * pTask->GetNumOfFunc()];
        outputSet.trials[i]->index = 0;
    }

    // Подсчет количества вычислений
    for (unsigned int i = 0; i < inputSet.trials.size(); i++)
    {
        for (int j = 0; j <= outputSet.trials[i]->index; j++)
            outputSet.countCalcTrials[j]++;
    }
}

// ------------------------------------------------------------------------------------------------
void CUDACalculation::Calculate(InformationForCalculation& inputSet,
    TResultForCalculation& outputSet)
{
    if (inputSet.trials.size() > 0)
    {
        outputSet.trials.clear();
        outputSet.trials.resize(inputSet.trials.size());
        for (unsigned i = 0; i < outputSet.trials.size(); i++)
            outputSet.trials[i] = inputSet.trials[i];


        outputSet.countCalcTrials.clear();
        outputSet.countCalcTrials.resize(pTask->GetNumOfFunc());
        for (int i = 0; i < pTask->GetNumOfFunc(); i++)
            outputSet.countCalcTrials[i] = 0;

    }

    // Запускать вычисления как только пришли данные
    if (isStartComputingAway)
    {
        StartCalculate(inputSet, outputSet);
    }
    else // Собрать данные в один блок, и потом вычислить все сразу
    {
        if (countCalculation > 0)
        {
            countCalculation--;

            for (unsigned i = 0; i < outputSet.trials.size(); i++)
            {
                inputCalculation.trials.push_back(inputSet.trials[i]);
            }

            if (countCalculation > 0)
                firstCalculation->ContinueComputing();

            for (unsigned int i = 0; i < inputSet.trials.size(); i++)
            {
                for (int j = 0; j <= outputSet.trials[i]->index; j++)
                    outputSet.countCalcTrials[j]++;
            }
        }

        if (countCalculation == 0)
        {
            countCalculation--;
            isStartComputingAway = true;


            for (unsigned i = 0; i < resultCalculation.trials.size(); i++)
                resultCalculation.trials[i] = inputCalculation.trials[i];

            StartCalculate(inputCalculation, resultCalculation);
        }
    }
}