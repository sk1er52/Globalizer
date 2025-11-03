/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2021 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      MPICalculationAsync.cpp                                     //
//                                                                         //
//  Purpose:   Source file for Async MPI calculation class                 //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "MPICalculationAsync.h"

#include <stdlib.h>
#include <string.h>
#include <cstring>

#include "TaskFactory.h"
#include "TrialFactory.h"

// ------------------------------------------------------------------------------------------------

void MPICalculationAsync::AsyncFinilize()
{
    // Необходимо собрать данные со всех процессов, кроме того, который
    // последним прислал ответ, так как он уже свободен.
    MPI_Status status;
    Trial OptimEstimation;
    int Child;

    if (parameters.DebugAsyncCalculation != 0) {
        std::ofstream fout;
        fout.open("../_build/async.txt");
        fout << parameters.GetProcNum() << "\n";
        fout << parameters.GetProcNum();
        fout.close();
    }

    for (int i = 0; i < parameters.GetProcNum() - 2; i++) {
        MPI_Recv(&OptimEstimation.index, 1, MPI_INT, MPI_ANY_SOURCE, TagChildSolved, MPI_COMM_WORLD, &status);
        Child = status.MPI_SOURCE;
        MPI_Recv(OptimEstimation.y, parameters.Dimension, MPI_DOUBLE, Child, TagChildSolved, MPI_COMM_WORLD, &status);
        MPI_Recv(OptimEstimation.FuncValues, MaxNumOfFunc, MPI_DOUBLE, Child, TagChildSolved, MPI_COMM_WORLD, &status);
    }
}

// ------------------------------------------------------------------------------------------------

void MPICalculationAsync::RecieveCalculatedFunctional()
{
    MPI_Status status;
    int index;
    int i;

    // Принимаем индекс точки
    MPI_Recv(&index, 1, MPI_INT, MPI_ANY_SOURCE, TagChildSolved, MPI_COMM_WORLD, &status);
    ChildNumRecv = status.MPI_SOURCE; // MPI-ранг ответившего процесса

    // Запоминаем внутренний индекс процесса
    ChildNum = ChildNumRecv - 1;
    vecTrials[ChildNum]->index = index;
    // Принимаем координаты точки
    MPI_Recv(vecTrials[ChildNum]->y, parameters.Dimension, MPI_DOUBLE, ChildNumRecv, TagChildSolved, MPI_COMM_WORLD, &status);
    // Принимаем значения функционалов
    MPI_Recv(vecTrials[ChildNum]->FuncValues, MaxNumOfFunc, MPI_DOUBLE, ChildNumRecv, TagChildSolved, MPI_COMM_WORLD, &status);

    // Определяем итоговый индекс
    int fNumber = 0;
    vecTrials[ChildNum]->index = -1;
    while ((vecTrials[ChildNum]->index == -1) && (fNumber < pTask->GetNumOfFunc()))
    {
        if ((fNumber == (pTask->GetNumOfFunc() - 1)) || (vecTrials[ChildNum]->FuncValues[fNumber] > 0))
        {
            vecTrials[ChildNum]->index = fNumber;
        }
        fNumber++;
    }
}

// ------------------------------------------------------------------------------------------------

void MPICalculationAsync::FirstStartCalculate(InformationForCalculation& inputSet,
    TResultForCalculation& outputSet)
{
    isFirst = false;

    vecTrials.reserve(parameters.GetProcNum() - 1);

    for (unsigned int i = 0; i < parameters.NumPoints; i++)
    {
        int isFinish = 0;
        // Отправляем флаг, что есть работа
        MPI_Send(&isFinish, 1, MPI_INT, i + 1, TagChildSolved, MPI_COMM_WORLD);

        vecTrials.push_back(inputSet.trials[i]);
        Trial* trail = inputSet.trials[i];
        trail->index = 0;
        vecTrials[i]->index = 0;
        // Отправляем координаты y
        MPI_Send(trail->y, parameters.Dimension, MPI_DOUBLE, i + 1, TagChildSolved, MPI_COMM_WORLD);
    }

    RecieveCalculatedFunctional();

    outputSet.trials[ChildNum] = vecTrials[ChildNum];
    for (int j = 0; j <= outputSet.trials[ChildNum]->index; j++)
        outputSet.countCalcTrials[j]++;
    vecTrials[ChildNum] = 0;
}


// ------------------------------------------------------------------------------------------------

void MPICalculationAsync::StartCalculate(InformationForCalculation& inputSet,
    TResultForCalculation& outputSet)
{
    int isFinish = 0;
    if (ChildNumRecv < 1 || ChildNumRecv >= parameters.GetProcNum()) {
        std::cout << "Error with CHILDNUMRECV!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    }
    if (ChildNum < 0 || ChildNum >= vecTrials.size()) {
        std::cout << "Error with CHILDNUM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    }
    if (vecTrials[ChildNum] != 0 || inputSet.trials[0] == 0) {
        std::cout << "Error with trial!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    }

    MPI_Send(&isFinish, 1, MPI_INT, ChildNumRecv, TagChildSolved, MPI_COMM_WORLD);

    vecTrials[ChildNum] = inputSet.trials[0];

    MPI_Send(vecTrials[ChildNum]->y, parameters.Dimension, MPI_DOUBLE, ChildNumRecv, TagChildSolved, MPI_COMM_WORLD);

    RecieveCalculatedFunctional();

    outputSet.trials[0] = vecTrials[ChildNum];
    for (int j = 0; j <= outputSet.trials[0]->index; j++)
        outputSet.countCalcTrials[j]++;
    vecTrials[ChildNum] = 0;
}


// ------------------------------------------------------------------------------------------------

void MPICalculationAsync::Calculate(InformationForCalculation& inputSet,
    TResultForCalculation& outputSet)
{
    if (inputSet.trials.size() > 0)
    {
        outputSet.trials.clear();
        outputSet.trials.resize(inputSet.trials.size());
        for (unsigned i = 0; i < outputSet.trials.size(); i++)
            outputSet.trials[i] = 0;


        outputSet.countCalcTrials.clear();
        outputSet.countCalcTrials.resize(pTask->GetNumOfFunc());
        for (int i = 0; i < pTask->GetNumOfFunc(); i++)
            outputSet.countCalcTrials[i] = 0;

    }
    if (isStartComputingAway)
    {
        if (isFirst) {
            if (parameters.DebugAsyncCalculation != 0) {
                std::ofstream fout;
                fout.open("../_build/async.txt");
                fout << parameters.GetProcNum() << "\n";
                fout << 1;
                fout.close();
            }
            FirstStartCalculate(inputSet, outputSet);
        }
        else {
            StartCalculate(inputSet, outputSet);
        }
    }
    else
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

            if (isFirst)
                FirstStartCalculate(inputCalculation, resultCalculation);
            else
                StartCalculate(inputCalculation, resultCalculation);
        }
    }

    if (parameters.DebugAsyncCalculation != 0) {
        std::ifstream fin("../_build/async.txt");
        int i;
        fin >> i;
        fin >> i;
        fin.close();

        i++;
        if (i == parameters.GetProcNum()) {
            i = 1;
        }

        std::ofstream fout;
        fout.open("../_build/async.txt");
        fout << parameters.GetProcNum() << "\n";
        fout << i;
        fout.close();
    }
}