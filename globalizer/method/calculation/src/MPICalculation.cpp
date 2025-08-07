/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2021 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      mpi_calculation.cpp                                         //
//                                                                         //
//  Purpose:   Source file for MPI calculation class                       //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "MPICalculation.h"

#include <stdlib.h>
#include <string.h>
#include <cstring>



#include "TaskFactory.h"
#include "TrialFactory.h"
#include "OmpCalculation.h"

void MPICalculation::StartCalculate(InformationForCalculation& inputSet,
  TResultForCalculation& outputSet)
{
  for (unsigned int i = 0; i < parameters.GetProcNum() - 1; i++)
  {
    int isFinish = 0;
    //Îòïðàâëÿåì â Solver ôëàã, ÷òî ìû ðàáîòàåì
    MPI_Send(&isFinish, 1, MPI_INT, i + 1, TagChildSolved, MPI_COMM_WORLD);

    //Îòïðàâëÿåì íåñêîëüêî òî÷åê íà ïðîöåññû
    for (unsigned int j = 0; j < parameters.mpiBlockSize; j++) {
      Trial* trail = inputSet.trials[i*(parameters.mpiBlockSize) + j];
      trail->index = -1;
      //Îòïðàâëÿåì êîîðäèíàòó y
      MPI_Send(trail->y, parameters.Dimension, MPI_DOUBLE, i + 1, TagChildSolved, MPI_COMM_WORLD);
    }
  }

  MPI_Status status;
  for (unsigned int i = 0; i < parameters.GetProcNum() - 1; i++)
  {
    //Ïðèíèìàåì âñå îòïðàâëåííûå òî÷êè îáðàòíî
    for (unsigned int j = 0; j < parameters.mpiBlockSize; j++) {
      Trial* trail = inputSet.trials[i*(parameters.mpiBlockSize) + j];
      trail->index = -1;

      //Ïðèíèìàåì âû÷èñëåííîå çíà÷åíèå ôóíêöèè èç Solver
      MPI_Recv(trail->FuncValues, MaxNumOfFunc, MPI_DOUBLE, i + 1, TagChildSolved, MPI_COMM_WORLD, &status);

      int fNumber = 0;
      while ((trail->index == -1) && (fNumber < pTask->GetNumOfFunc()))
      {
        if ((fNumber == (pTask->GetNumOfFunc() - 1)) || (trail->FuncValues[fNumber] > 0))
        {
          trail->index = fNumber;
        }
        fNumber++;
      }
    }
  }

  for (unsigned int i = 0; i < inputSet.trials.size(); i++)
  {
    for (int j = 0; j <= outputSet.trials[i]->index; j++)
      outputSet.countCalcTrials[j]++;
  }
}

// ------------------------------------------------------------------------------------------------
void MPICalculation::StartCalculateInBorder(InformationForCalculation& inputSet,
  TResultForCalculation& outputSet)
{
  //for (unsigned int i = 0; i < 2; i++)
  //{
  //  int isFinish = 0;
  //  //Îòïðàâëÿåì â Solver ôëàã, ÷òî ìû ðàáîòàåì
  //  MPI_Send(&isFinish, 1, MPI_INT, i + 1, TagChildSolved, MPI_COMM_WORLD);


  //  Trial* trail = inputSet.trials[i];
  //  trail->index = -1;
  //  //Îòïðàâëÿåì êîîðäèíàòó y
  //  MPI_Send(trail->y, parameters.Dimension, MPI_DOUBLE, i + 1, TagChildSolved, MPI_COMM_WORLD);
  //}

  //MPI_Status status;
  //for (unsigned int i = 0; i < 2; i++)
  //{
  //  //Ïðèíèìàåì âñå îòïðàâëåííûå òî÷êè îáðàòíî

  //  Trial* trail = inputSet.trials[i];
  //  trail->index = -1;

  //  //Ïðèíèìàåì âû÷èñëåííîå çíà÷åíèå ôóíêöèè èç Solver
  //  MPI_Recv(trail->FuncValues, MaxNumOfFunc, MPI_DOUBLE, i + 1, TagChildSolved, MPI_COMM_WORLD, &status);

  //  int fNumber = 0;
  //  while ((trail->index == -1) && (fNumber < pTask->GetNumOfFunc()))
  //  {
  //    if ((fNumber == (pTask->GetNumOfFunc() - 1)) || (trail->FuncValues[fNumber] > 0))
  //    {
  //      trail->index = fNumber;
  //    }
  //    fNumber++;
  //  }
  //}

  Calculation* calculation;
  calculation = new OMPCalculation(*pTask);
  calculation->Calculate(inputSet, outputSet);

  for (unsigned int i = 0; i < inputSet.trials.size(); i++)
  {
    for (int j = 0; j <= outputSet.trials[i]->index; j++)
      outputSet.countCalcTrials[j]++;
  }
}


// ------------------------------------------------------------------------------------------------
void MPICalculation::Calculate(InformationForCalculation& inputSet,
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

  // Çàïóñêàòü âû÷èñëåíèÿ êàê òîëüêî ïðèøëè äàííûå
  if (isStartComputingAway)
  {
    if ((isFirst) && ((parameters.isCalculationInBorderPoint == true) || (parameters.LocalTuningType != 0)))
    {
      isFirst = false;
      StartCalculateInBorder(inputSet, outputSet);
    }
    else
      StartCalculate(inputSet, outputSet);
  }
  else//ñîáðàòü äàííûå â îäèí áëîê, è ïîòîì âû÷èñëèòü âñå ñðàçó
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
