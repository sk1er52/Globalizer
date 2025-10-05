/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      process.cpp                                                 //
//                                                                         //
//  Purpose:   Source file for optimization process class                  //
//                                                                         //
//  Author(s): Sysoyev A., Barkalov K.                                     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
#include <mpi.h>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <stdio.h>
#include <stdlib.h>

#include "Common.h"
#include "GlobProcess.h"
#include "MethodFactory.h"
#include "MPICalculationAsync.h"


#include "SearcDataIterator.h"

#include "TaskFactory.h"


#define DEBUG_MPI

// ------------------------------------------------------------------------------------------------
void ShowIterResults(Process *pProcess)
{
}

// ------------------------------------------------------------------------------------------------
Process::Process(SearchData& data, Task& task) :
pData(&data), pTask(&task)
{
  isFirstRun = true;

  if (parameters.MapType == mpBase)
    evolvent = new Evolvent(parameters.Dimension - pTask->GetNumberOfDiscreteVariable(), parameters.m);
  else
    throw EXCEPTION("Unknown type of evolvent");

  calculation = CalculationFactory::CreateCalculation(*pTask, evolvent);

  pMethod = MethodFactory::CreateMethod(*pTask, *pData, *calculation, *evolvent);

  functionCalculationCount.resize(pTask->GetNumOfFunc());

  isPrintOptimEstimation = false;

  Neighbours.resize(parameters.GetProcNum() - 1);
  for (int i = 0; i < parameters.GetProcNum() - 1; i++)
    Neighbours[i] = i + 1;
  addPoints = nullptr;
}

// ------------------------------------------------------------------------------------------------
Process::~Process()
{
  delete calculation;
  delete pMethod;
  delete evolvent;
}

// ------------------------------------------------------------------------------------------------
double Process::GetSolveTime()
{
  return duration;
}

// ------------------------------------------------------------------------------------------------
void Process::Solve()
{
  Trial OptimEstimation, NewOptimEstimation;

  Timer.Start();
  if (isFirstRun)
  {
    BeginIterations();
    isFirstRun = false;
  }
  OptimEstimation = *(pMethod->GetOptimEstimation());
 
  while (!IsOptimumFound)
  {

    if (!pTask->IsLeaf())
    {
      if (!(pMethod->GetIterationCount() % parameters.StepPrintMessages))
      {
        printf("process 0, iteration %d\t point count = %d\n", pMethod->GetIterationCount(), pData->GetCount());
        Trial* p = (pMethod->GetOptimEstimation());
        if (p->index >= 0)
        {
          printf("  Cur min = %lf \n", p->FuncValues[p->index]);
          for (int i = 0; i < pTask->GetN(); i++)
          {
            printf("  Cur x[%d] = %lf \n", i, p->y[i]);
          }
        }
      }

      if (static_cast<std::string>(parameters.iterPointsSavePath).size())
      {
        if (!(pMethod->GetIterationCount() % parameters.StepSavePoint))
        {
          pMethod->PrintPoints(std::to_string(pMethod->GetIterationCount()) + "_" + parameters.iterPointsSavePath.ToString());
        }
      }
    }
    pMethod->GetIterationCount();

    DoIteration();

    NewOptimEstimation = *(pMethod->GetOptimEstimation());
    if (NewOptimEstimation.FuncValues[NewOptimEstimation.index] !=
      OptimEstimation.FuncValues[OptimEstimation.index])
    {
      OptimEstimation = NewOptimEstimation;
      duration = Timer.GetTime();
      PrintOptimEstimationToFile(OptimEstimation);
    }

  }



  EndIterations();

  duration = Timer.GetTime();

  OptimEstimation = *(pMethod->GetOptimEstimation());

  if (!pTask->IsLeaf() && isPrintOptimEstimation)
  {
    printf("\n");
    PrintOptimEstimationToFile(OptimEstimation);
    PrintOptimEstimationToConsole(OptimEstimation);
    PrintResultToFile(OptimEstimation);
  }
}

void Process::Reset(SearchData* data, Task* task)
{
  isFirstRun = true;
  pData = data;
  pTask = task;
  if (evolvent == 0)
  {
    if (parameters.MapType == mpBase)
      evolvent = new Evolvent(parameters.Dimension, parameters.m);
    else
      throw EXCEPTION("Unknown type of evolvent");
  }
  if (calculation == 0)
    calculation = CalculationFactory::CreateCalculation(*pTask, evolvent);
  else
    calculation->SetTask(pTask);

  if (pMethod != 0)
    delete pMethod;
  pMethod = MethodFactory::CreateMethod(*pTask, *pData, *calculation, *evolvent);

  functionCalculationCount.clear();
  functionCalculationCount.resize(pTask->GetNumOfFunc());
}

// ------------------------------------------------------------------------------------------------
void Process::PrintOptimEstimationToFile(Trial OptimEstimation)
{
  int NumberOfTrials;
  FILE* pf;

  if (parameters.IsPrintFile && !pTask->IsLeaf())
  {
    pf = fopen("optim.dat", "a");

    fprintf(pf, "Iteration = %d\n", pMethod->GetIterationCount());
    fprintf(pf, "min = %lf \n", OptimEstimation.FuncValues[OptimEstimation.index]);
    for (int i = 0; i < pTask->GetN(); i++)
      fprintf(pf, "x[%d] = %lf \n", i, OptimEstimation.y[i]);

    fprintf(pf, "\n");
    for (int i = 0; i < pTask->GetNumOfFunc(); i++)
      fprintf(pf, "M[%d] = %lf \n", i, (this->pData->M)[i]);

    NumberOfTrials = pMethod->GetNumberOfTrials();
    fprintf(pf, "NumberOfTrials = %d\n", NumberOfTrials);
    fprintf(pf, "Solve time = %f \n\n\n", duration);

    fclose(pf);


  }
}

// ------------------------------------------------------------------------------------------------
void Process::PrintOptimEstimationToConsole(Trial OptimEstimation)
{
  int NumberOfTrials;
  double ValueDifference = HUGE_VAL;
  double PointDifference = HUGE_VAL;
  int numOfOptima = 1;
  double* allOptimumPoints = new double[MAX_TRIAL_DIMENSION*MAX_NUM_MIN];
  double* allPointDifference = new double[MAX_NUM_MIN];

  printf("ProcLevel = %d\n", pTask->GetProcLevel());
  printf("Iteration = %d \n", pMethod->GetIterationCount());
  printf("Point = %d \n", pMethod->GetIterationCount() * parameters.NumPoints);

  printf("min = %lf \n", OptimEstimation.FuncValues[OptimEstimation.index]);
  for (int i = 0; i < pTask->GetN(); i++)
  {
    printf("x[%d] = %lf \n", i, OptimEstimation.y[i]);
  }

  printf("\n");
  printf("constants = ");
  for (int i = 0; i < pTask->GetN(); i++)
  {
    printf("%lf,_", OptimEstimation.y[i]);
  }

  printf("\n");

  if (pTask->GetIsOptimumValueDefined())
  {
    ValueDifference = OptimEstimation.FuncValues[OptimEstimation.index] - pTask->GetOptimumValue();
    printf("\nValue difference  = %lf \n", ValueDifference);
  }

  if (pTask->GetIsOptimumPointDefined())
  {
    PointDifference = 0.0;
    for (int i = 0; i < pTask->GetN(); i++)
    {
      PointDifference = GLOBALIZER_MAX(fabs(OptimEstimation.y[i] - pTask->GetOptimumPoint()[i]), PointDifference);
    }
    if (pTask->getProblem()->GetAllOptimumPoint(allOptimumPoints, numOfOptima) != IProblem::UNDEFINED)
    {
      PointDifference = 0.0;
      for (int j = 0; j < numOfOptima; j++)
      {
        allPointDifference[j] = 0.0;
        for (int i = 0; i < pTask->GetN(); i++)
        {
          allPointDifference[j] = GLOBALIZER_MAX(fabs(OptimEstimation.y[i] -
            allOptimumPoints[j*pTask->GetN() + i]), allPointDifference[j]);
        }
        if ((allPointDifference[j] < PointDifference) || (j == 0))
        {
          PointDifference = allPointDifference[j];
        }
      }
    }
    printf("Coordinates max difference = %lf \n", PointDifference);
  }

  //if (pTask->GetIsOptimumPointDefined() || pTask->GetIsOptimumValueDefined())
  //{
  //  printf("Global optimum %s\n", (PointDifference < parameters.Epsilon ||
  //                                 ValueDifference < parameters.Epsilon) ? "FOUND!" : "NOT FOUND" );
  //}


  bool res = false;

  switch (1)
  {
  case Accuracy:
    if (pMethod->GetAchievedAccuracy() < parameters.Epsilon)
      res = true;
    break;
  case OptimumVicinity:
    {
      res = true;
      //numOfOptima = ;
      if (pTask->getProblem()->GetAllOptimumPoint(allOptimumPoints, numOfOptima) == IProblem::UNDEFINED)
      {
        for (int i = 0; i < parameters.Dimension; i++)
        {

          double fabsx = fabs(OptimEstimation.y[i] - pTask->GetOptimumPoint()[i]);
          double fm = parameters.Epsilon * (pTask->GetB()[i] - pTask->GetA()[i]);
          if (fabsx > fm)
          {
            res = res && false;
          }
        }
      }
      else
      {
        for (int j = 0; j < numOfOptima; j++)
        {
          for (int i = 0; i < parameters.Dimension; i++)
          {
            double fabsx = fabs(OptimEstimation.y[i] - allOptimumPoints[parameters.Dimension*j + i]);
            double fm = parameters.Epsilon * (pTask->GetB()[i] - pTask->GetA()[i]);
            if (fabsx > fm)
            {
              res = res && false;
              break;
            }
            if (i == parameters.Dimension - 1)
            {
              res = true;
            }
          }
          if (res == true)
          {
            break;
          }
        }
      }
    }
    break;
  case OptimumVicinity2:
    {
      res = true;
      for (int i = 0; i < pTask->GetN(); i++)
      {
        if (fabs(OptimEstimation.y[i] - pTask->GetOptimumPoint()[i]) > parameters.Epsilon)
        {
          res = false;
          break;
        }
      }
    }
    break;
  case OptimumValue:
    if (OptimEstimation.index == pTask->GetNumOfFunc() - 1 &&
      OptimEstimation.FuncValues[OptimEstimation.index] - pTask->GetOptimumValue() <
      parameters.Epsilon)
      res = true;
    break;
  }
  delete[] allOptimumPoints;

  printf("Global optimum %s\n", (res) ? "FOUND!" : "NOT FOUND");

//  /************************************************/
//
//  std::string fileName;
//  if (parameters.TypeCalculation != OMP) {
//    fileName = "seria_full_result.txt";
//
//    std::ofstream vmdelet_out;                    //создаем поток 
//    vmdelet_out.open(fileName, std::ios::app);  // открываем файл для записи в конец
//
//    if (parameters.TypeCalculation != OMP) {
//      vmdelet_out << "Global: " << parameters.DimInTaskLevel[0] << "\nLocal: " << parameters.DimInTaskLevel[1] << "\n";
//
//      vmdelet_out << "Permutations vector: ";
//      for (auto const& element : TaskFactory::permutations) {
//        vmdelet_out << element << " ";
//      }
//      vmdelet_out << "\n\n";
//    }
//
//    vmdelet_out << "Min: " << OptimEstimation.FuncValues[OptimEstimation.index] << " [ ";
//    for (int i = 0; i < pTask->GetN(); i++)
//    {
//      vmdelet_out << OptimEstimation.y[i] << " ";
//    }
//    vmdelet_out << "]\nREAL Min: " << pTask->GetOptimumValue() << " [ ";
//    for (int i = 0; i < pTask->GetN(); i++)
//    {
//      vmdelet_out << pTask->GetOptimumPoint()[i] << " ";
//    }
//    vmdelet_out << "]\n";
//
//    vmdelet_out << "Coordinates max difference: " << PointDifference << "\n";
//    ValueDifference = OptimEstimation.FuncValues[OptimEstimation.index] - pTask->GetOptimumValue();
//    vmdelet_out << "Value difference: " << ValueDifference << "\n";
//
//    vmdelet_out << "Iteration: " << pMethod->GetIterationCount() << "\n";
//    vmdelet_out << "Point: " << pMethod->GetIterationCount() * parameters.NumPoints << "\n";
//
//    vmdelet_out << "Global optimum " << ((res) ? "FOUND!" : "NOT FOUND") << "\n-------------------------------------\n";
//    vmdelet_out.close();                          // закрываем файл
//  }
//  /************************************************/
//
  NumberOfTrials = pMethod->GetNumberOfTrials();
  pMethod->PrintSection();

  printf("\nNumberOfTrials = %d\n", NumberOfTrials);
  for (int i = 0; i < pTask->GetNumOfFunc(); i++)
    printf("Number of calculations function %d of = %d\n", i,
    pMethod->GetFunctionCalculationCount()[i]);

  printf("\nLocalPointCount = %d\n", pMethod->GetLocalPointCount());
  printf("\nNumberLocalMethodtStart = %d\n", pMethod->GetNumberLocalMethodtStart());
  //if (parameters.printAdvancedInfo)
  //{
  //  printf("\n");
  //  for (int i = 0; i < pTask->GetNumOfFunc(); i++)
  //    printf("M estimation for function %d = %f\n", i, pMethod->GetM()[i]);
  //  printf("\nOptimum index: %i\n\n", OptimEstimation.index);
  //  for (int i = 0; i <= OptimEstimation.index; i++)
  //    printf("Function %d value in estimated optimum: %f\n", i, OptimEstimation.FuncValues[i]);
  //}
  printf("\nSolve time = %lf\n\n\n", duration);
//
//
//  if (parameters.TypeCalculation != OMP) {
//    fileName = "seria_table_result.txt";
//
//    std::ofstream statOut;                    //создаем поток 
//    statOut.open(fileName, std::ios::app);  // открываем файл для записи в конец
//
//    statOut << ((res) ? "FOUND" : "NOT_FOUND") << " "
//      << pMethod->GetIterationCount() << " "
//      << NumberOfTrials << " "
//      << pMethod->GetLocalPointCount() << " "
//      << duration << " "
//      << parameters.DimInTaskLevel[0] << " ";
//    for (int ii = 0; ii < parameters.DimInTaskLevel[0]; ii++) {
//       statOut << TaskFactory::permutations[ii] + 1 << "_";
//    }
//    statOut << "\n";
//
//    statOut.close();                          // закрываем файл
//  }
}


// ------------------------------------------------------------------------------------------------
void Process::PrintResultToFile(Trial OptimEstimation)
{
  FILE* pf;

  if (parameters.ResulLog.ToString() != "000" && pTask->GetProcLevel() == 0)
  {
    pf = fopen(parameters.ResulLog.ToString().c_str(), "a");

    parameters.PrintParametersToFile(pf);

    int NumberOfTrials;
    double ValueDifference = HUGE_VAL;
    double PointDifference = HUGE_VAL;
    int numOfOptima = 1;
    double* allOptimumPoints = new double[MAX_TRIAL_DIMENSION * MAX_NUM_MIN];
    double* allPointDifference = new double[MAX_NUM_MIN];

    fprintf(pf, "ProcLevel = %d\n", pTask->GetProcLevel());
    fprintf(pf, "Iteration = %d \n", pMethod->GetIterationCount());
    fprintf(pf, "Point = %d \n", pMethod->GetIterationCount() * parameters.NumPoints);

    fprintf(pf, "min = %lf \n", OptimEstimation.FuncValues[OptimEstimation.index]);
    for (int i = 0; i < pTask->GetN(); i++)
    {
      fprintf(pf, "x[%d] = %lf \n", i, OptimEstimation.y[i]);
    }

    fprintf(pf, "\n");
    fprintf(pf, "constants = ");
    for (int i = 0; i < pTask->GetN(); i++)
    {
      fprintf(pf, "%lf,_", OptimEstimation.y[i]);
    }

    fprintf(pf, "\n");

    if (pTask->GetIsOptimumValueDefined())
    {
      ValueDifference = OptimEstimation.FuncValues[OptimEstimation.index] - pTask->GetOptimumValue();
      fprintf(pf, "\nValue difference  = %lf \n", ValueDifference);
    }

    if (pTask->GetIsOptimumPointDefined())
    {
      PointDifference = 0.0;
      for (int i = 0; i < pTask->GetN(); i++)
      {
        PointDifference = GLOBALIZER_MAX(fabs(OptimEstimation.y[i] - pTask->GetOptimumPoint()[i]), PointDifference);
      }
      if (pTask->getProblem()->GetAllOptimumPoint(allOptimumPoints, numOfOptima) != IProblem::UNDEFINED)
      {
        PointDifference = 0.0;
        for (int j = 0; j < numOfOptima; j++)
        {
          allPointDifference[j] = 0.0;
          for (int i = 0; i < pTask->GetN(); i++)
          {
            allPointDifference[j] = GLOBALIZER_MAX(fabs(OptimEstimation.y[i] -
              allOptimumPoints[j * pTask->GetN() + i]), allPointDifference[j]);
          }
          if ((allPointDifference[j] < PointDifference) || (j == 0))
          {
            PointDifference = allPointDifference[j];
          }
        }
      }
      fprintf(pf, "Coordinates max difference = %lf \n", PointDifference);
    }

    //if (pTask->GetIsOptimumPointDefined() || pTask->GetIsOptimumValueDefined())
    //{
    //  fprintf(pf, "Global optimum %s\n", (PointDifference < parameters.Epsilon ||
    //                                 ValueDifference < parameters.Epsilon) ? "FOUND!" : "NOT FOUND" );
    //}


    bool res = false;

    switch (1)
    {
    case Accuracy:
      if (pMethod->GetAchievedAccuracy() < parameters.Epsilon)
        res = true;
      break;
    case OptimumVicinity:
    {
      res = true;
      //numOfOptima = ;
      if (pTask->getProblem()->GetAllOptimumPoint(allOptimumPoints, numOfOptima) == IProblem::UNDEFINED)
      {
        for (int i = 0; i < parameters.Dimension; i++)
        {

          double fabsx = fabs(OptimEstimation.y[i] - pTask->GetOptimumPoint()[i]);
          double fm = parameters.Epsilon * (pTask->GetB()[i] - pTask->GetA()[i]);
          if (fabsx > fm)
          {
            res = res && false;
          }
        }
      }
      else
      {
        for (int j = 0; j < numOfOptima; j++)
        {
          for (int i = 0; i < parameters.Dimension; i++)
          {
            double fabsx = fabs(OptimEstimation.y[i] - allOptimumPoints[parameters.Dimension * j + i]);
            double fm = parameters.Epsilon * (pTask->GetB()[i] - pTask->GetA()[i]);
            if (fabsx > fm)
            {
              res = res && false;
              break;
            }
            if (i == parameters.Dimension - 1)
            {
              res = true;
            }
          }
          if (res == true)
          {
            break;
          }
        }
      }
    }
    break;
    case OptimumVicinity2:
    {
      res = true;
      for (int i = 0; i < pTask->GetN(); i++)
      {
        if (fabs(OptimEstimation.y[i] - pTask->GetOptimumPoint()[i]) > parameters.Epsilon)
        {
          res = false;
          break;
        }
      }
    }
    break;
    case OptimumValue:
      if (OptimEstimation.index == pTask->GetNumOfFunc() - 1 &&
        OptimEstimation.FuncValues[OptimEstimation.index] - pTask->GetOptimumValue() <
        parameters.Epsilon)
        res = true;
      break;
    }
    delete[] allOptimumPoints;

    fprintf(pf, "Global optimum %s\n", (res) ? "FOUND!" : "NOT FOUND");


    NumberOfTrials = pMethod->GetNumberOfTrials();
    pMethod->PrintSection();

    fprintf(pf, "\nNumberOfTrials = %d\n", NumberOfTrials);
    for (int i = 0; i < pTask->GetNumOfFunc(); i++)
      fprintf(pf, "Number of calculations function %d of = %d\n", i,
        pMethod->GetFunctionCalculationCount()[i]);

    fprintf(pf, "\nLocalPointCount = %d\n", pMethod->GetLocalPointCount());
    fprintf(pf, "\nNumberLocalMethodtStart = %d\n", pMethod->GetNumberLocalMethodtStart());
    //if (parameters.printAdvancedInfo)
    //{
    //  fprintf(pf, "\n");
    //  for (int i = 0; i < pTask->GetNumOfFunc(); i++)
    //    fprintf(pf, "M estimation for function %d = %f\n", i, pMethod->GetM()[i]);
    //  fprintf(pf, "\nOptimum index: %i\n\n", OptimEstimation.index);
    //  for (int i = 0; i <= OptimEstimation.index; i++)
    //    fprintf(pf, "Function %d value in estimated optimum: %f\n", i, OptimEstimation.FuncValues[i]);
    //}
    fprintf(pf, "\nSolve time = %lf\n\n\n", duration);
    fclose(pf);
  }
}


// ------------------------------------------------------------------------------------------------
void Process::BeginIterations()
{
  IsOptimumFound = false;
  isPrintOptimEstimation = false;

  pMethod->FirstIteration();

  if (addPoints != nullptr)
    pMethod->InsertPoints(*addPoints);
}

// ------------------------------------------------------------------------------------------------
void Process::DoIteration()
{
  bool IsStop;

  ////if (parameters.GetProcRank() == 0) printf("AA1! %d\n", parameters.GetProcRank());

  //  проверяем критерий остановки
  IsStop = pMethod->CheckStopCondition();
  IsOptimumFound = IsStop;

  if (IsOptimumFound)
    isPrintOptimEstimation = true;

  IsOptimumFound = CheckIsStop(IsOptimumFound);
  IsStop = IsOptimumFound;

  //if (pTask->GetProcLevel() == 0)  printf("CC-1A! Proc=%d iter=%d PL=%d IsStop=%d\n", parameters.GetProcRank(), this->pMethod->GetIterationCount(), (parameters.MapProcInLevel[pTask->GetProcLevel()] - 1), IsStop);

  try
  {
    if (!IsStop)
    {
      // вычисляем координаты испытаний
      pMethod->CalculateIterationPoints();

      IsStop = pMethod->CheckStopCondition();
      if (IsStop)
      {
        IsOptimumFound = true;
        isPrintOptimEstimation = true;
      }

      pMethod->CalculateFunctionals();

      //if (pTask->GetProcLevel() == 0)
      //{
      //  print << "IterationCount = " << pMethod->GetIterationCount() << "\n";
      //}

      for (int j = 0; j < pTask->GetNumOfFunc(); j++)
      {
        functionCalculationCount[j] += pMethod->GetFunctionCalculationCount()[j];
      }
      ////if (parameters.GetProcRank() == 0) printf("AA9! %d\n", parameters.GetProcRank());
      //Провести оценку оптимума нужно до обновления данных,  т.к. в этой функции может быть поднят флаг recalc
      pMethod->EstimateOptimum();
      //Все случаи поднятия флага recalc обработаны, можно обновлять базу
      ////if (parameters.GetProcRank() == 0) printf("AA10! %d\n", parameters.GetProcRank());
      pMethod->RenewSearchData();
      //на первой итерации проводим сепарабельную оптимизацию
      //if (pMethod->GetNumberOfTrials() == 1)
      //{
      //  if (parameters.sepS)
      //    pMethod->SeparableSearch();
      //  if (parameters.rndS)
      //    pMethod->RandomSearh();
      //}
      ////if (parameters.GetProcRank() == 0)     printf("AA20! %d\n", parameters.GetProcRank());
      //if ((parameters.localVerificationType == IntegratedOnePoint ||
      //  parameters.localVerificationType == IntegratedManyPoints) && isNewOptimumFound)
      //  pMethod->LocalSearch();

      pMethod->FinalizeIteration();

    }
    else
    {
      IsOptimumFound = true;
      if (parameters.TypeCalculation == AsyncMPI)
        MPICalculationAsync::AsyncFinilize();  
      pMethod->FinalizeIteration();
    }
  }
  catch (const Exception & e)
  {
    
    if (std::string(e.GetDescription()) == std::string("Point is outside the interval !"))
    {
      IsOptimumFound = true;
      isPrintOptimEstimation = true;
    }
    else
      throw e;
  }
}

// ------------------------------------------------------------------------------------------------
bool Process::CheckIsStop(bool IsStop)
{
  bool f = IsStop;

  MPI_Status status;

  if (pTask->GetProcLevel() == 0)
  {
      for (int k = 0; k < Neighbours.size(); k++)
      {
        MPI_Send(&f, 1, MPI_CHAR, Neighbours[k], TagChildSolved, MPI_COMM_WORLD);
      }
      //цикл по нашим соседям
      for (int k = 0; k < Neighbours.size(); k++)
      {
        MPI_Recv(&f, 1, MPI_CHAR, Neighbours[k], TagChildSolved, MPI_COMM_WORLD, &status);

        if (isPrintOptimEstimation && f && (parameters.GetProcRank() > Neighbours[k]))
          isPrintOptimEstimation = false;

        IsStop = IsStop || f;
      }
    

  }

  return IsStop;
}


// ------------------------------------------------------------------------------------------------
void Process::EndIterations()
{
  ////Локальное уточнение найденного решения
  pMethod->LocalSearch();
  //Запись в файл точек испытаний
  //if (ProcLevel == 0 && static_cast<std::string>(parameters.iterPointsSavePath).size())
  //  pMethod->PrintPoints(parameters.iterPointsSavePath);
  if (pTask->GetProcLevel() == 0 && static_cast<std::string>(parameters.iterPointsSavePath).size())
    pMethod->PrintPoints(parameters.iterPointsSavePath);
  else
    pMethod->SavePoints();
}

// ------------------------------------------------------------------------------------------------
void Process::InsertPoints(std::vector<Trial*>& points)
{
  addPoints = &points;
}

// - end of file ----------------------------------------------------------------------------------
