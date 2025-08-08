/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      task.cpp                                                    //
//                                                                         //
//  Purpose:   Source file for optimization task class                     //
//                                                                         //
//  Author(s): Sysoyev A., Barkalov K.                                     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "Exception.h"
#include "Task.h"
#include <cstring>

// ------------------------------------------------------------------------------------------------
Task::Task(IProblem* _problem, int _ProcLevel)
{
  if ((parameters.Dimension <= 0) || (parameters.Dimension > MaxDim))
  {
    throw EXCEPTION("N is out of range");
  }

  NumOfFunc = _problem->GetNumberOfFunctions();
  _problem->GetBounds(A, B);
  IsOptimumPointDefined = _problem->GetOptimumPoint(OptimumPoint) == IProblem::OK ? true : false;
  IsOptimumValueDefined = _problem->GetOptimumValue(OptimumValue) == IProblem::OK ? true : false;
  pProblem = _problem;

  ProcLevel = _ProcLevel;
  isInit = true;
}

Task::Task()
{
  NumOfFunc = 0;
  IsOptimumPointDefined = false;
  IsOptimumValueDefined = false;
  pProblem = 0;

  ProcLevel = 0;
  isInit = false;
}

// ------------------------------------------------------------------------------------------------
Task::~Task()
{
}

Task* Task::Clone()
{
  Task* res = 0;
  if (isInit)
    res = new Task(pProblem, ProcLevel);
  else
    res = new Task();

  res->isInit = isInit;
  return res;
}

void Task::Init(IProblem * _problem, int _ProcLevel)
{
  if ((parameters.Dimension <= 0) || (parameters.Dimension > MaxDim))
  {
    throw EXCEPTION("N is out of range");
  }

  NumOfFunc = _problem->GetNumberOfFunctions();
  _problem->GetBounds(A, B);
  IsOptimumPointDefined = _problem->GetOptimumPoint(OptimumPoint) == IProblem::OK ? true : false;
  IsOptimumValueDefined = _problem->GetOptimumValue(OptimumValue) == IProblem::OK ? true : false;
  pProblem = _problem;
  // По умолчанию фиксированные переменные отсутствуют

  ProcLevel = _ProcLevel;
  isInit = true;
}


// - end of file ----------------------------------------------------------------------------------
