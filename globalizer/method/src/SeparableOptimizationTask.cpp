#include "SeparableOptimizationTask.h"

// ------------------------------------------------------------------------------------------------
SeparableOptimizationTask::SeparableOptimizationTask(IProblem* _problem, int _ProcLevel) : Task::Task(_problem, _ProcLevel)
{
  startParameterNumber = 0;
}

// ------------------------------------------------------------------------------------------------
SeparableOptimizationTask::SeparableOptimizationTask() : Task::Task()
{
  startParameterNumber = 0;
}

Task* SeparableOptimizationTask::Clone()
{
  SeparableOptimizationTask* res = 0;
  if (isInit)
    res = new SeparableOptimizationTask(pProblem, ProcLevel);
  else
    res = new SeparableOptimizationTask();

  res->startParameterNumber = startParameterNumber;
  res->isInit = isInit;
  return res;
}

// ------------------------------------------------------------------------------------------------
const double* SeparableOptimizationTask::GetA() const
{ 
  return &(A[startParameterNumber]);
}

// ------------------------------------------------------------------------------------------------
const double* SeparableOptimizationTask::GetB() const
{
  return &(B[startParameterNumber]);
}

// ------------------------------------------------------------------------------------------------
const double* SeparableOptimizationTask::GetOptimumPoint() const
{ 
  return &(OptimumPoint[startParameterNumber]);
}

// ------------------------------------------------------------------------------------------------
double SeparableOptimizationTask::CalculateFuncs(const double* y, int fNumber)
{
  double* point = new double[parameters.startPoint.GetSize()];
  for (int i = 0; i < parameters.startPoint.GetSize(); i++)
  {
    point[i] = parameters.startPoint[i];
  }
  for (int i = 0; i < parameters.Dimension; i++)
  {
    point[i + startParameterNumber] = y[i];
  }
  double multInLevel = parameters.functionSignMultiplier[GetProcLevel()];
  double result = multInLevel * pProblem->CalculateFunctionals(point, fNumber);
  return result;
}

// ------------------------------------------------------------------------------------------------
void SeparableOptimizationTask::CalculateFuncsInManyPoints(double* y, int fNumber, int numPoints, double* values)
{
  throw "Not implemented";
  //IGPUProblem* newProblem = dynamic_cast<IGPUProblem*>(pProblem);
  //if (newProblem != 0)
  //{
  //  newProblem->CalculateFunctionals(y, fNumber, numPoints, values);
  //}
}

// ------------------------------------------------------------------------------------------------
void SeparableOptimizationTask::SetStartParameterNumber(int _startParameterNumber)
{
  startParameterNumber = _startParameterNumber;
}
