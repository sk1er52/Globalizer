#include "SeparableOptimizationSolver.h"
#include "TaskFactory.h"

// ------------------------------------------------------------------------------------------------
void SeparableOptimizationSolver::SetDimentions(std::vector<int> _dimentions)
{
  if (_dimentions.size() == 0)
  {
    dimensions.resize(parameters.Dimension);
    for (int i = 0; i < dimensions.size(); i++)
      dimensions[i] = 1;
  }
  else
  {
    dimensions = _dimentions;
  }
}

// ------------------------------------------------------------------------------------------------
void SeparableOptimizationSolver::CreateStartPoint()
{
  if (parameters.startPoint.GetIsChange() == false)
  {
    double* A = new double[parameters.Dimension];
    double* B = new double[parameters.Dimension];
    problem->GetBounds(A, B);
    parameters.startPoint.SetSize(parameters.Dimension);
    for (int i = 0; i < parameters.Dimension; i++)
    {
      parameters.startPoint[i] = A[i] + (B[i] - A[i]) / 2.0;
    }
  }
}

// ------------------------------------------------------------------------------------------------
void SeparableOptimizationSolver::Construct()
{
  solvers.resize(dimensions.size());
  tasks.resize(dimensions.size());
  for (int i = 0; i < dimensions.size(); i++)
  {
    solvers[i] = new Solver(problem);
    tasks[i] = nullptr;
  }
  solutionResult = nullptr;
  originalDimension = parameters.Dimension;
  CreateStartPoint();
  parameters.TypeSolver = SeparableSearch;
}



// ------------------------------------------------------------------------------------------------
SeparableOptimizationSolver::SeparableOptimizationSolver(IProblem* _problem, std::vector<int> _dimentions)
{
  problem = _problem;
  SetDimentions(_dimentions);
  Construct();
  finalSolver = nullptr;
}

// ------------------------------------------------------------------------------------------------
#ifdef _GLOBALIZER_BENCHMARKS
SeparableOptimizationSolver::SeparableOptimizationSolver(IGlobalOptimizationProblem* _problem, std::vector<int> _dimentions) 
  : SeparableOptimizationSolver::SeparableOptimizationSolver(new GlobalizerBenchmarksProblem(_problem), _dimentions)
{
}
#endif

// ------------------------------------------------------------------------------------------------
SeparableOptimizationSolver::~SeparableOptimizationSolver()
{
  if (solutionResult == nullptr)
    delete solutionResult;
  for (int i = 0; i < solvers.size(); i++)
  {
    if (solvers[i] != nullptr)
    {
      auto solver = solvers[i];
      delete solver;
    }
    if (tasks[i] != nullptr)
    {
      auto task = tasks[i];
      delete task;
    }
  }
}

// ------------------------------------------------------------------------------------------------
int SeparableOptimizationSolver::Solve()
{
  try
  {
    auto doLV = parameters.localVerificationType;
    parameters.localVerificationType = None;
    int startParameterNumber = 0;

    std::vector<Trial*> points;
    double bestVolue = MaxDouble;
    for (int i = 0; i < solvers.size(); i++)
    {
      parameters.Dimension = dimensions[i];

      Solver* solver = solvers[i];
      if (tasks[i] != nullptr)
        delete tasks[i];
      tasks[i] = dynamic_cast<SeparableOptimizationTask*>(TaskFactory::CreateTask(problem, 0));      

      tasks[i]->SetStartParameterNumber(startParameterNumber);

      solver->Solve(tasks[i]);

      auto solution = solver->GetSolutionResult();
      
      if (solution->BestTrial->index == problem->GetNumberOfConstraints())
      {
#ifdef WIN32
        if (_finite(solution->BestTrial->GetValue()) != 0)
#else
        if (std::isfinite(solution->BestTrial->GetValue()) != 0)
#endif
        {
          if (bestVolue >= solution->BestTrial->GetValue())
          {
            bestVolue = solution->BestTrial->GetValue();
            for (int j = 0; j < dimensions[i]; j++)
            {
              parameters.startPoint[j + startParameterNumber] = solution->BestTrial->y[j];
            }

            parameters.startPointValues.SetSize(tasks[i]->GetNumOfFunc());
            for (int j = 0; j < tasks[i]->GetNumOfFunc(); j++)
            {
              parameters.startPointValues[j] = solution->BestTrial->FuncValues[j];
            }
          }
        }
      }
      startParameterNumber = startParameterNumber + parameters.Dimension;

      parameters.Dimension = originalDimension;

      points.insert(points.end(), solver->GetAllPoint().begin(), solver->GetAllPoint().end());

      Calculation::leafCalculation = 0;  

      print << "bestVolue =\t" << bestVolue << "\t Dimension " << i << "\n";
      print << "best coordinate = \t" << parameters.startPoint.ToString() << "\n";
    }

    parameters.localVerificationType = doLV;
    if (finalSolver == nullptr)
      finalSolver = new Solver(problem);
    else
    {

      finalSolver = new Solver(problem);
    }
    parameters.MaxNumOfPoints[0] = 2;

    //finalSolver->SetPoint(points);

    finalSolver->Solve();

    if (solutionResult != nullptr)
      delete solutionResult;
    solutionResult = finalSolver->GetSolutionResult();

  }
  catch (const Exception& e)
  {
    std::string excFileName = std::string("exception_") +
      toString(parameters.GetProcRank()) + ".txt";
    e.Print(excFileName.c_str());

    for (int i = 0; i < parameters.GetProcNum(); i++)
      if (i != parameters.GetProcRank())
        MPI_Abort(MPI_COMM_WORLD, i);
    return 1;
  }
  catch (...)
  {
    print << "\nUNKNOWN EXCEPTION !!!\n";
    std::string excFileName = std::string("exception_") +
      toString(parameters.GetProcRank()) + ".txt";
    Exception e("UNKNOWN FILE", -1, "UNKNOWN FUCNTION", "UNKNOWN EXCEPTION");
    e.Print(excFileName.c_str());

    for (int i = 0; i < parameters.GetProcNum(); i++)
      if (i != parameters.GetProcRank())
        MPI_Abort(MPI_COMM_WORLD, i);
    return 1;
  }
  if (parameters.GetProcRank() == 0)
  {
    if (parameters.GetProcNum() > 1) {
      int childNum = parameters.GetProcNum() - 1;
      int curr_child = 0;
      for (unsigned int i = 0; i < childNum; i++) {
        ///curr_child = parameters.parallel_tree.ProcChild[i];!!!!!
        int finish = 1;
        MPI_Send(&finish, 1, MPI_INT, curr_child, TagChildSolved, MPI_COMM_WORLD);
      }
    }
  }

  return 0;
}

// ------------------------------------------------------------------------------------------------
SolutionResult* SeparableOptimizationSolver::GetSolutionResult()
{
  return solutionResult;
}

// ------------------------------------------------------------------------------------------------
void SeparableOptimizationSolver::SetPoint(std::vector<Trial*>& points)
{
  finalSolver->SetPoint(points);
}

// ------------------------------------------------------------------------------------------------
std::vector<Trial*>& SeparableOptimizationSolver::GetAllPoint()
{
  return finalSolver->GetAllPoint();
}
