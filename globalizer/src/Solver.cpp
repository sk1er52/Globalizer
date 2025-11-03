#include <cstring>

#include "Solver.h"



#include "MethodFactory.h"

#include "TaskFactory.h"

#include "TrialFactory.h"

#include "CalculationFactory.h"

#include "OMPCalculation.h"
#include "CudaCalculation.h"

#ifdef USE_PYTHON

#include "ml_wrapper.h"
#include <numpy/arrayobject.h>

#endif

// ------------------------------------------------------------------------------------------------
void Solver::ClearData()
{
  if (pTask != 0 && isExternalTask == false)
  {
    delete pTask;
    pTask = nullptr;
  }

  if (pData != 0)
  {
    delete pData;
    pData = nullptr;
  }

  mProcess = nullptr;
}

#ifdef USE_PYTHON
PyObject* makeFloatList(const double* array, int size)
{
  PyObject* l = PyList_New(size);
  for (int i = 0; i != size; i++)
    PyList_SET_ITEM(l, i, PyFloat_FromDouble(array[i]));
  return l;
}
#endif

// ------------------------------------------------------------------------------------------------
Solver::Solver(IProblem* problem)
{
  mProblem = problem;

  mProcess = 0;

  pTask = 0;
  pData = 0;

  result = 0;

  isExternalTask = false;

  addPoints = nullptr;
}

#ifdef _GLOBALIZER_BENCHMARKS
// ------------------------------------------------------------------------------------------------
Solver::Solver(IGlobalOptimizationProblem* problem) : Solver::Solver(new GlobalizerBenchmarksProblem(problem))
{
}
#endif

// ------------------------------------------------------------------------------------------------
int Solver::CheckParameters()
{
  double optimumValue;
  if (mProblem->GetOptimumValue(optimumValue) == IProblem::UNDEFINED &&
    parameters.stopCondition == OptimumValue)
  {
    print << "Stop by reaching optimum value is unsupported by this problem\n";
    return 1;
  }

  double optimumPoint[MaxDim * MaxNumOfGlobMinima];
  int n;
  if (mProblem->GetOptimumPoint(optimumPoint) == IProblem::UNDEFINED &&
    mProblem->GetAllOptimumPoint(optimumPoint, n) == IProblem::UNDEFINED &&
    parameters.stopCondition == OptimumVicinity)
  {
    print << "Stop by reaching optimum vicinity is unsupported by this problem\n";
    return 1;
  }

  return 0;
}


// ------------------------------------------------------------------------------------------------
void Solver::MpiCalculation()
{
  int isFinish = 0;
  while (isFinish == 0)
  {
    MPI_Status status;
    // Принимаем данные из mpi_calculation
    // Проверяем, что ещё работаем
    MPI_Recv(&isFinish, 1, MPI_INT, 0, TagChildSolved, MPI_COMM_WORLD, &status);
    if (isFinish == 1)
      break;

    /// Входные данные для вычислителя, формируются в CalculateFunctionals()
    InformationForCalculation inputSet;
    /// Выходные данные вычислителя, обрабатываются в CalculateFunctionals()
    TResultForCalculation outputSet;

    Trial* trail = TrialFactory::CreateTrial();

    inputSet.Resize(parameters.mpiBlockSize);
    outputSet.Resize(parameters.mpiBlockSize);

    for (unsigned int j = 0; j < parameters.mpiBlockSize; j++) {
      inputSet.trials[j] = TrialFactory::CreateTrial();
      // Получаем координаты точки
      MPI_Recv(trail->y, parameters.Dimension, MPI_DOUBLE, 0, TagChildSolved, MPI_COMM_WORLD, &status);
      trail->index = -1;

      for (int k = 0; k < parameters.Dimension; k++)
        inputSet.trials[j]->y[k] = trail->y[k];

    }

    IProblem* _problem = mProblem;
    Task* _pTask = TaskFactory::CreateTask(_problem, 0);;
    Calculation* calculation;
    if (parameters.calculationsArray[1] == OMP) {
      calculation = new OMPCalculation(*_pTask);
    }
    else if (parameters.calculationsArray[1] == CUDA) {
      calculation = new CUDACalculation(*_pTask);
    }

    calculation->Calculate(inputSet, outputSet);

    for (unsigned int j = 0; j < parameters.mpiBlockSize; j++) {
      // Отправляем обратно значение функции
      MPI_Send(inputSet.trials[j]->FuncValues, MaxNumOfFunc, MPI_DOUBLE, 0, TagChildSolved, MPI_COMM_WORLD);
    }
  }
}

// ------------------------------------------------------------------------------------------------
void Solver::AsyncCalculation()
{
  int isFinish = 0;
  while (isFinish == 0)
  {
    MPI_Status status;
    // Принимаем данные из mpi_calculation
    // Проверяем, что ещё работаем
    MPI_Recv(&isFinish, 1, MPI_INT, 0, TagChildSolved, MPI_COMM_WORLD, &status);
    if (isFinish == 1)
      break;

    Trial* trail = TrialFactory::CreateTrial();

    // Получаем координаты точки
    MPI_Recv(trail->y, parameters.Dimension, MPI_DOUBLE, 0, TagChildSolved, MPI_COMM_WORLD, &status);
    trail->index = -1;

    int fNumber = 0;

    IProblem* _problem = mProblem;
    Task* _pTask = TaskFactory::CreateTask(_problem, 0);

    if (parameters.DebugAsyncCalculation != 0) {
      while (true) {
#ifdef WIN32
        Sleep(5);
#endif
        std::ifstream fin("../_build/async.txt");
        int i;
        fin >> i;
        fin >> i;
        if (i == parameters.GetProcRank() || i == parameters.GetProcNum()) {
          fin.close();
          break;
        }
      }
    }

    // Вычисляем значение функции
    while ((trail->index == -1) && (fNumber < _pTask->GetNumOfFunc()))
    {
      trail->FuncValues[fNumber] = _pTask->CalculateFuncs(trail->y, fNumber);

#ifdef WIN32
      if (!_finite(trail->FuncValues[fNumber]))
#else
      if (!std::isfinite(trail->FuncValues[fNumber]))
#endif
      {
        //throw EXCEPTION("Infinite trail->FuncValues[fNumber]!");
        trail->index = -2;
        std::cout << " CalculateFuncs Error!!!\n";
      }
      else
        if ((fNumber == (_pTask->GetNumOfFunc() - 1)) || (trail->FuncValues[fNumber] > 0))
        {
          trail->index = fNumber;
        }
      fNumber++;
    }

    // Отправляем индекс точки
    MPI_Send(&(trail->index), 1, MPI_INT, 0, TagChildSolved, MPI_COMM_WORLD);
    // Отправляем точку
    MPI_Send(trail->y, parameters.Dimension, MPI_DOUBLE, 0, TagChildSolved, MPI_COMM_WORLD);
    // Отправляем обратно значение функции
    MPI_Send(trail->FuncValues, MaxNumOfFunc, MPI_DOUBLE, 0, TagChildSolved, MPI_COMM_WORLD);
  }
}


// ------------------------------------------------------------------------------------------------
int Solver::Solve()
{
  try
  {
    if (CheckParameters())
      return 1;
        
    if ((parameters.calculationsArray[0] == MPI_calc) && (parameters.GetProcNum() > 1) && (parameters.GetProcRank() > 0))
    {
      MpiCalculation();
    }
    else if ((parameters.TypeCalculation == AsyncMPI) && (parameters.GetProcNum() > 1) && (parameters.GetProcRank() > 0))
    {
      AsyncCalculation();
    }
    else
    {
      ClearData();
      CreateProcess();
      if (addPoints != nullptr)
        mProcess->InsertPoints(*addPoints);
      mProcess->Solve();
    }  

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
int Solver::Solve(Task* task)
{
  try
  {
    pTask = task;
    isExternalTask = true;

    CreateProcess();

    mProcess->Solve();
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
  return 0;
}

// ------------------------------------------------------------------------------------------------
Solver::~Solver()
{
  ClearData();
}

// ----------------------------------------------------------------------------
void Solver::InitAutoPrecision()
{
  // if user has not set the precision by the command line,
  // then set it to 1 / (2^((m + 1) * N) - 1)
  if (Extended::GetPrecision() == 0.01)
  {
    if (parameters.m * (parameters.Dimension - pTask->GetNumberOfDiscreteVariable()) <= 50)
    {
      Extended::SetTypeID(etDouble);
    }
    Extended::SetPrecision(1 / (::pow(2., parameters.m * (parameters.Dimension -
      pTask->GetNumberOfDiscreteVariable()))));
  }
}

// ------------------------------------------------------------------------------------------------
int Solver::CreateProcess()
{
  // В случае, если не совпадают (задача пришла снаружи — берём новые), иначе всё равно
  IProblem* _problem = mProblem;

  /// Создание задачи (Task) // перенести в фабрику
  if (pTask == 0)
  {
    pTask = TaskFactory::CreateTask(_problem, 0);
  }
  /// Создаём данные для поисковой информации




    if (pData == 0)
    {
      pData = new SearchData(_problem->GetNumberOfFunctions());
      int qSize = GLOBALIZER_MAX((int)pow(2.0, (int)(log((double)parameters.MaxNumOfPoints[pTask->GetProcLevel()])
        / log(2.0) - 2)) - 1, 1023);
      pData->ResizeQueue(qSize);
    }
    else
    {
      pData->Clear();
    }
  
  // Инициализируем числа с расширенной точностью
  InitAutoPrecision();

  if (mProcess != 0)
    //delete mProcess;
    mProcess->Reset(pData, pTask);
  else
    mProcess = new Process(*pData, *pTask);

  return 0;
}

// ------------------------------------------------------------------------------------------------
void Solver::SetProblem(IProblem* problem)
{
  mProblem = problem;
}

// ------------------------------------------------------------------------------------------------
IProblem* Solver::GetProblem()
{
  return mProblem;
}

// ------------------------------------------------------------------------------------------------
SolutionResult* Solver::GetSolutionResult()  /// best point
{
  if (result != 0)
    delete result;
  result = new SolutionResult();
  result->BestTrial = mProcess->GetOptimEstimation();
  result->IterationCount = mProcess->GetIterationCount();
  result->TrialCount = mProcess->GetNumberOfTrials();
  return result;
}

// ------------------------------------------------------------------------------------------------
void Solver::SetPoint(std::vector<Trial*>& points)
{
  addPoints = &points;

}

// ------------------------------------------------------------------------------------------------
std::vector<Trial*>& Solver::GetAllPoint()
{
  return pData->GetTrials();
}
