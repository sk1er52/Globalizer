/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      method_factory.cpp                                          //
//                                                                         //
//  Purpose:   Source file for method factory class                        //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "CalculationFactory.h"

#include "OmpCalculation.h"
#include "MPICalculationAsync.h"
#include "MPICalculation.h"
#include "CudaCalculation.h"


// ------------------------------------------------------------------------------------------------
Calculation* CalculationFactory::CreateCalculation2(Task& _pTask, Evolvent* evolvent)
{
  Calculation* calculation = 0;

  //if (parameters.calculationsArray.GetIsChange())
  //{
  //  if (parameters.calculationsArray[parameters.GetProcRank()] != -1)
  //  {
  //    printf("\nprocRank = %d\n", parameters.GetProcRank());
  //    char buf[256] = { 0 };
  //    sprintf(buf, "%d", parameters.calculationsArray[parameters.GetProcRank()]);
  //    parameters.SetVal("TypeCalculation", buf);
  //    parameters.PrintParameter("TypeCalculation");

  //    if (parameters.TypeCalculation == OMP)
  //      calculation = new OMPCalculation(_pTask);
  //    //else if (parameters.TypeCalculation == CUDA)
  //    //  calculation = new CUDACalculation(_parameters, _pTask);
  //    //else if (parameters.TypeCalculation == PHI)
  //    //  calculation = new TPHICalculation(_parameters, _pTask);
  //    if (parameters.TypeCalculation == BlockScheme)
  //      calculation = new TBlockSchemeCalculation(_pTask);
  //  }
  //}
  if (_pTask.IsLeaf()) //Если в листе
  {
    //else if (parameters.TypeCalculation == PHI)
    //  calculation = new TPHICalculation(_parameters, _pTask);

        // Выбор между OMP иCUDA
    if ((parameters.TypeCalculation == OMP) || (parameters.TypeCalculation == MPI_calc && parameters.GetProcNum() == 1))
    {
      if (Calculation::leafCalculation == 0)
      {
        calculation = new OMPCalculation(_pTask);
        Calculation::leafCalculation = calculation;
      }
      else
      {
        calculation = Calculation::leafCalculation;
      }
    }
    else if (parameters.TypeCalculation == CUDA)
    {
      if (Calculation::leafCalculation == 0)
      {
        calculation = new CUDACalculation(_pTask);
        Calculation::leafCalculation = calculation;
      }
      else
      {
        calculation = Calculation::leafCalculation;
      }
    }

    else if (parameters.TypeCalculation == MPI_calc)
    {
      if (Calculation::leafCalculation == 0)
      {
        calculation = new MPICalculation(_pTask);
        Calculation::leafCalculation = calculation;
      }
      else
      {
        calculation = Calculation::leafCalculation;
      }
    }
    else if (parameters.TypeCalculation == AsyncMPI)
    {
      if (Calculation::leafCalculation == 0)
      {
        calculation = new MPICalculationAsync(_pTask);
        Calculation::leafCalculation = calculation;
      }
      else
      {
        calculation = Calculation::leafCalculation;
      }
    }
    else
    {
      if (Calculation::leafCalculation == 0)
      {
        if (parameters.calculationsArray.GetSize() < _pTask.GetProcLevel())
          calculation = new CUDACalculation(_pTask);
        else
        {
          if (parameters.calculationsArray[_pTask.GetProcLevel()] == OMP)
            calculation = new OMPCalculation(_pTask);
          else
            calculation = new CUDACalculation(_pTask);

          Calculation::leafCalculation = calculation;
        }
      }
      else
      {
        calculation = Calculation::leafCalculation;
      }
    }
    //calculation = new OMPCalculation(_pTask);

    if (_pTask.GetProcLevel() == 0) //и одновременно в корне (если только корень)
    {

    }
  }


  return calculation;
}


// ------------------------------------------------------------------------------------------------
Calculation* CalculationFactory::CreateCalculation(Task& _pTask, Evolvent* evolvent)
{
  Calculation* calculation = 0;

  if (!_pTask.IsLeaf()) //Если в листе
  {    
    // Выбор между OMP иCUDA
    if ((parameters.TypeCalculation == OMP))
    {
      if (Calculation::leafCalculation == 0)
      {
        calculation = new OMPCalculation(_pTask);
        Calculation::leafCalculation = calculation;
      }
      else
      {
        calculation = Calculation::leafCalculation;
      }
    }

    if (parameters.TypeCalculation == CUDA)
    {
      if (Calculation::leafCalculation == 0)
      {
        calculation = new CUDACalculation(_pTask);
        Calculation::leafCalculation = calculation;
      }
      else
      {
        calculation = Calculation::leafCalculation;
      }
    }
    else if (parameters.TypeCalculation == MPI_calc)
    {
      if (Calculation::leafCalculation == 0)
      {
        calculation = new MPICalculation(_pTask);
        Calculation::leafCalculation = calculation;
      }
      else
      {
        calculation = Calculation::leafCalculation;
      }
    }
    else if (parameters.TypeCalculation == AsyncMPI)
    {
      if (Calculation::leafCalculation == 0)
      {
        calculation = new MPICalculationAsync(_pTask);
        Calculation::leafCalculation = calculation;
      }
      else
      {
        calculation = Calculation::leafCalculation;
      }
    }
    else
    {
      if (Calculation::leafCalculation == 0)
      {
        if (parameters.calculationsArray.GetSize() < _pTask.GetProcLevel())
          calculation = new CUDACalculation(_pTask);
        else
        {          
          if (parameters.calculationsArray[_pTask.GetProcLevel()] == OMP)
            calculation = new OMPCalculation(_pTask);
          else
            calculation = new CUDACalculation(_pTask);

          Calculation::leafCalculation = calculation;
        }
      }
      else
      {
        calculation = Calculation::leafCalculation;
      }
    }

  }


  return calculation;
}

Calculation* CalculationFactory::CreateNewCalculation(Task& _pTask, Evolvent* evolvent)
{
  Calculation* calculation = 0;

  if (_pTask.IsLeaf()) // Если в листе
  {    

    // Выбор между OMP и CUDA
    if ((parameters.TypeCalculation == OMP))
    {
      calculation = new OMPCalculation(_pTask);
    }
    else if (parameters.TypeCalculation == CUDA)
    {      
      calculation = new CUDACalculation(_pTask);        
    }
    else if (parameters.TypeCalculation == MPI_calc)
    {
      calculation = new MPICalculation(_pTask);      
    }
    else if (parameters.TypeCalculation == AsyncMPI)
    {      
      calculation = new MPICalculationAsync(_pTask);        
    }
    else
    {
        if (parameters.calculationsArray.GetSize() < _pTask.GetProcLevel())
          calculation = new CUDACalculation(_pTask);
        else
        {          
          if (parameters.calculationsArray[_pTask.GetProcLevel()] == OMP)
            calculation = new OMPCalculation(_pTask);
          else
            calculation = new CUDACalculation(_pTask);

          Calculation::leafCalculation = calculation;
        }
     
    }


  }
  return calculation;
}
