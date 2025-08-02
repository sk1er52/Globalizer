/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2016 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      problemWithMulctFunction.h                                  //
//                                                                         //
//  Purpose:   Header file for Globalizer problem interface                //
//                                                                         //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

/**
\file problemWithMulctFunction.h

\authors Ëåáåäåâ Â.
\date 2016
\copyright ÍÍÃÓ èì. Í.È. Ëîáà÷åâñêîãî

\brief Îáúÿâëåíèå áàçîâîãî êëàññà äëÿ çàäà÷ ñ îãðàíè÷åíèÿìè

\details
*/

#ifndef __PROBLEM_WITH_MULCT_FUNCTION_H__
#define __PROBLEM_WITH_MULCT_FUNCTION_H__

#include <omp.h>
#include "problemWithConstraints.h"

template <class Owner, class Func>
class ProblemWithMulctFunction: public ProblemWithConstraints<Owner, Func>
{
  #undef OWNER_NAME
#define OWNER_NAME ProblemWithMulctFunction
protected:

  virtual int CheckValue(int index = -1)
  {
    ProblemWithConstraints::CheckValue(index);

    if (mIsInit)
    {
      if (Mulcts.GetSize() < constraint_count)
      {
        if (Mulcts.GetSize() == 1)
        {
          double d = Mulcts[0];
          Mulcts.SetSize(constraint_count);
          for (int i = 0; i < constraint_count; i++)
            Mulcts[i] = d;
        }
        else
        {
          int oldSize = Mulcts.GetSize();
          double d = Mulcts[0];
          Mulcts.SetSize(constraint_count);
          for (int i = oldSize; i < constraint_count; i++)
            Mulcts[i] = d;
        }
      }
    }

    return 0;
  }

  /// Çàäàíèå çíà÷åíèé ïî óìîë÷àíèþ áàçîâûõ ïàðàìåòðîâ
  virtual void SetBaseDefaultParameters()
  {
    ProblemWithConstraints::SetBaseDefaultParameters();
    int sizeVal = 1;

    Owner* aqqq = (Owner*)this;

    Mulcts.InitializationParameterProperty(aqqq, &ProblemWithMulctFunction::CheckValue,
      mOptionsCount, Separator, sizeVal, "mulct", "Mulcts", "-Mulcts", "10000");
    AddOption((BaseProperty<Owner>*)(&Mulcts));
  }

public:

  ProblemWithMulctFunction() : ProblemWithConstraints()
  {
  }

  /// Êîýôèöèåíò øòðàôà
  TDoubles<Owner> Mulcts;

  virtual double CalculateFunctionals(const double* y, int fNumber)
  {
    double x[MAX_TRIAL_DIMENSION];
    //for (int i = 0; i < Dimension; i++)
    //{
    //  x[i] = (y[i] - mShift[GetRealNumberOfConstraints()][i]) / mZoomRatios[GetRealNumberOfConstraints()];
    //}
    double result = ProblemWithConstraints::CalculateFunctionals(y, GetRealNumberOfConstraints()); //mPFunction[GetRealNumberOfConstraints()]->Calculate(y);
    for (int i = 0; i < GetRealNumberOfConstraints(); i++)
    {
      for (int j = 0; j < Dimension; j++)
      {
        x[j] = (y[j] - mShift[i][j]) / mZoomRatios[i];
      }
      double v = GLOBALIZER_MAX(mPFunction[i]->Calculate(x) - Q[i], 0.0);
      result += Mulcts[i] * (v * v);
    }
    return result;
  }

    /** Ìåòîä âîçâðàùàåò ÷èñëî îáùåå ôóíêöèé â çàäà÷å (îíî ðàâíî ÷èñëî îãðàíè÷åíèé + ÷èñëî êðèòåðèåâ)
  \return ×èñëî ôóíêöèé
  */
  virtual int GetNumberOfFunctions() const
  {
    return 1;
  }

  /** Ìåòîä âîçâðàùàåò ÷èñëî îãðàíè÷åíèé â çàäà÷å
  \return ×èñëî îãðàíè÷åíèé
  */
  virtual int GetNumberOfConstraints() const
  {
    return 0;
  }

};


#endif /// __PROBLEM_WITH_MULCT_FUNCTION_H__
// - end of file ----------------------------------------------------------------------------------
