/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2016 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      problemWithConstraints.h                                    //
//                                                                         //
//  Purpose:   Header file for Globalizer problem interface                //
//                                                                         //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

/**
\file problemWithConstraints.h

\authors Ëåáåäåâ È.
\date 2016
\copyright ÍÍÃÓ èì. Í.È. Ëîáà÷åâñêîãî

\brief Îáúÿâëåíèå áàçîâîãî êëàññà äëÿ çàäà÷ ñ îãðàíè÷åíèÿìè

\details
*/

#ifndef __PROBLEM_WITH_CONSTRAINTS_H__
#define __PROBLEM_WITH_CONSTRAINTS_H__

#include <omp.h>
#include "Problem.h"
#include "ConstrainedProblem.hpp"
#include "ConstrainedProblemGenerator.hpp"


template <class Owner, class FType>
class ProblemWithConstraints :
  public BaseProblem<Owner>, public ConstrainedProblemGenerator<FType>, public ConstrainedProblem<FType>
{
#undef OWNER_NAME
#define OWNER_NAME ProblemWithConstraints
protected:
  /// Ïðîâåðêà ïðàâèëüíîñòè ïîñëå îêîí÷àíèÿ ÷òåíèÿ ïàðàìåòðîâ
  virtual int CheckValue(int index = -1);

  /// Çàäàíèå çíà÷åíèé ïî óìîë÷àíèþ áàçîâûõ ïàðàìåòðîâ
  virtual void SetBaseDefaultParameters();


  /** Ìåòîä âîçâðàùàåò ÷èñëî îáùåå ôóíêöèé â çàäà÷å (îíî ðàâíî ÷èñëî îãðàíè÷åíèé + ÷èñëî êðèòåðèåâ)
  \return ×èñëî ôóíêöèé
  */
  virtual int GetRealNumberOfFunctions() const;

  /** Ìåòîä âîçâðàùàåò ÷èñëî îãðàíè÷åíèé â çàäà÷å
  \return ×èñëî îãðàíè÷åíèé
  */
  virtual int GetRealNumberOfConstraints() const;

  template <class ClassType, class Type>
  void ProblemParameterCheckSize(ClassType& par, Type defVal, Type left, Type right)
  {
    if (par.GetSize() < GetRealNumberOfFunctions())
    {
      if (par.GetSize() == 1)
      {
        Type d = par[0];
        if ((d > right) || (d < left))
          d = defVal;
        par.SetSize(GetRealNumberOfFunctions());
        for (int i = 0; i < GetRealNumberOfFunctions(); i++)
          par[i] = d;
      }
      else
      {
        int oldSize = par.GetSize();
        Type d = par[0];
        if ((d > right) || (d < left))
          d = defVal;
        par.SetSize(GetRealNumberOfFunctions());
        for (int i = oldSize; i < GetRealNumberOfFunctions(); i++)
          (par.GetData())[i] = d;
      }
    }
  }
public:

  ProblemWithConstraints();

  /// Êîëè÷åñòâî îãðàíè÷åíèé
  TInt<Owner> constraint_count;
  /// Äîëÿ îáëàñòè ïîèñêà, åñëè çàäîíî, òî Q âû÷èñëÿåòñÿ
  TDoubles<Owner> Deltas;
  /// Ñäâèã îãðàíè÷åíèé, îíî æå RHS
  TDoubles<Owner> Q;
  /// Ìàñøòàáèðîâàòü èëè íåò îãðàíè÷åíèÿ, åñëè õîòü ó îäíîãî îãðàíè÷åíèÿ ãëîáàëüíûé ìèíèìóì íà ãðàíèöå, òî íå èñïîëüçóåòñÿ
  TFlag<Owner> IsZoom;
  /// Ñäâèãàòü èëè íåò ãëîáàëüíûé ìèíèìóì îãðàíè÷åíèé â êîîðäèíàòû ãëîáàëüíîãî ìèíèìóìà öåëåâîé ôóíêöèè
  TFlag<Owner> IsShift;
  /// Ñäâèãàòü èëè íåò ãëîáàëüíûé ìèíèìóì öåëåâîé ôóíêöèè íà ãðàíèöó äîïóñòèìîé îáëàñòè
  TFlag<Owner> IsBoundaryShift;
  /// Òî÷íîñòü ïîèñêà áëèæàéøåé òî÷êè íà ãðàíèöå äîïóñòèìîé îáëàñòè (ñòåïåíü 2^{-1})
  TInt<Owner> BoundarySearchPrecision;
  /// èçìåíÿòü ëè öåëåâóþ ôóíêöèþ ïóòåì ïðèáàëåíèÿ ôóíêöèîíàëà îò îãðàíè÷åíèé
  TFlag<Owner> IsImprovementOfTheObjective;
  /// Êîýôèöèåíòû èçìåíåíèÿ
  TDoubles<Owner> ImprovementCoefficients;

  /// Èíèöèàëèçèðóåò ôóíêöèþ ñ íîìåðîì index
  //virtual void InitFunc(FType* func, int index) = 0;

  /// Èíèöèàëèçàöèÿ ïàðàìåòðîâ
  virtual void Init(int argc, char* argv[], bool isMPIInit = false);

  virtual double CalculateFunctionals(const double* y, int fNumber);


  /** Ìåòîä âîçâðàùàåò ÷èñëî îáùåå ôóíêöèé â çàäà÷å (îíî ðàâíî ÷èñëî îãðàíè÷åíèé + ÷èñëî êðèòåðèåâ)
  \return ×èñëî ôóíêöèé
  */
  virtual int GetNumberOfFunctions() const;
  /** Ìåòîä âîçâðàùàåò ÷èñëî îãðàíè÷åíèé â çàäà÷å
  \return ×èñëî îãðàíè÷åíèé
  */
  virtual int GetNumberOfConstraints() const;

  /// Âîçâðàùàåò ÷èñëî îãðàíè÷åíèé
  virtual int GetConstraintsNumber() const
  {
    return GetNumberOfConstraints();
  }

  /// Âîçâðàùàåò, íóæíî ëè ìàñøòàáèðîâàòü çàäà÷ó
  bool GetIsZoom() const
  {
    return this->mIsZoom;
  }
  /// Çàäàåò, íóæíî ëè ìàñøòàáèðîâàòü çàäà÷ó
  void SetIsZoom(bool isZoom)
  {
    this->mIsZoom = isZoom;
  }
  /// Âîçâðàùàåò, íóæíî ëè ñäâèãàòü ôóíêöèè îãðàíè÷åíèé
  bool GetIsShift() const
  {
    return this->mIsShift;
  }
  /// Çàäàåò, íóæíî ëè ñäâèãàòü ôóíêöèè îãðàíè÷åíèé
  void SetIsShift(bool isShift)
  {
    this->mIsShift = isShift;
  }
  /// Âîçâðàùàåò, íóæíî ëè ñäâèãàòü îïòèìóì öåëåâîé ôóíêöèè íà ãðàíèöó
  bool GetIsBoundaryShift() const
  {
    return this->mIsBoundaryShift;
  }
  /// Çàäàåò, íóæíî ëè ñäâèãàòü îïòèìóì öåëåâîé ôóíêöèè íà ãðàíèöó
  void SetIsBoundaryShift(bool isBoundaryShift)
  {
    this->mIsBoundaryShift = isBoundaryShift;
  }
  /// Âîçâðàùàåò òî÷íîñòü ïîèñêà áëèæàéøåé òî÷êè íà ãðàíèöå îáëàñòè (ñòåïåíü 0.5)
  int GetBoundarySearchPrecision() const
  {
    return this->mBoundarySearchPrecision;
  }
  /// Çàäàåò òî÷íîñòü ïîèñêà áëèæàéøåé òî÷êè íà ãðàíèöå îáëàñòè (ñòåïåíü 0.5)
  void SetBoundarySearchPrecision(int boundarySearchPrecision)
  {
    this->mBoundarySearchPrecision = boundarySearchPrecision;
  }
  /// Âîçâðàùàåò, íóæíî ëè ìîäèôèöèðîâàòü öåëåâóþ ôóíêöèþ
  bool GetIsImprovementOfTheObjective() const
  {
    return this->mIsImprovementOfTheObjective;
  }
  /// Çàäàåò, íóæíî ëè ìîäèôèöèðîâàòü öåëåâóþ ôóíêöèþ
  void SetIsImprovementOfTheObjective(bool isImprovementOfTheObjective)
  {
    this->mIsImprovementOfTheObjective = isImprovementOfTheObjective;
  }

  /** Ìåòîä âîçâðàùàåò êîîðäèíàòû òî÷êè ãëîáàëüíîãî ìèíèìóìà öåëåâîé ôóíêöèè
  \param[out] y òî÷êà, â êîòîðîé äîñòèãàåòñÿ îïòèìàëüíîå çíà÷åíèå
  \return Êîä îøèáêè (#OK èëè #UNDEFINED)
  */
  virtual int GetOptimumPoint(double* y) const
  {
    return BaseProblem<Owner>::UNDEFINED;
  }

  /// Âîçâðàùàåò òî÷êó ãëîáàëüíîãî îïòèìóìà äëÿ ôóíêöèè fNumber
  virtual int GetConstraintOptimumPoint(double* point, int fNumber)
  {
    return BaseProblem<Owner>::UNDEFINED;
  }
};

// ------------------------------------------------------------------------------------------------
/// Ïðîâåðêà ïðàâèëüíîñòè ïîñëå îêîí÷àíèÿ ÷òåíèÿ ïàðàìåòðîâ
template <class Owner, class FType>
int ProblemWithConstraints<Owner, FType>::CheckValue(int index)
{
  BaseProblem<Owner>::CheckValue(index);

  if ((index < this->mOptionsCount) && (index >= 0))
  {
    if (this->mOptions[index] == &Deltas)
    {
      if (Deltas.GetSize() == 1)
        this->mIsTotalDelta = true;
      else
        this->mIsTotalDelta = false;
    }
  }

  if (this->mIsInit)
  {
    if ((constraint_count > 50) || (constraint_count < 0))
      constraint_count = 0;

    if (Deltas.GetSize() < constraint_count)
    {
      if (Deltas.GetSize() == 1)
      {
        double d = Deltas[0];
        if ((d > 1) || (d < 0))
          d = 0.5;
        Deltas.SetSize(constraint_count);
        for (int i = 0; i < constraint_count; i++)
          Deltas[i] = d;
      }
      else
      {
        int oldSize = Deltas.GetSize();
        double d = Deltas[0];
        if ((d > 1) || (d < 0))
          d = 0.5;
        Deltas.SetSize(constraint_count);
        for (int i = oldSize; i < constraint_count; i++)
          Deltas[i] = d;
      }
    }

    if (Q.GetSize() < GetRealNumberOfFunctions())
    {
      if (Q.GetSize() == 1)
      {
        double d = Q[0];
        if ((d > 1) || (d < 0))
          d = 0.5;
        Q.SetSize(GetRealNumberOfFunctions());
        for (int i = 0; i < GetRealNumberOfConstraints(); i++)
          Q[i] = d;
        for (int i = GetRealNumberOfConstraints(); i < GetRealNumberOfFunctions(); i++)
          Q[i] = 0;
      }
      else
      {
        int oldSize = Q.GetSize();
        double d = Q[0];
        if ((d > 1) || (d < 0))
          d = 0.5;
        Q.SetSize(GetRealNumberOfFunctions());
        for (int i = oldSize; i < GetRealNumberOfConstraints(); i++)
          Q[i] = d;
        for (int i = std::max(GetRealNumberOfConstraints(), oldSize); i < GetRealNumberOfFunctions(); i++)
          Q[i] = 0;
      }
    }


    if (ImprovementCoefficients.GetSize() < GetRealNumberOfFunctions())
    {
      if (ImprovementCoefficients.GetSize() == 1)
      {
        double d = ImprovementCoefficients[0];
        ImprovementCoefficients.SetSize(GetRealNumberOfFunctions());
        for (int i = 0; i < GetRealNumberOfConstraints(); i++)
          ImprovementCoefficients[i] = d;
        for (int i = GetRealNumberOfConstraints(); i < GetRealNumberOfFunctions(); i++)
          ImprovementCoefficients[i] = 0;
      }
      else
      {
        int oldSize = ImprovementCoefficients.GetSize();
        double d = ImprovementCoefficients[0];
        ImprovementCoefficients.SetSize(GetRealNumberOfFunctions());
        for (int i = oldSize; i < GetRealNumberOfConstraints(); i++)
          ImprovementCoefficients[i] = d;
        for (int i = std::max(GetRealNumberOfConstraints(), oldSize); i < GetRealNumberOfFunctions(); i++)
          ImprovementCoefficients[i] = 0;
      }
    }

  }
  return 0;
}

// ------------------------------------------------------------------------------------------------
/// Çàäàíèå çíà÷åíèé ïî óìîë÷àíèþ áàçîâûõ ïàðàìåòðîâ
template <class Owner, class FType>
void ProblemWithConstraints<Owner, FType>::SetBaseDefaultParameters()
{
  BaseProblem<Owner>::SetBaseDefaultParameters();
  int sizeVal = 1;

  Owner* aqqq = (Owner*)this;


  constraint_count.InitializationParameterProperty(
    aqqq, &ProblemWithConstraints::CheckValue, this->mOptionsCount, this->Separator, sizeVal,
    "constraint_count", "Constraint count", "-CC", "0");
  this->AddOption((BaseProperty<Owner>*)(&constraint_count));

  Deltas.InitializationParameterProperty(
    aqqq, &ProblemWithConstraints::CheckValue, this->mOptionsCount, this->Separator, sizeVal,
    "delta", "The share of the search area", "-delta", "0.5");
  this->AddOption((BaseProperty<Owner>*)(&Deltas));
  this->mIsTotalDelta = true;

  Q.InitializationParameterProperty(
    aqqq, &ProblemWithConstraints::CheckValue, this->mOptionsCount, this->Separator, sizeVal,
    "Q", "Shift constraints", "-Q", "0");
  this->AddOption((BaseProperty<Owner>*)(&Q));

  this->mQ = (double**)Q.GetValue();

  IsZoom.InitializationParameterProperty(
    aqqq, &ProblemWithConstraints::GetIsZoom, &ProblemWithConstraints::SetIsZoom, &ProblemWithConstraints::CheckValue,
    this->mOptionsCount, this->Separator, sizeVal,
    "IsZoom", "Is Zoom", "-IZ", "false");
  this->AddOption((BaseProperty<Owner>*)(&IsZoom));

  IsShift.InitializationParameterProperty(
    aqqq, &ProblemWithConstraints::GetIsShift, &ProblemWithConstraints::SetIsShift, &ProblemWithConstraints::CheckValue,
    this->mOptionsCount, this->Separator, sizeVal,
    "IsShift", "Is Shift", "-ISH", "false");
  this->AddOption((BaseProperty<Owner>*)(&IsShift));

  IsBoundaryShift.InitializationParameterProperty(
    aqqq, &ProblemWithConstraints::GetIsBoundaryShift, &ProblemWithConstraints::SetIsBoundaryShift, &ProblemWithConstraints::CheckValue,
    this->mOptionsCount, this->Separator, sizeVal,
    "IsBoundaryShift", "Is Boundary Shift", "-IBSH", "false");
  this->AddOption((BaseProperty<Owner>*)(&IsBoundaryShift));

  BoundarySearchPrecision.InitializationParameterProperty(
    aqqq, &ProblemWithConstraints::GetBoundarySearchPrecision, &ProblemWithConstraints::SetBoundarySearchPrecision, &ProblemWithConstraints::CheckValue,
    this->mOptionsCount, this->Separator, sizeVal,
    "BoundarySearchPrecision", "Boundary Search Precision", "-bsp", "20");
  this->AddOption((BaseProperty<Owner>*)(&BoundarySearchPrecision));

  IsImprovementOfTheObjective.InitializationParameterProperty(aqqq,
    &ProblemWithConstraints::GetIsImprovementOfTheObjective, &ProblemWithConstraints::SetIsImprovementOfTheObjective, &ProblemWithConstraints::CheckValue,
    this->mOptionsCount, this->Separator, sizeVal,
    "IsImprovementOfTheObjective", "Is improvement of the objective", "-IIO", "false");
  this->AddOption((BaseProperty<Owner>*)(&IsImprovementOfTheObjective));

  ImprovementCoefficients.InitializationParameterProperty(
    aqqq, &ProblemWithConstraints::CheckValue, this->mOptionsCount, this->Separator, sizeVal,
    "ImprovementCoefficients", "Improvement of the objective coefficients", "-IC", "100_100");
  this->AddOption((BaseProperty<Owner>*)(&ImprovementCoefficients));

  this->mImprovementCoefficients = (double**)ImprovementCoefficients.GetValue();
}

// ------------------------------------------------------------------------------------------------
/** Ìåòîä âîçâðàùàåò ÷èñëî îáùåå ôóíêöèé â çàäà÷å (îíî ðàâíî ÷èñëî îãðàíè÷åíèé + ÷èñëî êðèòåðèåâ)
\return ×èñëî ôóíêöèé
*/
template <class Owner, class FType>
int ProblemWithConstraints<Owner, FType>::GetRealNumberOfFunctions() const
{
  return this->mNumberOfCriterions + constraint_count.GetAvailableData();
}

// ------------------------------------------------------------------------------------------------
/** Ìåòîä âîçâðàùàåò ÷èñëî îãðàíè÷åíèé â çàäà÷å
\return ×èñëî îãðàíè÷åíèé
*/
template <class Owner, class FType>
int ProblemWithConstraints<Owner, FType>::GetRealNumberOfConstraints() const
{
  return constraint_count.GetAvailableData();
}

// ------------------------------------------------------------------------------------------------
template <class Owner, class FType>
ProblemWithConstraints<Owner, FType>::ProblemWithConstraints() :
  BaseProblem<Owner>(), ConstrainedProblemGenerator<FType>(), ConstrainedProblem<FType>()
{
  this->mPFunction = 0;
}

// ------------------------------------------------------------------------------------------------
/// Èíèöèàëèçàöèÿ ïàðàìåòðîâ
template <class Owner, class FType>
void ProblemWithConstraints<Owner, FType>::Init(int argc, char* argv[], bool isMPIInit)
{
  BaseProblem<Owner>::Init(argc, argv, false);
  this->mPFunction = new FType*[GetRealNumberOfFunctions()];
  this->mZoomRatios = new double[GetRealNumberOfFunctions()];
  this->mShift = new double*[GetRealNumberOfFunctions()];
  this->mBoundaryShift = new double[this->Dimension];
  this->mvZoomRatios = &this->mZoomRatios;
  this->mvShift = &this->mShift;
  this->mvBoundaryShift = this->mBoundaryShift;

  for (int i = 0; i < GetRealNumberOfFunctions(); i++)
  {
    this->mZoomRatios[i] = 1.0;

    this->mShift[i] = new double[this->Dimension];
    for (int j = 0; j < this->Dimension; j++)
    {
      this->mShift[i][j] = 0.0;
    }
  }
  for (int i = 0; i < this->Dimension; i++)
  {
    this->mBoundaryShift[i] = 0.0;
  }
  for (int i = 0; i < GetRealNumberOfFunctions(); i++)
  {
    this->mPFunction[i] = new FType();
    this->mPConstraints.push_back(this->mPFunction[i]);
    this->InitFunc(this->mPFunction[i], i);
  }
  this->mPObjective = this->mPFunction[GetRealNumberOfFunctions() - 1];

  if (IsZoom)
  {
    this->SetZoom();
  }

  if (IsShift)
  {
    this->SetShift();
  }

  this->mmIsImprovementOfTheObjective = IsImprovementOfTheObjective;
  this->mmQ = new double*[1];
  (*this->mmQ) = new double[GetRealNumberOfConstraints() + 1];

  this->improvementCoefficients = new double*[1];
  *this->improvementCoefficients = new double[GetRealNumberOfConstraints() + 1];
  for (int i = 0; i < GetRealNumberOfConstraints(); i++)
    (*this->improvementCoefficients)[i] = ImprovementCoefficients[i];
  (*this->improvementCoefficients)[GetRealNumberOfConstraints()] = 0;

  this->currentFunctionNumber = 0;
  if (Deltas.GetIsReadValue())
  {
    /*for (int j = 0; j < Dimension; j++)
      objectiveMin[j] = 0;*/
      //GetOptimumPoint(objectiveMin);
    if (this->mIsTotalDelta)
    {
      double q = this->CalculateRHS(Deltas[0]);
      for (int i = 0; i < GetRealNumberOfConstraints(); i++)
      {
        Q[i] = q;
        (*this->mmQ)[i] = Q[i];
      }
    }
    else
    {
      for (int i = 0; i < GetRealNumberOfConstraints(); i++)
      {
        this->currentFunctionNumber = i;
        this->currentFunction = this->mPFunction[i];
        Q[i] = this->CalculateRHS(Deltas[i]);
        (*this->mmQ)[i] = Q[i];
      }
    }
  }

  if (IsBoundaryShift)
  {
    this->SetBoundaryShift();
  }

}

// ------------------------------------------------------------------------------------------------
template <class Owner, class FType>
double ProblemWithConstraints<Owner, FType>::CalculateFunctionals(const double* y, int fNumber)
{

  double koef = 1.0;
  double sum = 0.0;
  //for (int j = 0; j < this->Dimension; j++)
  //{
  //  for (int k = 0; k < 1000; k++)
  //    koef = (koef + exp(1.0) * y[j] / 4.0) * (y[j] / 4.0);
  //}
  sum += koef * this->CalculateFunction(y, fNumber) / koef;
  return sum;
}

// ------------------------------------------------------------------------------------------------
/** Ìåòîä âîçâðàùàåò ÷èñëî îáùåå ôóíêöèé â çàäà÷å (îíî ðàâíî ÷èñëî îãðàíè÷åíèé + ÷èñëî êðèòåðèåâ)
\return ×èñëî ôóíêöèé
*/
template <class Owner, class FType>
int ProblemWithConstraints<Owner, FType>::GetNumberOfFunctions() const
{
  return this->mNumberOfCriterions + constraint_count.GetAvailableData();
}

// ------------------------------------------------------------------------------------------------
/** Ìåòîä âîçâðàùàåò ÷èñëî îãðàíè÷åíèé â çàäà÷å
\return ×èñëî îãðàíè÷åíèé
*/
template <class Owner, class FType>
int ProblemWithConstraints<Owner, FType>::GetNumberOfConstraints() const
{
  return constraint_count.GetAvailableData();
}

#endif /// __PROBLEM_WITH_CONSTRAINTS_H__
// - end of file ----------------------------------------------------------------------------------
