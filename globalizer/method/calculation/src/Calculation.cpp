#include "Calculation.h"

int Calculation::countCalculation = 0;

Calculation* Calculation::leafCalculation = 0;
Calculation* Calculation::firstCalculation = 0;

bool Calculation::isStartComputingAway = true;
TResultForCalculation Calculation::resultCalculation;
InformationForCalculation Calculation::inputCalculation;


Calculation::Calculation(Task& _pTask) : pTask(&_pTask)
{  }

void Calculation::SetCountCalculation(int c)
{
  countCalculation = c;
  isStartComputingAway = false;
}


void Calculation::SetTask(Task* _pTask)
{
  pTask = _pTask;
}

void Calculation::SetSearchData(TSearchData* _pData)
{
  pData = _pData;
}

void Calculation::Reset()
{
}
