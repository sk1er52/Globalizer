#include "Trial.h"

#include "SearchInterval.h"

Trial* Trial::GetLeftPoint()
{
  if (this->leftInterval != 0)
  {
    if (this->leftInterval->LeftPoint != 0)
    {
      return this->leftInterval->LeftPoint;
    }
  }
  return 0;
}

Trial* Trial::GetRightPoint()
{
  if (this->rightInterval != 0)
  {
    if (this->rightInterval->RightPoint != 0)
    {
      return this->rightInterval->RightPoint;
    }
  }
  return 0;
}