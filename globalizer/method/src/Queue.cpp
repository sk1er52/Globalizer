/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      Queue.cpp                                                   //
//                                                                         //
//  Purpose:   Source file for priority queue class                        //
//                                                                         //
//  Author(s): Sysoyev A., Barkalov K., Sovrasov V                         //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "Queue.h"
#include "BaseInterval.h"
// ------------------------------------------------------------------------------------------------
PriorityQueue::PriorityQueue(int _MaxSize)
{
  int tmpPow2 = _MaxSize + 1;
  if (tmpPow2&(tmpPow2 - 1))
  {
    throw EXCEPTION("Max size of queue not divisible by power of two");
  }
  MaxSize = _MaxSize;
  CurSize = 0;
  pMem = new MinMaxHeap< QueueElement, _less >(MaxSize);
}

// ------------------------------------------------------------------------------------------------
PriorityQueue::~PriorityQueue()
{
  delete pMem;
}

// ------------------------------------------------------------------------------------------------
bool PriorityQueue::IsEmpty() const
{
  return CurSize == 0;
}
// ------------------------------------------------------------------------------------------------
int PriorityQueue::GetSize() const
{
  return CurSize;
}
// ------------------------------------------------------------------------------------------------
int PriorityQueue::GetMaxSize() const
{
  return MaxSize;
}

// ------------------------------------------------------------------------------------------------
bool PriorityQueue::IsFull() const
{
  return CurSize == MaxSize;
}

// ------------------------------------------------------------------------------------------------
QueueElement* PriorityQueue::Push(double globalKey, double localKey, void *value)
{
  QueueElement* a = 0;
  if (!IsFull()) {
    CurSize++;
    a = pMem->push(QueueElement(globalKey, value));
  }
  else {
    if (globalKey > pMem->findMin().Key)
      DeleteMinElem();
    else
      return a;
    CurSize++;
    a = pMem->push(QueueElement(globalKey, value));
  }
  return a;
}

// ------------------------------------------------------------------------------------------------
QueueElement* PriorityQueue::PushWithPriority(double globalKey, double localKey, void *value)
{
  QueueElement* a = 0;
  if (!IsEmpty()) {
    if (globalKey >= pMem->findMin().Key) {
      if (IsFull())
        DeleteMinElem();
      CurSize++;
      a = pMem->push(QueueElement(globalKey, value));
    }
  }
  else {
    CurSize++;
    a = pMem->push(QueueElement(globalKey, value));
  }
  return a;
}

// ------------------------------------------------------------------------------------------------
void PriorityQueue::Pop(double *key, void **value)
{
  if (CurSize != 0)
  {
    QueueElement tmp = pMem->popMax();
    *key = tmp.Key;
    *value = tmp.pValue;
    CurSize--;
  }
  else
  {
    throw EXCEPTION("Cannot pop element from empty queue");
  }
}

// ------------------------------------------------------------------------------------------------
void PriorityQueue::DeleteByValue(void *value)
{
  QueueElement* heapMem = pMem->getHeapMemPtr();
  for (int i = 0; i < CurSize; i++)
  {
    if (heapMem[i].pValue == value)
    {
      pMem->deleteElement(heapMem + i);
      CurSize--;

      break;
    }
  }
}

// ------------------------------------------------------------------------------------------------
void PriorityQueue::DeleteElement(QueueElement* item)
{
  pMem->deleteElement(item);
  CurSize--;
}

// ------------------------------------------------------------------------------------------------
void PriorityQueue::TrickleUp(QueueElement * item)
{
  pMem->TrickleUp(item);
}

// ------------------------------------------------------------------------------------------------
void PriorityQueue::DeleteMinElem()
{
  QueueElement tmp = pMem->popMin();
  CurSize--;
}

// ------------------------------------------------------------------------------------------------
void PriorityQueue::Clear()
{
  pMem->clear();
  CurSize = 0;
}

// ------------------------------------------------------------------------------------------------
void PriorityQueue::Resize(int size)
{
  CurSize = 0;
  MaxSize = size;
  delete pMem;
  pMem = new MinMaxHeap< QueueElement, _less >(MaxSize);
}

// ------------------------------------------------------------------------------------------------
QueueElement& PriorityQueue::FindMax()
{
  return pMem->findMax();
}
// - end of file ----------------------------------------------------------------------------------
