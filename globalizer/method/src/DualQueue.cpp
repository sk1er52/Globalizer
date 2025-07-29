/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      dual_queue.cpp                                              //
//                                                                         //
//  Purpose:   Source file for priority queue class                        //
//                                                                         //
//  Author(s): Sysoyev A., Barkalov K., Sovrasov V.                        //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "DualQueue.h"

// ------------------------------------------------------------------------------------------------
PriorityDualQueue::PriorityDualQueue(int _MaxSize)
{
  int tmpPow2 = _MaxSize + 1;
  if (tmpPow2&(tmpPow2 - 1))
  {
    throw EXCEPTION("Max size of queue not divisible by power of two");
  }
  MaxSize = _MaxSize;
  CurGlobalSize = CurLocalSize = 0;
  pLocalHeap = new MinMaxHeap< QueueElement, _less >(MaxSize);
  pGlobalHeap = new MinMaxHeap< QueueElement, _less >(MaxSize);
}

// ------------------------------------------------------------------------------------------------
PriorityDualQueue::~PriorityDualQueue()
{
  delete pLocalHeap;
  delete pGlobalHeap;
}

// ------------------------------------------------------------------------------------------------
void PriorityDualQueue::DeleteMinLocalElem()
{
  QueueElement tmp = pLocalHeap->popMin();
  CurLocalSize--;

  //update linked element in the global queue
  if (tmp.pLinkedElement != NULL)
    tmp.pLinkedElement->pLinkedElement = NULL;
}

// ------------------------------------------------------------------------------------------------
void PriorityDualQueue::DeleteMinGlobalElem()
{
  QueueElement tmp = pGlobalHeap->popMin();
  CurGlobalSize--;

  //update linked element in the local queue
  if (tmp.pLinkedElement != NULL)
    tmp.pLinkedElement->pLinkedElement = NULL;
}

// ------------------------------------------------------------------------------------------------
int PriorityDualQueue::GetLocalSize() const
{
  return CurLocalSize;
}

// ------------------------------------------------------------------------------------------------
int PriorityDualQueue::GetSize() const
{
  return CurGlobalSize;
}
// ------------------------------------------------------------------------------------------------
int PriorityDualQueue::GetMaxSize() const
{
  return MaxSize;
}

// ------------------------------------------------------------------------------------------------
bool PriorityDualQueue::IsLocalEmpty() const
{
  return CurLocalSize == 0;
}

// ------------------------------------------------------------------------------------------------
bool PriorityDualQueue::IsLocalFull() const
{
  return CurLocalSize == MaxSize;
}

// ------------------------------------------------------------------------------------------------
bool PriorityDualQueue::IsEmpty() const
{
  return CurGlobalSize == 0;
}

// ------------------------------------------------------------------------------------------------
bool PriorityDualQueue::IsFull() const
{
  return CurGlobalSize == MaxSize;
}

// ------------------------------------------------------------------------------------------------
QueueElement* PriorityDualQueue::Push(double globalKey, double localKey, void * value)
{
  QueueElement* pGlobalElem = NULL, *pLocalElem = NULL;
  //push to a global queue
  if (!IsFull()) {
    CurGlobalSize++;
    pGlobalElem = pGlobalHeap->push(QueueElement(globalKey, value));
  }
  else {
    if (globalKey > pGlobalHeap->findMin().Key) {
      DeleteMinGlobalElem();
      CurGlobalSize++;
      pGlobalElem = pGlobalHeap->push(QueueElement(globalKey, value));
    }
  }
  //push to a local queue
  if (!IsLocalFull()) {
    CurLocalSize++;
    pLocalElem = pLocalHeap->push(QueueElement(localKey, value));
  }
  else {
    if (localKey > pLocalHeap->findMin().Key) {
      DeleteMinLocalElem();
      CurLocalSize++;
      pLocalElem = pLocalHeap->push(QueueElement(localKey, value));
    }
  }
  //link elements
  if (pGlobalElem != NULL && pLocalElem != NULL) {
    pGlobalElem->pLinkedElement = pLocalElem;
    pLocalElem->pLinkedElement = pGlobalElem;
  }

  return pGlobalElem;
}

// ------------------------------------------------------------------------------------------------
QueueElement* PriorityDualQueue::PushWithPriority(double globalKey, double localKey, void * value)
{
  QueueElement* pGlobalElem = NULL, *pLocalElem = NULL;
  //push to a global queue
  if (!IsEmpty()) {
    if (globalKey >= pGlobalHeap->findMin().Key) {
      if (IsFull())
        DeleteMinGlobalElem();
      CurGlobalSize++;
      pGlobalElem = pGlobalHeap->push(QueueElement(globalKey, value));
    }
  }
  else {
    CurGlobalSize++;
    pGlobalElem = pGlobalHeap->push(QueueElement(globalKey, value));
  }
  //push to a local queue
  if (!IsLocalEmpty()) {
    if (localKey >= pLocalHeap->findMin().Key) {
      if (IsLocalFull())
        DeleteMinLocalElem();
      CurLocalSize++;
      pLocalElem = pLocalHeap->push(QueueElement(localKey, value));
    }
  }
  else {
    CurLocalSize++;
    pLocalElem = pLocalHeap->push(QueueElement(localKey, value));
  }
  //link elements
  if (pGlobalElem != NULL && pLocalElem != NULL) {
    pGlobalElem->pLinkedElement = pLocalElem;
    pLocalElem->pLinkedElement = pGlobalElem;
  }

  return pGlobalElem;
}

// ------------------------------------------------------------------------------------------------
void PriorityDualQueue::PopFromLocal(double * key, void ** value)
{
  if (CurLocalSize != 0)
  {
    QueueElement tmp = pLocalHeap->popMax();
    *key = tmp.Key;
    *value = tmp.pValue;
    CurLocalSize--;

    //delete linked element from the global queue
    if (tmp.pLinkedElement != NULL) {
      pGlobalHeap->deleteElement(tmp.pLinkedElement);
      CurGlobalSize--;
    }
  }
  else
  {
    throw EXCEPTION("Cannot pop element from empty queue");
  }
}

// ------------------------------------------------------------------------------------------------
void PriorityDualQueue::Pop(double * key, void ** value)
{
  if (CurGlobalSize != 0)
  {
    QueueElement tmp = pGlobalHeap->popMax();
    *key = tmp.Key;
    *value = tmp.pValue;
    CurGlobalSize--;

    //delete linked element from the local queue
    if (tmp.pLinkedElement != NULL)
    {
      pLocalHeap->deleteElement(tmp.pLinkedElement);
      CurLocalSize--;
    }
  }
  else
  {
    throw EXCEPTION("Cannot pop element from empty queue");
  }
}

// ------------------------------------------------------------------------------------------------
void PriorityDualQueue::DeleteByValue(void *value)
{
  QueueElement* globalHeapMem = pGlobalHeap->getHeapMemPtr();
  for (int i = 0; i < CurGlobalSize; i++)
    if (globalHeapMem[i].pValue == value)
    {
      //delete linked element from the local queue
      if (globalHeapMem[i].pLinkedElement != NULL)
      {
        pLocalHeap->deleteElement(globalHeapMem[i].pLinkedElement);
        CurLocalSize--;
      }
      CurGlobalSize--;
      pGlobalHeap->deleteElement(globalHeapMem + i);
      return;
    }

  //if the value exists only in local queue
  QueueElement* localHeapMem = pLocalHeap->getHeapMemPtr();
  for (int i = 0; i < CurLocalSize; i++)
    if (localHeapMem[i].pValue == value)
    {
      CurLocalSize--;
      pLocalHeap->deleteElement(localHeapMem + i);
      break;
    }
}

// ------------------------------------------------------------------------------------------------
void PriorityDualQueue::DeleteElement(QueueElement* item)
{
  DeleteByValue(item->pValue);
  //QueueElement* heapMem = pMem->getHeapMemPtr();
  //pMem->deleteElement(item);
  //CurSize--;
}

// ------------------------------------------------------------------------------------------------
void PriorityDualQueue::Clear()
{
  ClearLocal();
  ClearGlobal();
}

// ------------------------------------------------------------------------------------------------
void PriorityDualQueue::Resize(int size)
{
  MaxSize = size;
  CurGlobalSize = CurLocalSize = 0;
  delete pLocalHeap;
  delete pGlobalHeap;
  pLocalHeap = new MinMaxHeap< QueueElement, _less >(MaxSize);
  pGlobalHeap = new MinMaxHeap< QueueElement, _less >(MaxSize);
}

// ------------------------------------------------------------------------------------------------
void PriorityDualQueue::ClearLocal()
{
  pLocalHeap->clear();
  CurLocalSize = 0;
}

// ------------------------------------------------------------------------------------------------
void PriorityDualQueue::ClearGlobal()
{
  pGlobalHeap->clear();
  CurGlobalSize = 0;
}
// - end of file ----------------------------------------------------------------------------------