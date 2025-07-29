/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      Queue.h                                                     //
//                                                                         //
//  Purpose:   Header file for priority queue class                        //
//                                                                         //
//  Author(s): Sysoyev A., Barkalov K., Sovrasov V.                        //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __QUEUE_H__
#define __QUEUE_H__

#include "MinMaxHeap.h"
#include "Common.h"
#include "QueueCommon.h"

class PriorityQueue : public PriorityQueueCommon
{
protected:
  int MaxSize;
  int CurSize;
  MinMaxHeap< QueueElement, _less >* pMem;

  int GetIndOfMinElem();
  void DeleteMinElem();

public:
  PriorityQueue(int _MaxSize = DefaultQueueSize); // _MaxSize must be qual to 2^k - 1
  ~PriorityQueue();

  int GetSize() const;
  int GetMaxSize() const;
  bool IsEmpty() const;
  bool IsFull() const;

  //localKey value is not really used by the Push and PushWithPriority methods
  QueueElement* Push(double globalKey, double localKey, void *value);
  QueueElement* PushWithPriority(double globalKey, double localKey, void *value);
  void Pop(double *key, void **value);
  void DeleteByValue(void *value);

  /// Óäàëÿåò ýëåìåíò
  void DeleteElement(QueueElement * item);

  void Clear();
  void Resize(int size);

  QueueElement& FindMax();
  void TrickleUp(QueueElement * item);
  void TrickleDown(QueueElement * item)
  {
    pMem->TrickleDown(item);
  }
};
#endif
// - end of file ----------------------------------------------------------------------------------
