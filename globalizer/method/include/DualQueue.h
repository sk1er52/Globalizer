/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      dual_queue.h                                                //
//                                                                         //
//  Purpose:   Header file for priority dual queue class                   //
//                                                                         //
//  Author(s): Sysoyev A., Barkalov K., Sovrasov V.                        //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __DUAL_QUEUE_H__
#define __DUAL_QUEUE_H__

#include "MinMaxHeap.h"
#include "Common.h"
#include "QueueCommon.h"

class PriorityDualQueue : public PriorityQueueCommon
{
protected:
  int MaxSize;
  int CurLocalSize;
  int CurGlobalSize;

  MinMaxHeap< QueueElement, _less >* pGlobalHeap;
  MinMaxHeap< QueueElement, _less >* pLocalHeap;

  void DeleteMinLocalElem();
  void DeleteMinGlobalElem();
  void ClearLocal();
  void ClearGlobal();
public:

  PriorityDualQueue(int _MaxSize = DefaultQueueSize); // _MaxSize must be qual to 2^k - 1
  ~PriorityDualQueue();

  int GetLocalSize() const;
  int GetSize() const;
  int GetMaxSize() const;
  bool IsLocalEmpty() const;
  bool IsLocalFull() const;
  bool IsEmpty() const;
  bool IsFull() const;

  QueueElement* Push(double globalKey, double localKey, void *value);
  QueueElement* PushWithPriority(double globalKey, double localKey, void *value);
  void Pop(double *key, void **value);
  void DeleteByValue(void *value);
  /// Удаляет элемент
  virtual void DeleteElement(QueueElement * item);
  void PopFromLocal(double *key, void **value);

  void Clear();
  void Resize(int size);
  virtual QueueElement& FindMax()
  {
    return pGlobalHeap->findMax();
  }
  void TrickleUp(QueueElement * item)
  {
    pGlobalHeap->TrickleUp(item);
  }

  void TrickleDown(QueueElement * item)
  {
    pGlobalHeap->TrickleDown(item);
  }
};
#endif
// - end of file ----------------------------------------------------------------------------------