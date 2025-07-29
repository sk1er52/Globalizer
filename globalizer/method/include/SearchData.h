/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      data.h                                                      //
//                                                                         //
//  Purpose:   Header file for search data classes                         //
//                                                                         //
//  Author(s): Sysoyev A., Barkalov K., Sovrasov V.                        //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __SEARCH_DATA_H__
#define __SEARCH_DATA_H__

#include "Common.h"
#include "Extended.h"
#include "DualQueue.h"
#include "Queue.h"
#include <stack>
#include "Parameters.h"

#include <string.h>

#include <list>
#include <vector>


struct TreeNode;
class SearchInterval;
class SearcDataIterator;
class Trial;
// ------------------------------------------------------------------------------------------------
class TSearchData
{
  friend class SearcDataIterator;
protected:
  /// число функций задачи
  int NumOfFuncs;
  /// максимальный размер МСП = максимальному числу итераций метода
  int MaxSize;
  /// текущее число интервалов в дереве
  int Count;
  /// текущий индекс, используется в итераторе
  int CurIndex;
  /// Корень дерева
  TreeNode *pRoot;
  /// Текущая вершина дерева
  TreeNode *pCur;

  //TreeNode *pCurIter;
  /// стек для итератора
  std::stack<TreeNode*> Stack;
  ///очередь характеристик
  PriorityQueueCommon *pQueue;

  /// список всех точек, для их последующего удаления
  std::vector<Trial*> trials;

  /// истина, если нужен пересчет характеристик
  bool recalc;

  /// Лучшая точка, полученная для данной поисковой информации
  Trial* BestTrial;

  void DeleteTree(TreeNode *pNode);
  unsigned char GetHeight(TreeNode *p);
  int GetBalance(TreeNode *p);
  void FixHeight(TreeNode *p);
  TreeNode* RotateRight(TreeNode *p); // правый поворот вокруг p
  TreeNode* RotateLeft(TreeNode *p);  // левый поворот вокруг p
  TreeNode* Balance(TreeNode *p);     // балансировка узла p
  TreeNode* Maximum(TreeNode *p) const; //поиск самого левого интервала в поддереве
  TreeNode* Minimum(TreeNode *p) const; //поиск самого правого интервала в поддереве
  TreeNode* Previous(TreeNode *p) const; //получение предыдущего и следующего за p интервалов
  TreeNode* Next(TreeNode *p) const;
  // вставка в дерево с корнем p (рекурсивная)
  TreeNode* Insert(TreeNode *p, SearchInterval &pInterval);
  // поиск узла с нужным x в дереве с корнем p по левой границе интервала (рекурсивный)
  TreeNode* Find(TreeNode *p, Trial* x) const;
  // поиск узла по правой границе интервала
  TreeNode* FindR(TreeNode *p, Trial* x) const;
  // поиск узла с нужным x по левой и правой границам интервала (рекурсивный)
  //   xl() < x < xr
  TreeNode* FindIn(TreeNode *p, Trial* x) const;
public:
  /// Вектор указателей на матрицы состояния поиска, для которых нужно произвести пересчет
  static std::vector<TSearchData*> pRecalcDatas;

  TSearchData(int _NumOfFuncs, int _MaxSize = DefaultSearchDataSize);
  TSearchData(int _NumOfFuncs, int _MaxSize, int _queueSize);
  ~TSearchData();

  /// Очищает и дерево и очередь интервалов
  void Clear();
  /// новый интервал (по xl)
  SearchInterval* InsertInterval(SearchInterval &pInterval); 
  /// обновление интервала (по xl)
  void UpdateInterval(SearchInterval &pInterval); 
  /// Ищет интервал у которого левой точкой является x
  SearchInterval* GetIntervalByX(Trial* x);
  /** Поиск интервала, в котором содержится x, т.е. xl() < x < xr
     нужен для вставки прообразов при использовании множественной развертки */
  SearchInterval* FindCoveringInterval(Trial* x);
  /** Получение интервала с максимальной хар-кой. Интервал берется из очереди. Если очередь пуста,
     то сначала будет вызван Refill() */
  SearchInterval* GetIntervalWithMaxR();
  /** Получение интервала с максимальной локалььной хар-кой. Интервал берется из очереди. Если очередь пуста,
  то сначала будет вызван Refill() */
  SearchInterval* GetIntervalWithMaxLocalR();

  /** Вставка испытания в заданный интервал. Нужна для множественной развертки и
  добавления поисковой информации локального метода
  возвращает указатель на интервал с левым концом в newPoint
  */
  SearchInterval* InsertPoint(SearchInterval* coveringInterval, Trial& newPoint,
    int iteration, int methodDimension);

  /// Итератор
  SearcDataIterator GetIterator(SearchInterval* p);
  SearcDataIterator GetBeginIterator();

  /** Получить интервал, предыдущий к указанному
     не относится к итератору, не меняет текущий узел в итераторе
    SearchInterval* GetPrev(SearchInterval &pInterval);

   Для работы с очередью характеристик
   вставка, если новый элемент больше минимального в очереди
     если при этом очередь полна, то замещение минимального
  */
  void PushToQueue(SearchInterval *pInterval);
  /// Перезаполнение очереди (при ее опустошении или при смене оценки константы Липшица)
  void RefillQueue();
  /// Удалить интервал из очереди
  void DeleteIntervalFromQueue(SearchInterval* i);

  /// Берет из очереди один интервал
  void PopFromGlobalQueue(SearchInterval **pInterval);
  /// Берет из очереди локальных характеристик один интервал
  void PopFromLocalQueue(SearchInterval **pInterval);
  /// Очистить очередь интервалов
  void ClearQueue();
  /// Изменить размер очереди интервалов
  void ResizeQueue(int size);

  /// Возвращает текущее число интервалов в дереве 
  int GetCount();

  /// оценки констант Липшица
  double M[MaxNumOfFunc];
  /// минимальные значения функций задачи (для индексного метода)
  double Z[MaxNumOfFunc];

  void GetBestIntervals(SearchInterval** intervals, int count);
  void GetBestLocalIntervals(SearchInterval** intervals, int count);
  std::vector<Trial*>& GetTrials()
  {
    return trials;
  }

  /// возвращает максимальный элемент без извлечения
  SearchInterval& FindMax();

  /// истина, если нужен пересчет характеристик
  bool IsRecalc()
  {
    return recalc;
  }
  /// Задает нужно ли пересчитывать характеристики
  void SetRecalc(bool f)
  {    
    if (recalc == false)
      pRecalcDatas.push_back(this);
    recalc = f;    
  }
  /// Лучшая точка, полученная для данной поисковой информации
  Trial* GetBestTrial()
  {
    return BestTrial;
  }
  /// Задает лучшую точку
  void SetBestTrial(Trial* trial);

  /// Всплытие для интервала
  void TrickleUp(SearchInterval* intervals);

  /// Возвращает размер очереди
  int GetQueueSize();


  double local_r;//вычисляемое r
};


#endif
// - end of file ----------------------------------------------------------------------------------
