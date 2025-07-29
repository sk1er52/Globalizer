/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2013 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      Data.cpp                                                    //
//                                                                         //
//  Purpose:   Source file for search data classes                         //
//                                                                         //
//  Author(s): Sysoyev A., Barkalov K., Sovrasov V.                        //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////


#include "Exception.h"
#include "SearchData.h"
#include "SearchIntervalFactory.h"
#include "SearcDataIterator.h"
#include "Trial.h"
#include "TreeNode.h"
#include "TrialFactory.h"

//bool operator<(const SearchInterval& i1, const SearchInterval& i2)
//{
//  return (i1 < i2);
//}

  /// Вектор указателей на матрицы состояния поиска, для которых нужно произвести пересчет
std::vector<TSearchData*> TSearchData::pRecalcDatas;

// ------------------------------------------------------------------------------------------------
// TSearchData Methods
// ------------------------------------------------------------------------------------------------
TSearchData::TSearchData(int _NumOfFuncs, int _MaxSize)
{
  if (_MaxSize <= 0)
  {
    throw EXCEPTION("MaxSize of SearchData is out of range");
  }
  MaxSize = _MaxSize;
  Count = 0;
  pRoot = pCur = NULL;

  if (parameters.localMix == 0)
    pQueue = new PriorityQueue(DefaultQueueSize);
  else
    pQueue = new PriorityDualQueue(DefaultQueueSize);

  if ((_NumOfFuncs > MaxNumOfFunc) || (_NumOfFuncs <= 0))
  {
    throw EXCEPTION("NumOfFunc is out of range");
  }
  NumOfFuncs = _NumOfFuncs;
  for(int i=0;i<NumOfFuncs;i++)
  {
    M[i] = i < parameters.M_constant.GetSize() ? parameters.M_constant[i] : 1.; //Начальное значение оценок констант Липшица равно 1
    Z[i] = MaxDouble;  //Начальное значение минимумов - максимальное
  }

  BestTrial = 0;
  pRecalcDatas.clear();
}

TSearchData::TSearchData(int _NumOfFuncs, int _MaxSize, int _queueSize)
{
  if (_MaxSize <= 0)
  {
    throw EXCEPTION("MaxSize of SearchData is out of range");
  }
  MaxSize = _MaxSize;
  Count = 0;
  pRoot = pCur = NULL;

  if (parameters.localMix == 0)
    pQueue = new PriorityQueue(_queueSize);
  else
    pQueue = new PriorityDualQueue(_queueSize);

  if ((_NumOfFuncs > MaxNumOfFunc) || (_NumOfFuncs <= 0))
  {
    throw EXCEPTION("NumOfFunc is out of range");
  }
  NumOfFuncs = _NumOfFuncs;
  for (int i = 0;i<NumOfFuncs;i++)
  {
    M[i] = i < parameters.M_constant.GetSize() ? parameters.M_constant[i] : 1.; //Начальное значение оценок констант Липшица равно 1
    Z[i] = MaxDouble;  //Начальное значение минимумов - максимальное
  }

  BestTrial = 0;
  pRecalcDatas.clear();
}

// ------------------------------------------------------------------------------------------------
TSearchData::~TSearchData()
{
  DeleteTree(pRoot);
  delete pQueue;

  for (unsigned int i = 0; i < trials.size(); i++)
  {
    delete ((Trial*)(trials[i]));
  }
}

// ------------------------------------------------------------------------------------------------
void TSearchData::Clear()
{
  DeleteTree(pRoot);
  Count = 0;
  pRoot = pCur = NULL;
  ClearQueue();
  for (int i = 0; i < NumOfFuncs; i++)
  {
    M[i] = i < parameters.M_constant.GetSize() ? parameters.M_constant[i] : 1.;      //Начальное значение оценок констант Липшица равно 1
    Z[i] = MaxDouble;  //Начальное значение минимумов - максимальное
  }
  for (unsigned int i = 0; i < trials.size(); i++)
  {
    delete ((Trial*)(trials[i]));
  }

  trials.clear();

  BestTrial = 0;
}

// ------------------------------------------------------------------------------------------------
void TSearchData::GetBestIntervals(SearchInterval** intervals, int count)
{
  //Текущая длина очереди характеристик должна быть больше, чем число интервалов, которые из нее надо извлечь
  //Иначе очередь надо перезаполнить
  if (pQueue->GetSize() <= count)
    RefillQueue();

  //Здесь надо взять count лучших характеристик из глобальной очереди
  for (int i = 0; i < count; i++)
    intervals[i] = GetIntervalWithMaxR();
}

// ------------------------------------------------------------------------------------------------
void TSearchData::GetBestLocalIntervals(SearchInterval** intervals, int count)
{
  //Вариант решения проблемы - текущая длина очереди характеристик должна быть больше, чем число интервалов, которые из нее надо извлечь!
  if (pQueue->GetSize() <= count)
    RefillQueue();

  //Здесь надо взять count лучших характеристик из локальной очереди
  for (int i = 0; i < count; i++)
    intervals[i] = GetIntervalWithMaxLocalR();
}

// ------------------------------------------------------------------------------------------------
void TSearchData::DeleteTree(TreeNode *pNode)
{
  if (pNode != NULL)
  {
    DeleteTree(pNode->pLeft);
    DeleteTree(pNode->pRight);
    delete pNode;
  }
}

// ------------------------------------------------------------------------------------------------
unsigned char TSearchData::GetHeight(TreeNode *p)
{
  return p ? p->Height : 0;
}

// ------------------------------------------------------------------------------------------------
int TSearchData::GetBalance(TreeNode *p)
{
  return GetHeight(p->pRight) - GetHeight(p->pLeft);
}

// ------------------------------------------------------------------------------------------------
void TSearchData::FixHeight(TreeNode *p)
{
  unsigned char hl = GetHeight(p->pLeft);
  unsigned char hr = GetHeight(p->pRight);
  p->Height = (hl > hr ? hl : hr) + 1;
}

void TSearchData::SetBestTrial(Trial* trial)
{
  BestTrial = trial;

  // Обновить текущие значение минимумов
  for (int v = 0; v <= BestTrial->index; v++)
  {
    if (v < BestTrial->index)
    {
      Z[v] = -M[v] * parameters.rEps;
    }
    else
    {
      if (BestTrial->FuncValues[v] != MaxDouble)
        Z[v] = BestTrial->FuncValues[v];
      else
        Z[v] = 0;
    }
  }
}

// ------------------------------------------------------------------------------------------------
TreeNode* TSearchData::RotateRight(TreeNode *p)
{
  TreeNode *pTmp = p->pLeft;
  p->pLeft = pTmp->pRight;
  if (p->pLeft)
    p->pLeft->pParent = p;
  pTmp->pParent = p->pParent;
  pTmp->pRight = p;
  p->pParent = pTmp;

  FixHeight(p);
  FixHeight(pTmp);

  return pTmp;
}

// ------------------------------------------------------------------------------------------------
TreeNode* TSearchData::RotateLeft(TreeNode *p)
{
  TreeNode *pTmp = p->pRight;
  p->pRight = pTmp->pLeft;
  if (p->pRight)
    p->pRight->pParent = p;
  pTmp->pParent = p->pParent;
  pTmp->pParent = p->pParent;
  pTmp->pLeft = p;
  p->pParent = pTmp;

  FixHeight(p);
  FixHeight(pTmp);

  return pTmp;
}

// ------------------------------------------------------------------------------------------------
TreeNode* TSearchData::Maximum(TreeNode *p) const
{
  while (p->pRight != NULL)
    p = p->pRight;
  return p;
}

// ------------------------------------------------------------------------------------------------
TreeNode* TSearchData::Minimum(TreeNode *p) const
{
  while (p->pLeft != NULL)
    p = p->pLeft;
  return p;
}

// ------------------------------------------------------------------------------------------------
TreeNode* TSearchData::Balance(TreeNode *p)
{
  FixHeight(p);
  if (GetBalance(p) == 2)
  {
    if (GetBalance(p->pRight) < 0)
    {
      p->pRight = RotateRight(p->pRight);
      //p->pRight->pParent = p;
    }
    return RotateLeft(p);
  }
  if (GetBalance(p) == -2)
  {
    if (GetBalance(p->pLeft) > 0)
    {
      p->pLeft = RotateLeft(p->pLeft);
      //p->pLeft->pParent = p;
    }
    return RotateRight(p);
  }
  return p; // балансировка не нужна
}

// ------------------------------------------------------------------------------------------------
TreeNode* TSearchData::Insert(TreeNode *p, SearchInterval &pInterval)
{
  if (!p)
  {
    // не безопасно, если будет несколько потоков
    pCur = new TreeNode(pInterval);
    return pCur;
  }

  if (pInterval == *p->pInterval)
    throw EXCEPTION("Interval already exists!");
  if (pInterval < *p->pInterval)
  {
    p->pLeft = Insert(p->pLeft, pInterval);
    p->pLeft->pParent = p;
  }
  else //  if (pInterval.xl() > p->pInterval->xl)
  {
    p->pRight = Insert(p->pRight, pInterval);
    p->pRight->pParent = p;
  }

  return Balance(p);
}

// ------------------------------------------------------------------------------------------------
TreeNode* TSearchData::Find(TreeNode *p, Trial* x) const
{
  TreeNode *res;
  if (!p)
    return p;
  if (*x < *(p->pInterval->LeftPoint))
    res = Find(p->pLeft, x);
  else if (*x > *(p->pInterval->LeftPoint))
    res = Find(p->pRight, x);
  else
    res = p;

  return res;
}

// ------------------------------------------------------------------------------------------------
TreeNode* TSearchData::FindR(TreeNode *p, Trial* x) const
{
  TreeNode *res;
  if (!p)
    return p;
  if (*x < *(p->pInterval->LeftPoint))
    res = FindR(p->pLeft, x);
  else if (*x > * (p->pInterval->RightPoint))
    res = FindR(p->pRight, x);
  else
    res = p;

  return res;
}

// ------------------------------------------------------------------------------------------------
TreeNode * TSearchData::FindIn(TreeNode * p, Trial* x) const
{
  TreeNode *res;
  if (!p)
    return p;
  if (*x < *(p->pInterval->LeftPoint))
    res = FindIn(p->pLeft, x);
  else if (*x > * (p->pInterval->RightPoint))
    res = FindIn(p->pRight, x);
  else
    res = p;

  return res;
}

// ------------------------------------------------------------------------------------------------
SearchInterval* TSearchData::InsertInterval(SearchInterval &pInterval)
{
  if (pInterval.xl() >= pInterval.xr())
  {
    throw EXCEPTION("Cannot insertinterval with not positive length.");
  }
  pRoot = Insert(pRoot, pInterval);
  Count++;
  return pCur->pInterval;
}

// ------------------------------------------------------------------------------------------------
void TSearchData::UpdateInterval(SearchInterval &pInterval)
{
  pCur = Find(pRoot, pInterval.LeftPoint);
  if (pCur)
  {
    (*pCur->pInterval) = pInterval;
  }
}

// ------------------------------------------------------------------------------------------------
SearchInterval* TSearchData::GetIntervalByX(Trial* x)
{
  pCur = Find(pRoot, x);
  if (pCur)
    return pCur->pInterval;
  else
    return NULL;
}

// ------------------------------------------------------------------------------------------------
SearchInterval * TSearchData::FindCoveringInterval(Trial* x)
{
  pCur = FindIn(pRoot, x);
  if (pCur)
    return pCur->pInterval;
  else
    return NULL;
}

// ------------------------------------------------------------------------------------------------
TreeNode* TSearchData::Previous(TreeNode *p) const
{
  if (p == NULL)
    return NULL;
  if (p->pLeft != NULL)
    return Maximum(p->pLeft);
  else {
    TreeNode* tmp = p;
    TreeNode* parentTmp = p->pParent;
    while (parentTmp != NULL) {
      if (parentTmp->pLeft != tmp)
        break;
      tmp = parentTmp;
      parentTmp = tmp->pParent;
    }
    return parentTmp;
  }
}

// ------------------------------------------------------------------------------------------------
TreeNode* TSearchData::Next(TreeNode *p) const
{
  if (p == NULL)
    return NULL;
  if (p->pRight != NULL)
    return Minimum(p->pRight);
  else {
    TreeNode* tmp = p;
    TreeNode* parentTmp = p->pParent;
    while (parentTmp != NULL) {
      if (parentTmp->pRight != tmp)
        break;
      tmp = parentTmp;
      parentTmp = tmp->pParent;
    }
    return parentTmp;
  }
}

// ------------------------------------------------------------------------------------------------
SearchInterval* TSearchData::InsertPoint(SearchInterval* coveringInterval,
  Trial& newPoint, int iteration, int methodDimension)
{
  // Если точка уже есть в базе, то ничего не делаем
  if (newPoint == *(coveringInterval->LeftPoint) || newPoint == *(coveringInterval->RightPoint))
    return NULL;

  //SearchInterval NewInterval;

  // правый подинтервал
  SearchInterval* NewInterval = SearchIntervalFactory::CreateSearchInterval();

  NewInterval->ind = iteration;
  NewInterval->K = newPoint.K;

  NewInterval->LeftPoint = &newPoint;

  NewInterval->RightPoint = coveringInterval->RightPoint;

  NewInterval->delta = root(NewInterval->xr() - NewInterval->xl(), methodDimension);
  // Корректируем существующий интервал
  coveringInterval->RightPoint = NewInterval->LeftPoint;

  coveringInterval->delta = root(coveringInterval->xr() - coveringInterval->xl(), methodDimension);
  // Интервал сформирован - можно добавлять
  SearchInterval* p = InsertInterval(*NewInterval);

  newPoint.leftInterval = coveringInterval;
  newPoint.rightInterval = p;

  newPoint.leftInterval->LeftPoint->rightInterval = coveringInterval;
  newPoint.rightInterval->RightPoint->leftInterval = p;

  return p;
}

// ------------------------------------------------------------------------------------------------
SearcDataIterator TSearchData::GetIterator(SearchInterval* p)
{
  SearcDataIterator iter;
  iter.pContainer = this;
  iter.pObject = Find(pRoot, p->LeftPoint);
  return iter;
}

// ------------------------------------------------------------------------------------------------
SearcDataIterator TSearchData::GetBeginIterator()
{
  SearcDataIterator iter;
  iter.pContainer = this;
  iter.pObject = Minimum(pRoot);
  return iter;
}

// ------------------------------------------------------------------------------------------------
SearchInterval* TSearchData::GetIntervalWithMaxR()
{
  SearchInterval *pRes = NULL;
  //Если нужно вынуть два элемента, и текущий размер очереди - 2, то первый поток вынет максимальный элемент,
  //а потом очередь перестроится по предыдущей МСП, и следующим будет вынут тот же самый элемент!!!
  // Эта ситуация учтена в функции GetBestIntervals: очередь не будет короче числа интервалов, которые из нее надо извлечь на 1-й итерации метода
  if (pQueue->IsEmpty() || pQueue->GetSize() == 1)
    RefillQueue();

  PopFromGlobalQueue(&pRes);
  return pRes;
}

// ------------------------------------------------------------------------------------------------
SearchInterval* TSearchData::GetIntervalWithMaxLocalR()
{
  SearchInterval *pRes = NULL;
  if (((PriorityDualQueue*)pQueue)->IsLocalEmpty()
    || ((PriorityDualQueue*)pQueue)->GetLocalSize() == 1)
    RefillQueue();

  PopFromLocalQueue(&pRes);
  return pRes;
}

// ------------------------------------------------------------------------------------------------
void TSearchData::RefillQueue()
{
  pQueue->Clear(); // Если очередь уже была пуста, то вызов лишний, но он из одного действия, так
                   //   что ставить проверку нет смысла
  // Обход всех интервалов, вставка их в очередь
  //   Push() сначала заполняет очередь, потом вставляет с замещением

  QueueElement* queueElementa = 0;
  for (SearcDataIterator it = GetBeginIterator(); it; ++it)
  {
    queueElementa = pQueue->Push(it->R, it->locR, *it);
    if (queueElementa != 0)
    {
      it->SetQueueElementa(queueElementa);
    }
  }
}


// ------------------------------------------------------------------------------------------------
void TSearchData::DeleteIntervalFromQueue(SearchInterval * i)
{
  if (i->GetQueueElementa() == 0)
  {
    pQueue->DeleteByValue(i);
    i->SetQueueElementa(0);
  }
  else
  {
    pQueue->DeleteElement(i->GetQueueElementa());
    i->SetQueueElementa(0);
  }
}

// ------------------------------------------------------------------------------------------------
SearchInterval& TSearchData::FindMax()
{
  if (pQueue->IsEmpty())
  {
    RefillQueue();
  }
  QueueElement* tmp = &pQueue->FindMax();
  return *((SearchInterval*)tmp->pValue);
}

// ------------------------------------------------------------------------------------------------
void TSearchData::TrickleUp(SearchInterval * intervals)
{
  if (intervals->GetQueueElementa() != 0)
  {
    pQueue->TrickleUp(intervals->GetQueueElementa());

    QueueElement* qe = intervals->GetQueueElementa();
    double oldKey = qe->Key;
    qe->Key = intervals->R;
    if (qe->Key >= oldKey)
      pQueue->TrickleUp(intervals->GetQueueElementa());
    else
      pQueue->TrickleDown(intervals->GetQueueElementa());
  }
}

int TSearchData::GetQueueSize()
{
  return pQueue->GetMaxSize();
}

// ------------------------------------------------------------------------------------------------
void TSearchData::PushToQueue(SearchInterval *pInterval)
{
  if (pInterval == 0)
  {
    throw EXCEPTION("Cannot push NULL pointer to queue.");
  }
  //На начальном этапе в очередь записываем все данные
  QueueElement* queueElementa = 0;
  if (Count <= pQueue->GetMaxSize())
    queueElementa = pQueue->Push(pInterval->R, pInterval->locR, pInterval);
  else
    queueElementa = pQueue->PushWithPriority(pInterval->R, pInterval->locR, pInterval);

  if (queueElementa != 0)
  {
    pInterval->SetQueueElementa(queueElementa);
  }
}

// ------------------------------------------------------------------------------------------------
void TSearchData::PopFromGlobalQueue(SearchInterval **pInterval)
{
  void *p;
  double R;
  pQueue->Pop(&R, &p);
  *pInterval = (SearchInterval*)p;
  (*pInterval)->SetQueueElementa(0);
}

// ------------------------------------------------------------------------------------------------
void TSearchData::PopFromLocalQueue(SearchInterval **pInterval)
{
  void *p;
  double R;
  ((PriorityDualQueue*)pQueue)->PopFromLocal(&R, &p);
  *pInterval = (SearchInterval*)p;
  (*pInterval)->SetQueueElementa(0);
}

// ------------------------------------------------------------------------------------------------
void TSearchData::ClearQueue()
{
  pQueue->Clear();
}
// ------------------------------------------------------------------------------------------------
void TSearchData::ResizeQueue(int size)
{
  pQueue->Resize(size);
}
// ------------------------------------------------------------------------------------------------
int TSearchData::GetCount()
{
  return Count;
}

// ------------------------------------------------------------------------------------------------
SearcDataIterator::SearcDataIterator() :
  pContainer(NULL), pObject(NULL)
{}

// ------------------------------------------------------------------------------------------------
SearcDataIterator & SearcDataIterator::operator++()
{
  pObject = pContainer->Next(pObject);
  return *this;
}

// ------------------------------------------------------------------------------------------------
SearcDataIterator SearcDataIterator::operator++(int)
{
  SearcDataIterator tmp = *this;
  ++(*this);
  return tmp;
}

// ------------------------------------------------------------------------------------------------
SearcDataIterator & SearcDataIterator::operator--()
{
  pObject = pContainer->Previous(pObject);
  return *this;
}

// ------------------------------------------------------------------------------------------------
SearcDataIterator SearcDataIterator::operator--(int)
{
  SearcDataIterator tmp = *this;
  --(*this);
  return tmp;
}

// ------------------------------------------------------------------------------------------------
SearcDataIterator::operator void*() const
{
  if (pContainer && pObject)
    return pObject;
  else
    return NULL;
}

// ------------------------------------------------------------------------------------------------
SearchInterval * SearcDataIterator::operator->()
{
  return pObject->pInterval;
}

// ------------------------------------------------------------------------------------------------
SearchInterval* SearcDataIterator::operator*() const
{
  if (pObject)
    return pObject->pInterval;
  else
    return NULL;
}

// ------------------------------------------------------------------------------------------------
SearchInterval* SearchIntervalFactory::CreateSearchInterval(SearchInterval& interval)
{
  return new SearchInterval(interval);
}

// ------------------------------------------------------------------------------------------------
SearchInterval* SearchIntervalFactory::CreateSearchInterval()
{
  return new SearchInterval();
}

// ------------------------------------------------------------------------------------------------
SearchInterval::SearchInterval()
{
  LeftPoint = 0;
  RightPoint = 0;
  R = locR = 0;
  ind = 0;
  // Выделение памяти происходит в момент добавления в матрицу
  K = 0;
  delta = 0;
  queueElementa = 0;
  treeNode = 0;

  status = SearchInterval::usual;
}

// ------------------------------------------------------------------------------------------------
SearchInterval::SearchInterval(const SearchInterval &p)
{

  LeftPoint = p.LeftPoint;//->Clone();
  RightPoint = p.RightPoint;//->Clone();

  ind = p.ind;
  K = p.K;
  R = p.R;
  locR = p.locR;
  delta = p.delta;

  //МКО


  queueElementa = p.queueElementa;
  treeNode = p.treeNode;
  status = p.status;
}

// ------------------------------------------------------------------------------------------------
SearchInterval::~SearchInterval()
{
  //Конфликт с функцией RenewSearchData
  //Решение 1 - статический массив z() в TSearchInterval
  //Решение 2 - выделять память перед записью в массив

}

// ------------------------------------------------------------------------------------------------
void SearchInterval::CreatePoint()
{
  LeftPoint = TrialFactory::CreateTrial();
  RightPoint = TrialFactory::CreateTrial();
  RightPoint->SetX(1);
}
