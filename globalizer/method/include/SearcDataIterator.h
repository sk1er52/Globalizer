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

#ifndef __SEARC_DATA_ITERATOR_H__
#define __SEARC_DATA_ITERATOR_H__

#include "SearchData.h"
#include "TreeNode.h"



// ------------------------------------------------------------------------------------------------
class SearcDataIterator
{
  friend class SearchData;
protected:
  SearchData *pContainer;
  TreeNode *pObject;

public:
  SearcDataIterator();

  SearcDataIterator& operator++();
  SearcDataIterator operator++(int);
  SearcDataIterator& operator--();
  SearcDataIterator operator--(int);
  operator void*() const;
  SearchInterval* operator->();
  SearchInterval* operator *() const;
};



#endif
// - end of file ----------------------------------------------------------------------------------
