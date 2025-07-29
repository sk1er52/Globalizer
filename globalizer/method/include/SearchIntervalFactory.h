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

#ifndef __SEARCH_INTERVAL_FACTORY_H__
#define __SEARCH_INTERVAL_FACTORY_H__

#include "SearchInterval.h"

// ------------------------------------------------------------------------------------------------
class SearchInterval;
class SearchIntervalFactory
{
public:
  static SearchInterval* CreateSearchInterval(SearchInterval& interval);
  static SearchInterval* CreateSearchInterval();
};



#endif
// - end of file ----------------------------------------------------------------------------------
