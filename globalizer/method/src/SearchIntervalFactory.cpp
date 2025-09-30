/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2025 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      SearchIntervalFactory.cpp                                    //
//                                                                         //
//  Purpose:   Source file for search interval factory classes             //
//                                                                         //
//  Author(s): Sysoyev A., Barkalov K., Sovrasov V., Zaitsev A.            //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "SearchIntervalFactory.h"

SearchInterval* SearchIntervalFactory::CreateSearchInterval(SearchInterval& interval)
{
    return new SearchInterval(interval);
}

// ------------------------------------------------------------------------------------------------
SearchInterval* SearchIntervalFactory::CreateSearchInterval()
{
    return new SearchInterval();
}


// - end of file ----------------------------------------------------------------------------------
