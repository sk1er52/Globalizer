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

#ifndef __TRIAL_FACTORY_H__
#define __TRIAL_FACTORY_H__

#include "Parameters.h"
//#include "Trial.h"


class Trial;
class TMultievolventsTrial;
class TrialFactory
{
public:
  //static Trial* CreateTrial(Trial& interval);
  static Trial* CreateTrial()
  {
    if ((parameters.TypeMethod == MultievolventsMethod) || (parameters.TypeMethod == ParallelMultievolventsMethod))
      return new TMultievolventsTrial();
    else
      return new Trial();
  }

  static Trial* CreateTrial(const OBJECTIV_TYPE* startPoint)
  {
    Trial* res;
    if ((parameters.TypeMethod == MultievolventsMethod) || (parameters.TypeMethod == ParallelMultievolventsMethod))
      res = new TMultievolventsTrial();
    else
      res = new Trial();
    memcpy(res->y, startPoint, parameters.Dimension * sizeof(double));

    return res;
  }

  static Trial* CreateTrial(TMultievolventsTrial* point)
  {    
    return new TMultievolventsTrial(*point);
  }

  static Trial* CreateTrial(Trial* point)
  {
    return new Trial(*point);
  }
};


//// ------------------------------------------------------------------------------------------------
//Trial* TrialFactory::CreateTrial(Trial& interval)
//{
//  if (parameters.TypeMethod == MultievolventsMethod)
//    return new TMultievolventsTrial((SearchInterval&)(interval));
//  else
//    return new SearchInterval(interval);
//}

//// ------------------------------------------------------------------------------------------------
//Trial* TrialFactory::CreateTrial()
//{
//  if (parameters.TypeMethod == MultievolventsMethod)
//    return new TMultievolventsTrial();
//  else
//    return new Trial();
//}


#endif
// - end of file ----------------------------------------------------------------------------------
