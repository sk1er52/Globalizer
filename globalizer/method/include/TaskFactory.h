/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      method_factory.h                                            //
//                                                                         //
//  Purpose:   Header file for method factory class                        //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __TASK_FACTORY_H__
#define __TASK_FACTORY_H__

#include "Task.h"

class TaskFactory
{
public:
  static int num;
  static std::vector<int> permutations;
  static Task* CreateTask(IProblem* _problem, int _ProcLevel);
  static Task* CreateTask();
  static Task* CreateTask(Task* t);
};

#endif
// - end of file ----------------------------------------------------------------------------------
