/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      performance.h                                               //
//                                                                         //
//  Purpose:   Performance measurement class                               //
//                                                                         //
//  Author(s): Sysoyev A.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __PERFORMANCE_H__
#define __PERFORMANCE_H__

#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
 //TODO: add Linux specific things
#endif

// ------------------------------------------------------------------------------------------------
class Performance
{
private:
#ifdef WIN32
  LARGE_INTEGER startCount;
  LARGE_INTEGER s_freq;
#else
// TODO: add Linux specific things
  struct timeval tv1,tv2,dtv;
  struct timezone tz;
#endif

public:
  Performance()
  {
#ifdef WIN32
    s_freq.QuadPart = 0;
#else
// TODO: add Linux specific things
#endif
  }

  void Start()
  {
#ifdef WIN32
    QueryPerformanceCounter(&startCount);
#else
// TODO: add Linux specific things
    gettimeofday(&tv1, &tz);
#endif
  }

  double GetTime()
  {
#ifdef WIN32
    LARGE_INTEGER finishCount;
    QueryPerformanceCounter(&finishCount);

    if(s_freq.QuadPart == 0)
      QueryPerformanceFrequency(&s_freq);

    return ((double)(finishCount.QuadPart - startCount.QuadPart) /
      (double)s_freq.QuadPart);// * 1000.0;
#else
// TODO: add Linux specific things
    gettimeofday(&tv2, &tz);
  dtv.tv_sec= tv2.tv_sec -tv1.tv_sec;
  dtv.tv_usec=tv2.tv_usec-tv1.tv_usec;
  if(dtv.tv_usec<0) { dtv.tv_sec--; dtv.tv_usec+=1000000; }
  long res = dtv.tv_sec*1000+dtv.tv_usec/1000;
  return (double)res / 1000.0;
#endif
  }
};

#endif // __PERFORMANCE_H__
// - end of file ----------------------------------------------------------------------------------
