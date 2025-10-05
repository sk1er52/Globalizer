/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2025 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      SearchInterval.cpp                                           //
//                                                                         //
//  Purpose:   Source file for search interval classes                     //
//                                                                         //
//  Author(s): Sysoyev A., Barkalov K., Sovrasov V., Zaitsev A.            //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "SearchInterval.h"
#include "TrialFactory.h"

SearchInterval::SearchInterval()
{
    LeftPoint = nullptr;
    RightPoint = nullptr;

    R = locR = 0.0;
    ind = 0;
    // Выделение памяти происходит в момент добавления в матрицу
    K = 0;
    delta = 0.0;

    queueElementa = nullptr;
    treeNode = nullptr;
}

// ------------------------------------------------------------------------------------------------
SearchInterval::SearchInterval(const SearchInterval& p)
{
    LeftPoint = p.LeftPoint;
    RightPoint = p.RightPoint;

    ind = p.ind;
    K = p.K;
    R = p.R;
    locR = p.locR;
    delta = p.delta;

    //МКО

    queueElementa = p.queueElementa;
    treeNode = p.treeNode;
}

// ------------------------------------------------------------------------------------------------
SearchInterval::~SearchInterval()
{
}

// ------------------------------------------------------------------------------------------------
void SearchInterval::CreatePoint()
{
    LeftPoint = TrialFactory::CreateTrial();
    RightPoint = TrialFactory::CreateTrial();
    RightPoint->SetX(1);
}


// - end of file ----------------------------------------------------------------------------------
