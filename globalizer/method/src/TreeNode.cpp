/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2025 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      TreeNode.cpp                                                //
//                                                                         //
//  Purpose:   Source file for search data classes                         //
//                                                                         //
//  Author(s): Sysoyev A., Barkalov K., Sovrasov V., Zaitsev A.            //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "SearchInterval.h"
#include "SearchIntervalFactory.h"
#include "TreeNode.h"

TreeNode::TreeNode(SearchInterval& p)
{
    pInterval = SearchIntervalFactory::CreateSearchInterval(p);
    Height = 1;
    pParent = pLeft = pRight = nullptr;
}

TreeNode::~TreeNode()
{
    delete pInterval;
}
