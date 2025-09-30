/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2025 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      SearchDataIterator.cpp                                       //
//                                                                         //
//  Purpose:   Source file for search data iterator classes                //
//                                                                         //
//  Author(s): Sysoyev A., Barkalov K., Sovrasov V., Zaitsev A.            //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#include "SearcDataIterator.h"

SearcDataIterator::SearcDataIterator() : pContainer(NULL), pObject(NULL)
{
}

// ------------------------------------------------------------------------------------------------
SearcDataIterator& SearcDataIterator::operator++()
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
SearcDataIterator& SearcDataIterator::operator--()
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
SearcDataIterator::operator void* () const
{
    if (pContainer && pObject)
        return pObject;
    else
        return NULL;
}

// ------------------------------------------------------------------------------------------------
SearchInterval* SearcDataIterator::operator->()
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


// - end of file ----------------------------------------------------------------------------------
