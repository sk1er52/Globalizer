/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2016 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      ParameterProperty.h                                                  //
//                                                                         //
//  Purpose:   Header file for random generator class                      //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __PARAMETER_PROPERTY_H__
#define __PARAMETER_PROPERTY_H__

#include <cctype>

/**
\file ParameterProperty.h

\authors Лебедев И.
\date 2015-2016
\copyright ННГУ им. Н.И. Лобачевского

\brief Объявление базоых классов для свойств

*/

#include "TypedProperty.h"

//

/* ======================================================================== *\
**  Объявление классов                                                      **
\* ======================================================================== */

/**
Базовый класс для свойств параметров, имеет полную реализацию
*/
template <class Type, class Owner>
class ParameterProperty : public TypedProperty<Type, Owner>
{
protected:
  typedef typename BaseProperty<Owner>::tCheckValue tCheckValue;
  typedef typename TypedProperty<Type, Owner>::tGetter tGetter;
  typedef typename TypedProperty<Type, Owner>::tSetter tSetter;
  /// Разделитель массива
  std::string mSeparator;
  /// Размер массива
  int mSize;
  /// Имя свойства
  std::string mName;
  /// Текст справки
  std::string mHelp;
  /// Короткое имя для домандной строки
  std::string mLink;

  /// Полное консольное имя из имени
  std::string GetFullLink();

public:

  /// Перегрузка оператора копрования двух объектов
  virtual void operator = (ParameterProperty<Type, Owner>&  data);
  /// Перегрузка оператора копрования для указателя
  virtual void operator = (ParameterProperty<Type, Owner>* data);

  /// Создает объект data
  virtual void Clone(BaseProperty<Owner>** data);

  /// Задает разделитель элементов массива
  virtual void SetSeparator(std::string separator);
  /// Возвращает разделитель элементов массива
  virtual std::string GetSeparator();

  /// Задает размер массива элементов
  virtual void SetSize(int size);
  /// Возвращает Размер массива элементов
  virtual int GetSize();

  /// Имя свойства
  virtual std::string GetName();
  /// Текст справки
  virtual std::string GetHelp();
  /// Короткое имя для домандной строки
  virtual std::string GetLink();

  /// Имя свойства
  virtual void SetName(std::string name);
  /// Текст справки
  virtual void SetHelp(std::string help);
  /// Короткое имя для домандной строки
  virtual void SetLink(std::string link);
  /// Возвращает справку по свойству
  virtual std::string GetHelpString();
  /// Возвращает текущего состояния параметра
  virtual std::string GetCurrentStringValue();

  /// Приводит к строке
  virtual std::string ToString();
  /// Получение значения из строки
  virtual void FromString(std::string val);

  /// Совпадает ли имя параметра со введенной строкой
  virtual bool IsNameEqual(std::string name);

  /** Инициализация свойства
  \param[in] owner - класс владелец свойства
  \param[in] checkMethod - метод проверки правильности введенных данных
  \param[in] index - номер свойства
  \param[in] separator - разделитель элементов массива
  \param[in] size - размер массива значений, для типов данных не являющимися масивами - всегда равен 1
  \param[in] name - имя свойства
  \param[in] help - выводимая на консоль справка
  \param[in] link - короткая строка для запуска
  \param[in] defValue - значение по умолчанию
  */
  virtual void InitializationParameterProperty(Owner * owner,
    tCheckValue checkMethod, int index, std::string separator, int size, std::string name,
    std::string help, std::string link, std::string defValue);

  /** Инициализация свойства
  \param[in] owner - класс владелец свойства
  \param[in] getMethod - геттер
  \param[in] setMethod - сеттер
  \param[in] checkMethod - метод проверки правильности введенных данных
  \param[in] index - номер свойства
  \param[in] separator - разделитель элементов массива
  \param[in] size - размер массива значений, для типов данных не являющимися масивами - всегда равен 1
  \param[in] name - имя свойства
  \param[in] help - выводимая на консоль справка
  \param[in] link - короткая строка для запуска
  \param[in] defValue - значение по умолчанию
  */
  virtual void InitializationParameterProperty(Owner * owner, tGetter getMethod, tSetter setMethod,
    tCheckValue checkMethod, int index, std::string separator, int size, std::string name,
    std::string help, std::string link, std::string defValue);

  ParameterProperty();
  ParameterProperty(Type value);
  ParameterProperty(Owner * owner, tGetter getMethod, tSetter setMethod,
    tCheckValue checkMethod, Type value);




};


/* ======================================================================== *\
**  Реализация методов класса     ParameterProperty                                 **
\* ======================================================================== */

// ------------------------------------------------------------------------------------------------
/// Перегрузка оператора копрования двух объектов
template <class Type, class Owner>
void ParameterProperty<Type, Owner>::operator = (ParameterProperty<Type, Owner>&  data)
{
  mSeparator = data.mSeparator;
  mSize = data.mSize;
  mName = data.mName;
  mHelp = data.mHelp;
  mLink = data.mLink;
  this->mIsFlag = data.mIsFlag;
  TypedProperty<Type, Owner>::operator=(data);
}

// ------------------------------------------------------------------------------------------------
/// Перегрузка оператора копрования для указателя
template <class Type, class Owner>
void ParameterProperty<Type, Owner>::operator = (ParameterProperty<Type, Owner>* data)
{
  mSeparator = data->mSeparator;
  mSize = data->mSize;
  mName = data->mName;
  mHelp = data->mHelp;
  mLink = data->mLink;
  this->mIsFlag = data->mIsFlag;
  TypedProperty<Type, Owner>::operator=(data);
}

// ------------------------------------------------------------------------------------------------
/// Создает объект data
template <class Type, class Owner>
void ParameterProperty<Type, Owner>::Clone(BaseProperty<Owner>** data)
{
  *data = new ParameterProperty<Type, Owner>();

  *((ParameterProperty<Type, Owner>*)*data) = *this;
}

// ------------------------------------------------------------------------------------------------
/// Задает разделитель элементов массива
template <class Type, class Owner>
void ParameterProperty<Type, Owner>::SetSeparator(std::string separator)
{
  mSeparator = separator;
}

// ------------------------------------------------------------------------------------------------
/// Возвращает разделитель элементов массива
template <class Type, class Owner>
std::string ParameterProperty<Type, Owner>::GetSeparator()
{
  return mSeparator;
}

// ------------------------------------------------------------------------------------------------
/// Задает размер массива элементов
template <class Type, class Owner>
void ParameterProperty<Type, Owner>::SetSize(int size)
{
  //if (size > 0)
  //  mSize = size;
}

// ------------------------------------------------------------------------------------------------
/// Возвращает Размер массива элементов
template <class Type, class Owner>
int ParameterProperty<Type, Owner>::GetSize()
{
  return mSize;
}

// ------------------------------------------------------------------------------------------------
/// Имя свойства
template <class Type, class Owner>
std::string ParameterProperty<Type, Owner>::GetName()
{
  return mName;
}

// ------------------------------------------------------------------------------------------------
/// Текст справки
template <class Type, class Owner>
std::string ParameterProperty<Type, Owner>::GetHelp()
{
  return mHelp;
}

// ------------------------------------------------------------------------------------------------
/// Короткое имя для домандной строки
template <class Type, class Owner>
std::string ParameterProperty<Type, Owner>::GetLink()
{
  return mLink;
}

// ------------------------------------------------------------------------------------------------
/// Имя свойства
template <class Type, class Owner>
void ParameterProperty<Type, Owner>::SetName(std::string name)
{
  mName = name;
}

// ------------------------------------------------------------------------------------------------
/// Текст справки
template <class Type, class Owner>
void ParameterProperty<Type, Owner>::SetHelp(std::string help)
{
  mHelp = help;
}

// ------------------------------------------------------------------------------------------------
/// Короткое имя для домандной строки
template <class Type, class Owner>
void ParameterProperty<Type, Owner>::SetLink(std::string link)
{
  mLink = link;
}

// ------------------------------------------------------------------------------------------------
/// Возвращает справку по свойству
template <class Type, class Owner>
std::string ParameterProperty<Type, Owner>::GetHelpString()
{
  std::string result = "";
  result = result + GetName() + " (" + GetLink() + ") - \'" + GetHelp() + "\' default:\t" + ToString();
  return result;
}

// ------------------------------------------------------------------------------------------------
/// Возвращает текущего состояния параметра
template <class Type, class Owner>
std::string ParameterProperty<Type, Owner>::GetCurrentStringValue()
{
  std::string result = "";
  result = result + GetName() + " = " + ToString();
  return result;
}

// ------------------------------------------------------------------------------------------------
/// Приводит к строке
template <class Type, class Owner>
std::string ParameterProperty<Type, Owner>::ToString()
{
  return "";
}

// ------------------------------------------------------------------------------------------------
/// Получение значения из строки
template <class Type, class Owner>
void ParameterProperty<Type, Owner>::FromString(std::string val)
{}

// ------------------------------------------------------------------------------------------------
/// Полное консольное имя из имени
template <class Type, class Owner>
std::string ParameterProperty<Type, Owner>::GetFullLink()
{
  return "-" + GetName();
}

// ------------------------------------------------------------------------------------------------
/// Совпадает ли имя параметра со введенной строкой
template <class Type, class Owner>
bool ParameterProperty<Type, Owner>::IsNameEqual(std::string name)
{
  if ((name.length() == GetFullLink().length()) || (name.length() == GetLink().length()) || (name.length() == GetName().length()))
  {
    std::string lName = name;
    std::string LGetFullLink = GetFullLink();
    std::string LGetLink = GetLink();
    std::string LGetName = GetName();

    for (int i = 0; i < lName.length(); i++) 
      lName[i] = std::tolower(lName[i]);

    for (int i = 0; i < LGetFullLink.length(); i++)
      LGetFullLink[i] = std::tolower(LGetFullLink[i]);

    for (int i = 0; i < LGetLink.length(); i++)
      LGetLink[i] = std::tolower(LGetLink[i]);

    for (int i = 0; i < LGetName.length(); i++)
      LGetName[i] = std::tolower(LGetName[i]);


    if ((lName == LGetFullLink) || (lName == LGetLink) || (lName == LGetName))
      return true;  
    else
      return false;
  }
  else
    return false;
}

// ------------------------------------------------------------------------------------------------
template <class Type, class Owner>
void ParameterProperty<Type, Owner>::InitializationParameterProperty(Owner * owner,
  tCheckValue checkMethod, int index, std::string separator, int size, std::string name,
  std::string help, std::string link, std::string defValue)
{
  TypedProperty<Type, Owner>::Init(owner, 0, 0, checkMethod);

  this->SetIndex(index);
  SetSeparator(separator);
  SetSize(size);
  SetName(name);
  SetHelp(help);
  SetLink(link);

  FromString(defValue);
  this->mIsChange = false;
}

// ------------------------------------------------------------------------------------------------
template <class Type, class Owner>
void ParameterProperty<Type, Owner>::InitializationParameterProperty(Owner * owner, tGetter getMethod, tSetter setMethod,
  tCheckValue checkMethod, int index, std::string separator, int size, std::string name,
  std::string help, std::string link, std::string defValue)
{
  TypedProperty<Type, Owner>::Init(owner, getMethod, setMethod, checkMethod);

  this->SetIndex(index);
  SetSeparator(separator);
  SetSize(size);
  SetName(name);
  SetHelp(help);
  SetLink(link);

  FromString(defValue);
  this->mIsChange = false;
}

// ------------------------------------------------------------------------------------------------
template <class Type, class Owner>
ParameterProperty<Type, Owner>::ParameterProperty() : mSize(1), TypedProperty<Type, Owner>()
{
  this->mIsEdit = true;

  mSeparator = "_";
  mName = "";
  mHelp = "";
  mLink = "";
}

// ------------------------------------------------------------------------------------------------
template <class Type, class Owner>
ParameterProperty<Type, Owner>::ParameterProperty(Type value) : mSize(1), TypedProperty<Type, Owner>(value)
{
  this->mIsEdit = true;

  mSeparator = "_";
  mName = "";
  mHelp = "";
  mLink = "";
}

// ------------------------------------------------------------------------------------------------
template <class Type, class Owner>
ParameterProperty<Type, Owner>::ParameterProperty(Owner * owner, tGetter getMethod,
  tSetter setMethod, tCheckValue checkMethod, Type value) :
  mSize(1), TypedProperty<Type, Owner>(this->mOwner, this->mGetter, this->mSetter, this->mCheckValue, value)
{
  this->mIsEdit = true;

  mSeparator = "_";
  mName = "";
  mHelp = "";
  mLink = "";
}

#endif //__PARAMETER_PROPERTY_H__
// - end of file ----------------------------------------------------------------------------------
