/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2016 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      property.h                                                  //
//                                                                         //
//  Purpose:   Header file for random generator class                      //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __BASE_PARAMETERS_H__
#define __BASE_PARAMETERS_H__

/**
\file baseParameters.h

\authors Ëåáåäåâ È.
\date 2015-2016
\copyright ÍÍÃÓ èì. Í.È. Ëîáà÷åâñêîãî

\brief Îáúÿâëåíèå áàçîûõ êëàññîâ äëÿ ñâîéñòâ

*/

#include "Types.h"
#include "CombinableBaseParameters.h"

/**
Áàçîâûé êëàññ ïàðàìåòðîâ
Ïðè ñîçäàíèå íàñëåäíèêà íåîáõîäèìî ïåðåîïðåäåëèòü êîíñòàíó OWNER_NAME
è â êîíñòðóêòîðå çàäàòü mOwner
(mOwner = this;)
*/
template <class Owner>
class BaseParameters : public CombinableBaseParameters
{
#undef OWNER_NAME
#define OWNER_NAME Owner

protected:
  /// Ìàññèâ ñâîéñòâ îïðåäåëåííûõ â äðóãèõ êëàññàõ
  IBaseValueClass** mOtherOptions;
  /// Ìàññèâ ñâîéñòâ îïðåäåëåííûõ â ýòîì êëàññå
  BaseProperty<Owner>** mOptions;
  /// Ìàññèâ ñâîéñòâ îïðåäåëåííûõ â ýòîì êëàññå ïðèâåäåííûå ê áàçîâîìó òèïó
  IBaseValueClass** mBaseTypeOptions;
  /// Êîëè÷åñòâî îïöèé
  int mOptionsCount;
  /// Êîëè÷åñòâî îïöèé èç äðóãèõ êëàññîâ
  int mOtherOptionsCount;
  /// Ìàêñèìàëüíûé ðàçìåð ìàññèâà ïàðàìåòðîâ
  int mOptionsSize;
  /// Ìàêñèìàëüíûé ðàçìåð ìàññèâà ïàðàìåòðîâ  îïðåäåëåííûõ â äðóãèõ êëàññàõ
  int mOtherOptionsSize;

  /// Ïå÷àòàòü ñïðàâêó
  bool mIsPrintHelp;
  /// Èìåþòñÿ àðãóìåíòû êîìàíäíîé ñòðîêè
  bool mIsHaveArguments;
  /// Âëàäåëåö ýòîãî êëàññà
  Owner* mOwner;
  /// Ïðîèíèöèàëèçèðîâàíû ëè ïàðàìåòðû
  bool mIsInit;
  /// Ïóòü ïî óìîë÷àíèþ äî êîôèãóðàöèîííîãî ôàéëà
  std::string mConfigPath;

  /// Ïðîâåðêà ïðàâèëüíîñòè ïîñëå îêîí÷àíèÿ ÷òåíèÿ ïàðàìåòðîâ
  virtual int CheckValue(int index = -1);

  /**
  Èíèöèàëèçàöèÿ ïàðàìåòðà
  \param[in] option - ïàðàìåòð êîòîðûé èíèöèàëèçèðóåì
  \param[in] sizeVal - ðàçìåð ìàññèâà çíà÷åíèé, äëÿ òèïîâ äàííûõ íå ÿâëÿþùèìèñÿ ìàñèâàìè - âñåãäà ðàâåí 1
  \param[in] name - èìÿ ñâîéñòâà
  \param[in] help - âûâîäèìàÿ íà êîíñîëü ñïðàâêà
  \param[in] com - êîðîòêàÿ ñòðîêà äëÿ çàïóñêà
  \param[in] defVal - çíà÷åíèå ïî óìîë÷àíèþ
  */
  virtual void InitializationOption(BaseProperty<Owner>* option, std::string name, std::string defVal,
    std::string com, std::string help, int sizeVal);
  /// Äîáàâëÿåò îïöèþ â îáùèé ñïèñîê
  virtual void AddOption(BaseProperty<Owner>* option);
  /// Çàäàíèå çíà÷åíèé ïî óìîë÷àíèþ áàçîâûõ ïàðàìåòðîâ
  virtual void SetBaseDefaultParameters();
  /**
  Çàäàíèå çíà÷åíèé ïî óìîë÷àíèþ äëÿ âñåõ ïàðàìåòðîâ
  Ïðèìåð:
  InitOption(èìÿ ïàðàìåòðà, çíà÷åíèå ïî óìîë÷àíèþ, "êîðîòêàÿ êîìàíäà", "ñïðàâêà ïî ïàðàìåòðó", êîë-âî ýëåìåíòîâ);
  *êîë-âî ýëåìåíòîâ äëÿ íå ìàññèâîâ âñåãäà ðàâíî 1.
  InitOption(Separator,_, "-Separator", "eparator", 1);
  */
  virtual void SetDefaultParameters() = 0;
  /// ×òåíèå ïàðàìåòðîâ èç ôàéëà ConfigPath
  //virtual void ReadConfigFile();
  /// ×òåíèå ïàðàìåòðîâ êîìàíäíîé ñòðîêè
  virtual void ReadParameters(int argc, char* argv[]);
  /// Çàäàòü ðàçäåëèòåëü ìàññèâà äëÿ âñåõ îïöèé
  void SetSeparator();


public:
  ///Ïå÷àòàòü èëè íåò ñïðàâêó ïðè ïóñòîé êîìàíäíîé ñòðîêå
  TBool<BaseParameters<Owner>> IsPrintHelpWithoutArguments;
  /// Çàïóñêàòü ëè ïðè ïóñòîé êîìàíäíîé ñòðîêå
  TBool<BaseParameters<Owner>> IsStartWithoutArguments;
  /// Ðàçäåëèòåëü ýëåìåíòîâ ìàññèâà
  TString<BaseParameters<Owner>> Separator;
  /// Ïóòü äî êîíôèã ôàéëà ïðîãðàììû
  TString<BaseParameters<Owner>> ConfigPath;
  /// Ïå÷àòü ñïðàâêè ïî ïàðàìåòðàì
  void PrintHelp();
  /// Ïå÷àòü òåêóùèõ çíà÷åíèé ïàðàìåòðîâ
  void PrintParameters();
  /// Çàïóñêàòü ëè ðàáîòó ïðîãðàììû
  bool IsStart();

  /**
  Ïðîâåðêà ïðàâèëüíîñòè ïðè èçìåíåíèå ïàðàìåòðîâ
  Ïðè ïåðåîïðåäåëåíèå íåîáõîäèìî âûçâàòü âûçâàòü ìåòîä áàçîâîãî êëàññà!
  */
  virtual int CheckValueParameters(int index = 0);

  /// Ïå÷àòü çíà÷åíèÿ ïàðàìåòðà ñ èìåíåì name
  void PrintParameter(std::string name);
  /// Çàäàíèå ïàðàìåòðó ñ èìåíåì name çíà÷åíèÿ val
  void SetVal(std::string name, std::string val);
  /// Çàäàíèå ïàðàìåòðó ñ èìåíåì name çíà÷åíèÿ val
  void SetVal(std::string name, void* val);
  /// Âîçâðàùàåò ñòðîêó ñ çíà÷åíèåì ïàðàìåòðà ñ èìåíåì name
  std::string GetStringVal(std::string name);
  /// Âîçâðàùàåò çíà÷åíèå ïàðàìåòðà ñ èìåíåì name
  void* GetVal(std::string name);

  /**
  Èíèöèàëèçàöèÿ ïàðàìåòðà
  \param[in] pt - òèï ïàðàìåòðà
  \param[in] sizeVal - ðàçìåð ìàññèâà çíà÷åíèé, äëÿ òèïîâ äàííûõ íå ÿâëÿþùèìèñÿ ìàñèâàìè - âñåãäà ðàâåí 1
  \param[in] name - èìÿ ñâîéñòâà
  \param[in] help - âûâîäèìàÿ íà êîíñîëü ñïðàâêà
  \param[in] com - êîðîòêàÿ ñòðîêà äëÿ çàïóñêà
  \param[in] defVal - çíà÷åíèå ïî óìîë÷àíèþ
  */
  virtual void AddOption(EParameterType pt, std::string name, std::string defVal,
    std::string com, std::string help, int sizeVal);

  /// Îáúåäèíÿåò âñå ïàðàìåòðû äâóõ êëàññîâ â îáîèõ
  virtual void CombineOptions(IBaseValueClass** otherOptions, int count);
  /// Âîçâðàùàåò èìåþùèåñÿ ñâîéñòâà
  virtual IBaseValueClass** GetOptions();
  /// Âîçâðàùàåò äîïîëíèòåëüíûå ñâîéñòâÿ âçÿòûå èç äðóãèõ êëàññîâ ïàðàìåòðîâ
  virtual IBaseValueClass** GetOtherOptions();
  /// Âîçâðàùàåò êîëè÷åñòâî îïöèé
  virtual int GetOptionsCount();
  /// Âîçâðàùàåò êîëè÷åñòâî äîïîëíèòåëüíûõ ñâîéñòâ âçÿòûå èç äðóãèõ êëàññîâ ïàðàìåòðîâ
  virtual int GetOtherOptionsCount();

  /// Èíèöèàëèçàöèÿ ïàðàìåòðîâ
  virtual void Init(int argc, char* argv[], bool isMPIInit = false);
  BaseParameters();
  BaseParameters(BaseParameters& _parameters);
  virtual ~BaseParameters();
  /// ßâëÿåòñÿ ëè êëàññ çàäà÷åé
  virtual bool IsProblem();
};

// ------------------------------------------------------------------------------------------------
/// Ïðîâåðêà ïðàâèëüíîñòè
template <class Owner>
int BaseParameters<Owner>::CheckValue(int index)
{

  if (Separator.ToString().length() < 1)
  {
    Separator = std::string("_");
  }

  return 0;
}

// ------------------------------------------------------------------------------------------------
/// Èíèöèàëèçàöèÿ ïàðàìåòðà
template <class Owner>
void BaseParameters<Owner>::InitializationOption(BaseProperty<Owner>* option, std::string name, std::string defVal,
  std::string com, std::string help, int sizeVal)
{
  option->InitializationParameterProperty(mOwner, &Owner::CheckValueParameters, mOptionsCount, Separator, sizeVal, name, help, com, defVal);
  AddOption(option);

  if (IsProblem())
    option->mIsEdit = true;
}

// ------------------------------------------------------------------------------------------------
/// Äîáàâëÿåò îïöèþ â îáùèé ñïèñîê
template <class Owner>
void BaseParameters<Owner>::AddOption(BaseProperty<Owner>* option)
{
  for (int i = 0; i < mOptionsCount; i++)
    if (mOptions[i]->IsNameEqual(option->GetName()))
      return;

  mOptions[mOptionsCount] = option;
  mBaseTypeOptions[mOptionsCount] = (IBaseValueClass*)option;
  mOptionsCount++;


  if (mOptionsCount >= mOptionsSize)
  {
    BaseProperty<Owner>** bufOptions = new BaseProperty<Owner>*[mOptionsSize];
    for (int i = 0; i < mOptionsSize; i++)
    {
      bufOptions[i] = mOptions[i];
    }

    delete[] mOptions;
    delete[] mBaseTypeOptions;

    mOptions = new BaseProperty<Owner>*[mOptionsSize * 2];
    mBaseTypeOptions = new IBaseValueClass*[mOptionsSize * 2];

    for (int i = 0; i < mOptionsSize * 2; i++)
    {
      mOptions[i] = 0;
      mBaseTypeOptions[i] = 0;
    }

    for (int i = 0; i < mOptionsSize; i++)
    {
      mOptions[i] = bufOptions[i];
      mBaseTypeOptions[i] = (IBaseValueClass*)mOptions[i];
    }
    mOptionsSize = mOptionsSize * 2;

    delete[] bufOptions;
  }
}

// ------------------------------------------------------------------------------------------------
/// Çàäàíèå çíà÷åíèé ïî óìîë÷àíèþ áàçîâûõ ïàðàìåòðîâ
template <class Owner>
void BaseParameters<Owner>::SetBaseDefaultParameters()
{
  int sizeVal = 1;

  Separator.InitializationParameterProperty(this, &BaseParameters::CheckValue, mOptionsCount, Separator, sizeVal,
    "Separator", "eparator", "-Separator", "_");
  AddOption((BaseProperty<Owner>*)(&Separator));

  IsPrintHelpWithoutArguments.InitializationParameterProperty(this, &BaseParameters::CheckValue, mOptionsCount, Separator, sizeVal,
    "IsPrintHelpWithoutArguments", "Is print help without console arguments", "-PHWA", "false");
  AddOption((BaseProperty<Owner>*)(&IsPrintHelpWithoutArguments));

  IsStartWithoutArguments.InitializationParameterProperty(this, &BaseParameters::CheckValue, mOptionsCount, Separator, sizeVal,
    "IsStartWithoutArguments", "Is start without console arguments", "-SWA", "true");
  AddOption((BaseProperty<Owner>*)(&IsStartWithoutArguments));

  ConfigPath.InitializationParameterProperty(this, &BaseParameters::CheckValue, mOptionsCount, Separator, sizeVal,
    "ConfigPath", "The path to the configuration file of the program", "-CP", mConfigPath);
  AddOption((BaseProperty<Owner>*)(&ConfigPath));

}

//// ------------------------------------------------------------------------------------------------
///// ×òåíèå ïàðàìåòðîâ èç ôàéëà ConfigPath
//template <class Owner>
//void BaseParameters<Owner>::ReadConfigFile()
//{
//  if (mConfigPath != "")
//    ConfigPath = mConfigPath;
//  if (ConfigPath.operator std::string() != "")
//  {
//    pugi::xml_document doc;
//    pugi::xml_parse_result result = doc.load_file(ConfigPath.ToString().c_str());
//    if (result.status != pugi::status_ok)
//      return;
//    pugi::xml_node config = doc.child("config");
//    for (pugi::xml_node iter = config.first_child(); iter != 0; iter = iter.next_sibling())
//    {
//      std::string name = iter.name();
//      for (int i = 0; i < mOptionsCount; i++)
//      {
//        if (mOptions[i]->IsNameEqual(name))
//        {
//          if (mOptions[i]->mIsEdit)
//          {
//            if (!mOptions[i]->IsFlag())
//            {
//              std::string value = iter.child_value();
//              mOptions[i]->FromString(value);
//              mOptions[i]->SetIsReadValue(true);
//              break;
//
//            }
//            else
//            {
//              mOptions[i]->FromString("1");
//              mOptions[i]->SetIsReadValue(true);
//              break;
//            }
//          }
//        }
//      }
//    }
//  }
//}

// ------------------------------------------------------------------------------------------------
/// ×òåíèå ïàðàìåòðîâ êîìàíäíîé ñòðîêè
template <class Owner>
void BaseParameters<Owner>::ReadParameters(int argc, char* argv[])
{
  for (int i = 1; i < argc; i++)
  {
    std::string argument = argv[i];
    for (int j = 0; j < mOptionsCount; j++)
    {
      if (mOptions[j]->IsNameEqual(argument))
      {
        if (!mOptions[j]->IsPreChange())
        {
          if (mOptions[j]->mIsEdit)
          {
            if (!mOptions[j]->IsFlag())
            {
              i++;
              if (i < argc)
              {
                std::string value = argv[i];
                mOptions[j]->FromString(value);
                mOptions[j]->SetIsReadValue(true);
                break;
              }
            }
            else
            {
              mOptions[j]->FromString("1");
              mOptions[j]->SetIsReadValue(true);
              break;
            }
          }
        }
      }
    }
  }
}

// ------------------------------------------------------------------------------------------------
/// Çàäàòü ðàçäåëèòåëü ìàññèâà äëÿ âñåõ îïöèé
template <class Owner>
void BaseParameters<Owner>::SetSeparator()
{
  for (int i = 0; i < mOptionsCount; i++)
  {
    mOptions[i]->SetSeparator(Separator);
  }
}

// ------------------------------------------------------------------------------------------------
/// Ïå÷àòü ñïðàâêè ïî ïàðàìåòðàì
template <class Owner>
void BaseParameters<Owner>::PrintHelp()
{
  printf("\n\nHelp:\n");
  for (int i = 0; i < mOptionsCount; i++)
  {
    if (mOptions[i]->mIsEdit)
      printf("%s\n", mOptions[i]->GetHelpString().c_str());
  }
  for (int i = 0; i < mOtherOptionsCount; i++)
  {
    if (mOtherOptions[i]->mIsEdit)
      printf("%s\n", mOtherOptions[i]->GetHelpString().c_str());
  }
}

// ------------------------------------------------------------------------------------------------
template <class Owner>
void BaseParameters<Owner>::PrintParameters()
{
  for (int i = 0; i < mOptionsCount; i++)
  {
    if (mOptions[i]->mIsEdit)
      printf("%s\n", mOptions[i]->GetCurrentStringValue().c_str());
  }
  for (int i = 0; i < mOtherOptionsCount; i++)
  {
    if (mOtherOptions[i]->mIsEdit)
      printf("%s\n", mOtherOptions[i]->GetCurrentStringValue().c_str());
  }
}

// ------------------------------------------------------------------------------------------------
/// Çàïóñêàòü ëè ðàáîòó ïðîãðàììû
template <class Owner>
bool BaseParameters<Owner>::IsStart()
{
  if ((mIsHaveArguments == false) && (IsStartWithoutArguments == false))
    return false;
  return true;
}

// ------------------------------------------------------------------------------------------------
/**
Ïðîâåðêà ïðàâèëüíîñòè ïðè èçìåíåíèå ïàðàìåòðîâ
Ïðè ïåðåîïðåäåëåíèå íåîáõîäèìî âûçâàòü âûçâàòü ìåòîä áàçîâîãî êëàññà!
*/
template <class Owner>
int BaseParameters<Owner>::CheckValueParameters(int index)
{
  if (mIsInit)
  {
    CheckValue(index);
    mIsInit = false;
    for (int i = 0; i < mOtherOptionsCount; i++)
    {
      for (int j = 0; j < mOptionsCount; j++)
      {
        if (mOtherOptions[i]->IsNameEqual(mOptions[j]->GetName()))
        {
          std::string oldBaseParameterValue = mOptions[j]->ToString();
          std::string oldOtherParameterValue = mOtherOptions[i]->ToString();

          if (oldBaseParameterValue != oldOtherParameterValue)
          {
            mOtherOptions[i]->FromString(oldBaseParameterValue);

            std::string newOtherParameterValue = mOtherOptions[i]->ToString();

            if (newOtherParameterValue == oldOtherParameterValue)
            {
              mOptions[j]->FromString(oldOtherParameterValue);
            }
          }
          break;
        }
      }
    }
    mIsInit = true;
  }
  return 0;
}

// ------------------------------------------------------------------------------------------------
/// Ïå÷àòü çíà÷åíèÿ ïàðàìåòðà ñ èìåíåì name
template <class Owner>
void BaseParameters<Owner>::PrintParameter(std::string name)
{
  for (int i = 0; i < mOptionsCount; i++)
  {
    if (mOptions[i]->IsNameEqual(name))
    {
      printf("%s\n", mOptions[i]->GetCurrentStringValue().c_str());
      break;
    }
  }
  for (int i = 0; i < mOtherOptionsCount; i++)
  {
    if (mOtherOptions[i]->IsNameEqual(name))
    {
      printf("%s\n", mOtherOptions[i]->GetCurrentStringValue().c_str());
      break;
    }
  }
}

// ------------------------------------------------------------------------------------------------
/// Çàäàíèå ïàðàìåòðó ñ èìåíåì name çíà÷åíèÿ val
template <class Owner>
void BaseParameters<Owner>::SetVal(std::string name, std::string val)
{
  for (int i = 0; i < mOptionsCount; i++)
  {
    if (mOptions[i]->IsNameEqual(name))
    {
      mOptions[i]->FromString(val);
      break;
    }
  }

  for (int i = 0; i < mOtherOptionsCount; i++)
  {
    if (mOtherOptions[i]->IsNameEqual(name))
    {
      mOtherOptions[i]->FromString(val);
      break;
    }
  }
}

// ------------------------------------------------------------------------------------------------
/// Çàäàíèå ïàðàìåòðó ñ èìåíåì name çíà÷åíèÿ val
template <class Owner>
void BaseParameters<Owner>::SetVal(std::string name, void* val)
{
  for (int i = 0; i < mOptionsCount; i++)
  {
    if (mOptions[i]->IsNameEqual(name))
    {
      mOptions[i]->SetValue(val);
      break;
    }
  }
  for (int i = 0; i < mOtherOptionsCount; i++)
  {
    if (mOtherOptions[i]->IsNameEqual(name))
    {
      mOtherOptions[i]->SetValue(val);
      break;
    }
  }
}

// ------------------------------------------------------------------------------------------------
/// Âîçâðàùàåò ñòðîêó ñ çíà÷åíèåì ïàðàìåòðà ñ èìåíåì name
template <class Owner>
std::string BaseParameters<Owner>::GetStringVal(std::string name)
{
  for (int i = 0; i < mOptionsCount; i++)
  {
    if (mOptions[i]->IsNameEqual(name))
    {
      return mOptions[i]->ToString();
    }
  }

  for (int i = 0; i < mOtherOptionsCount; i++)
  {
    if (mOtherOptions[i]->IsNameEqual(name))
    {
      return mOtherOptions[i]->ToString();
    }
  }

  return std::string("");
}

// ------------------------------------------------------------------------------------------------
/// Âîçâðàùàåò çíà÷åíèå ïàðàìåòðà ñ èìåíåì name
template <class Owner>
void* BaseParameters<Owner>::GetVal(std::string name)
{
  for (int i = 0; i < mOptionsCount; i++)
  {
    if (mOptions[i]->IsNameEqual(name))
    {
      return mOptions[i]->GetValue();
    }
  }

  for (int i = 0; i < mOtherOptionsCount; i++)
  {
    if (mOtherOptions[i]->IsNameEqual(name))
    {
      return mOtherOptions[i]->GetValue();
    }
  }

  return NULL;
}

// ------------------------------------------------------------------------------------------------
/**
Èíèöèàëèçàöèÿ ïàðàìåòðà
\param[in] pt - òèï ïàðàìåòðà
\param[in] sizeVal - ðàçìåð ìàññèâà çíà÷åíèé, äëÿ òèïîâ äàííûõ íå ÿâëÿþùèìèñÿ ìàñèâàìè - âñåãäà ðàâåí 1
\param[in] name - èìÿ ñâîéñòâà
\param[in] help - âûâîäèìàÿ íà êîíñîëü ñïðàâêà
\param[in] com - êîðîòêàÿ ñòðîêà äëÿ çàïóñêà
\param[in] defVal - çíà÷åíèå ïî óìîë÷àíèþ
*/
template <class Owner>
void BaseParameters<Owner>::AddOption(EParameterType pt, std::string name, std::string defVal,
  std::string com, std::string help, int sizeVal)
{
  BaseProperty<Owner>* option = 0;

  if (pt == Pint)
    option = new TInt<Owner>();
  else if (pt == Pdouble)
    option = new TDouble<Owner>();
  else if (pt == Pstring)
    option = new TString<Owner>();
  else if (pt == PETypeMethod)
    option = new TETypeMethod<Owner>();
  else if (pt == PETypeProcess)
    option = new TETypeProcess<Owner>();
  else if (pt == PETypeCalculation)
    option = new TETypeCalculation<Owner>();
  else if (pt == PELocalMethodScheme)
    option = new TELocalMethodScheme<Owner>();
  else if (pt == PESeparableMethodType)
    option = new TESeparableMethodType<Owner>();
  else if (pt == PEStopCondition)
    option = new TEStopCondition<Owner>();
  else if (pt == Pbool)
    option = new TBool<Owner>();
  else if (pt == Pints)
    option = new TInts<Owner>();
  else if (pt == Pdoubles)
    option = new TDoubles<Owner>();
  else if (pt == PEMapType)
    option = new TEMapType<Owner>();

  option->InitializationParameterProperty(mOwner, 0, mOptionsCount, Separator, sizeVal,
    name, help, com, defVal);
  AddOption(option);
}

// ------------------------------------------------------------------------------------------------
/// Îáúåäèíÿåò âñå ïàðàìåòðû äâóõ êëàññîâ â îáîèõ
template <class Owner>
void BaseParameters<Owner>::CombineOptions(IBaseValueClass** otherOptions, int count)
{
  int newOtherOptionsCount = mOtherOptionsCount + count;
  if (newOtherOptionsCount >= mOtherOptionsSize)
  {
    IBaseValueClass** bufOptions = new IBaseValueClass*[mOtherOptionsSize];
    for (int i = 0; i < mOtherOptionsSize; i++)
    {
      bufOptions[i] = mOtherOptions[i];
    }

    delete[] mOtherOptions;

    mOtherOptions = new IBaseValueClass*[mOtherOptionsSize * 2];

    for (int i = 0; i < mOtherOptionsSize * 2; i++)
    {
      mOtherOptions[i] = 0;
    }

    for (int i = 0; i < mOtherOptionsSize; i++)
    {
      mOtherOptions[i] = bufOptions[i];
    }
    mOtherOptionsSize = mOtherOptionsSize * 2;

    delete[] bufOptions;
  }

  for (int i = 0; i < count; i++)
  {
    mOtherOptions[mOtherOptionsCount] = otherOptions[i];

    for (int j = 0; j < mOptionsCount; j++)
    {
      if (mOtherOptions[mOtherOptionsCount]->IsNameEqual(mOptions[j]->GetName()))
      {
        std::string oldBaseParameterValue = mOptions[j]->ToString();
        std::string oldOtherParameterValue = mOtherOptions[mOtherOptionsCount]->ToString();
        if (!mOptions[j]->IsPreChange() && !mOtherOptions[mOtherOptionsCount]->IsPreChange())
        {
          if (oldBaseParameterValue != oldOtherParameterValue)
          {
            mOptions[j]->FromString(oldOtherParameterValue);
            std::string newBaseParameterValue = mOptions[j]->ToString();

            if (newBaseParameterValue == oldBaseParameterValue)
            {
              mOtherOptions[mOtherOptionsCount]->FromString(oldBaseParameterValue);
            }
          }
        }
        break;
      }
    }

    mOtherOptionsCount++;
  }
}

// ------------------------------------------------------------------------------------------------
/// Âîçâðàùàåò èìåþùèåñÿ ñâîéñòâà
template <class Owner>
IBaseValueClass** BaseParameters<Owner>::GetOptions()
{
  return mBaseTypeOptions;
}

// ------------------------------------------------------------------------------------------------
/// Âîçâðàùàåò äîïîëíèòåëüíûå ñâîéñòâà âçÿòûå èç äðóãèõ êëàññîâ ïàðàìåòðîâ
template <class Owner>
IBaseValueClass** BaseParameters<Owner>::GetOtherOptions()
{
  return mOtherOptions;
}

// ------------------------------------------------------------------------------------------------
/// Âîçâðàùàåò êîëè÷åñòâî îïöèé
template <class Owner>
int BaseParameters<Owner>::GetOptionsCount()
{
  return mOptionsCount;
}

// ------------------------------------------------------------------------------------------------
/// Âîçâðàùàåò êîëè÷åñòâî äîïîëíèòåëüíûõ ñâîéñòâ âçÿòûå èç äðóãèõ êëàññîâ ïàðàìåòðîâ
template <class Owner>
int BaseParameters<Owner>::GetOtherOptionsCount()
{
  return mOtherOptionsCount;
}

// ------------------------------------------------------------------------------------------------
template <class Owner>
void BaseParameters<Owner>::Init(int argc, char* argv[], bool isMPIInit)
{
  mIsMPIInit = isMPIInit;
  mIsInit = false;
  mOptionsCount = 0;
  mOtherOptionsCount = 0;
  mOptions = 0;

  mOptionsSize = 200;
  mOtherOptionsSize = 200;

  mOptions = new BaseProperty<Owner>*[mOptionsSize];
  mOtherOptions = new IBaseValueClass*[mOtherOptionsSize];
  mBaseTypeOptions = new IBaseValueClass*[mOptionsSize];
  for (int i = 0; i < mOptionsSize; i++)
  {
    mOptions[i] = 0;
    mBaseTypeOptions[i] = 0;
  }

  for (int i = 0; i < mOtherOptionsSize; i++)
  {
    mOtherOptions[i] = 0;
  }

  // Îïðåäåëÿåì åñòü ëè ïàðàìåòðû êîíñîëè

  if (argc > 0)
  {
    mArgumentCount = argc;
    mAargumentValue = argv;
  }

  if (mArgumentCount <= 1)
    mIsHaveArguments = false;
  else
    mIsHaveArguments = true;

  // Èíèöèàëèçàöèÿ áàçîâûõ ïàðàìåòðîâ ïî óìîë÷àíèþ
  SetBaseDefaultParameters();
  // Èíèöèàëèçàöèÿ ðàáî÷èõ ïàðàìåòðîâ
  SetDefaultParameters();
  // Çàäàòü ðàçäåëèòåëü äëÿ ìàññèâîâ
  SetSeparator();

  // Îïðåäåëÿåì ïàðàìåòðû èç êîíñîëè
  ReadParameters(mArgumentCount, mAargumentValue);


  // Îïðåäåëÿåì ïàðàìåòðû èç ôàéëîâ
  //ReadConfigFile();
  // Ïðîâåðêà ïàðàìåòðîâ
  CheckValue();

  if ((mIsHaveArguments == false) && (IsPrintHelpWithoutArguments == false))
    mIsPrintHelp = false;
  else if (mIsHaveArguments == true)
    mIsPrintHelp = false;
  else
    mIsPrintHelp = true;
  mIsInit = true;

  CheckValueParameters();
}

// ------------------------------------------------------------------------------------------------
template <class Owner>
BaseParameters<Owner>::BaseParameters()
{
  mConfigPath = "";

  mOwner = 0;
  mIsInit = false;
}

// ------------------------------------------------------------------------------------------------
template <class Owner>
BaseParameters<Owner>::BaseParameters(BaseParameters& _parameters) : CombinableBaseParameters()
{
  mIsInit = false;
  mOptionsCount = 0;
  mOtherOptionsCount = 0;
  mIsMPIInit = _parameters.mIsMPIInit;
  mOptionsSize = _parameters.mOptionsSize;
  mOtherOptionsSize = _parameters.mOtherOptionsSize;
  mOwner = 0;
  mIsHaveArguments = _parameters.mIsHaveArguments;
  IsPrintHelpWithoutArguments = _parameters.IsPrintHelpWithoutArguments;
  IsStartWithoutArguments = _parameters.IsStartWithoutArguments;
  mIsPrintHelp = _parameters.mIsPrintHelp;
  mOtherOptionsCount = _parameters.mOtherOptionsCount;

  mOptions = new BaseProperty<Owner>*[mOptionsSize];
  mOtherOptions = new IBaseValueClass*[mOtherOptionsSize];
  mBaseTypeOptions = new IBaseValueClass*[mOptionsSize];
  for (int i = 0; i < mOptionsSize; i++)
  {
    mOptions[i] = 0;

    mBaseTypeOptions[i] = 0;
  }

  for (int i = 0; i < mOtherOptionsSize; i++)
  {
    mOtherOptions[i] = 0;
  }

  // Èíèöèàëèçàöèÿ áàçîâûõ ïàðàìåòðîâ ïî óìîë÷àíèþ
  SetBaseDefaultParameters();

  for (int i = 0; i < mOptionsCount; i++)
  {
    *mOptions[i] = *_parameters.mOptions[i];
  }

  for (int i = 0; i < mOtherOptionsCount; i++)
  {
    mOtherOptions[i] = _parameters.mOtherOptions[i];
  }

}

// ------------------------------------------------------------------------------------------------
template <class Owner>
BaseParameters<Owner>::~BaseParameters()
{
  delete[] mOptions;
  mOptions = 0;
}

/// ßâëÿåòñÿ ëè êëàññ çàäà÷åé
template <class Owner>
bool BaseParameters<Owner>::IsProblem()
{
  return true;
}

#endif //__BASE_PARAMETERS_H__
// - end of file ----------------------------------------------------------------------------------
