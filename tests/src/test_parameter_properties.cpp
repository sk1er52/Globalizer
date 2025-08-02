//#include "types.h"
//#include "CombinableBaseParameters.h"

#include "Parameters.h"

#include <gtest/gtest.h>
#include <string>
#include <cstdlib>

using namespace std;

int PropertiesIndex = 11;

class PropertiesTest : public ::testing::Test
{

public:

  bool boolVal;

  bool GetBool() const
  {
    return boolVal;
  }
  void SetBool(bool val)
  {
    boolVal = val;
  }

  virtual int CheckValue(int index = -1)
  {
    boolVal = !boolVal;

    return 0;
  }

  virtual void TestBody()
  {  }

  PropertiesTest()
  {
    boolVal = true;
  }
};

/**
Проверка класса TBool
*/

TEST(Properties_TBool, can_create_default_Bool)
{
  ASSERT_NO_THROW(TBool<PropertiesTest> a);
}

TEST(Properties_TBool, can_create_Bool)
{
  bool val = true;
  ASSERT_NO_THROW(TBool<PropertiesTest> a(val));
}

TEST(Properties_TBool, is_init_Bool_value)
{
  bool val = true;
  TBool<PropertiesTest> b(val);
  ASSERT_EQ(val, (bool)b);
}

TEST(Properties_TBool, is_getter_and_setter_working_Bool)
{
  bool val = true;
  TBool<PropertiesTest> b(!val);

  b = val;

  ASSERT_EQ(val, (bool)b);
}

TEST(Properties_TBool, is_SetIndex_and_GetIndex_working_Bool)
{
  int val = PropertiesIndex++;
  TBool<PropertiesTest> b;

  b.SetIndex(val);

  ASSERT_EQ(val, b.GetIndex());
}

TEST(Properties_TBool, is_Bool_GetData_working)
{
  bool val = true;
  TBool<PropertiesTest> b(val);
  ASSERT_EQ(val, b.GetData());
}

TEST(Properties_TBool, is_Bool_GetValue_working)
{
  bool val = true;
  TBool<PropertiesTest> b(val);
  ASSERT_EQ(val, *((bool*) b.GetValue()));
}

TEST(Properties_TBool, is_Bool_Clone_working)
{
  bool val = true;
  TBool<PropertiesTest> b(val);
  TBool<PropertiesTest>* c;

  b.Clone((BaseProperty<PropertiesTest>**)&c);

  ASSERT_EQ(val, c->GetData());
}

TEST(Properties_TBool, is_Bool_GetIsChange_working)
{
  bool val = true;
  TBool<PropertiesTest> b(!val);

  ASSERT_EQ(false, b.GetIsChange());

  b = val;

  ASSERT_EQ(true, b.GetIsChange());
}

TEST(Properties_TBool, is_Bool_Copy_working)
{
  bool val = true;
  TBool<PropertiesTest> b(!val);
  TBool<PropertiesTest> c(val);

  b.Copy((void*)&c);

  ASSERT_EQ(val, b.GetData());
}

TEST(Properties_TBool, is_Bool_GetCurrentStringValue_working)
{
  bool val = true;
  string name = "n";
  string result = "n = true";
  TBool<PropertiesTest> b(val);

  b.SetName(name);

  ASSERT_EQ(result, b.GetCurrentStringValue());
}

TEST(Properties_TBool, is_Bool_SetName_and_GetName_working)
{
  string name = "n";
  TBool<PropertiesTest> b;

  ASSERT_NO_THROW(b.SetName(name));

  ASSERT_EQ(name, b.GetName());
}

TEST(Properties_TBool, is_Bool_IsNameEqual_working)
{
  string name = "n";
  TBool<PropertiesTest> b;

  ASSERT_NO_THROW(b.SetName(name));

  ASSERT_EQ(true, b.IsNameEqual(name));
}

TEST(Properties_TBool, is_Bool_IsFlag_working)
{
  TBool<PropertiesTest> b;

  ASSERT_EQ(false, b.IsFlag());
}

TEST(Properties_TBool, is_Bool_SetIsReadValue_and_GetIsReadValue_working)
{
  TBool<PropertiesTest> b;

  ASSERT_NO_THROW(b.SetIsReadValue(true));

  ASSERT_EQ(true, b.GetIsReadValue());
}

TEST(Properties_TBool, is_Bool_SetIsPreChange_and_IsPreChange_working)
{
  TBool<PropertiesTest> b;

  ASSERT_NO_THROW(b.SetIsPreChange(true));

  ASSERT_EQ(true, b.IsPreChange());
}

TEST(Properties_TBool, is_Bool_SetHelp_and_GetHelp_working)
{
  string help = "n";
  TBool<PropertiesTest> b;

  ASSERT_NO_THROW(b.SetHelp(help));

  ASSERT_EQ(help, b.GetHelp());
}

TEST(Properties_TBool, is_Bool_SetLink_and_GetLink_working)
{
  string link = "n";
  TBool<PropertiesTest> b;

  ASSERT_NO_THROW(b.SetLink(link));

  ASSERT_EQ(link, b.GetLink());
}

TEST(Properties_TBool, is_Bool_GetHelpString_working)
{
  string result = "b (-b) - 'This is B' default:\ttrue";
  string link = "-b";
  string help = "This is B";
  string name = "b";
  bool val = true;
  TBool<PropertiesTest> b(val);

  b.SetLink(link);
  b.SetHelp(help);
  b.SetName(name);

  ASSERT_EQ(result, b.GetHelpString());
}

TEST(Properties_TBool, is_init_function_Bool)
{
  PropertiesTest a;
  TBool<PropertiesTest> b;

  ASSERT_NO_THROW(b.Init(&a, &PropertiesTest::GetBool, &PropertiesTest::SetBool,
    &PropertiesTest::CheckValue));
}


TEST(Properties_TBool, is_InitializationParameterProperty_function_Bool)
{
  int index = 12;
  string result = "b (-b) - 'This is B' default:\ttrue";
  string link = "-b";
  string help = "This is B";
  string name = "b";
  string sep = "_";
  string defVal = "true";
  bool val = true;

  PropertiesTest a;
  TBool<PropertiesTest> b;

  ASSERT_NO_THROW(b.InitializationParameterProperty(&a, &PropertiesTest::CheckValue, index,
    sep, 1, name, help, link, defVal));

  ASSERT_EQ(val, b);
  ASSERT_EQ(link, b.GetLink());
  ASSERT_EQ(help, b.GetHelp());
  ASSERT_EQ(name, b.GetName());
  ASSERT_EQ(result, b.GetHelpString());
}

TEST(Properties_TBool, is_owner_getter_working_Bool)
{
  bool val = true;
  PropertiesTest a;
  TBool<PropertiesTest> b(!val);

  b.Init(&a, &PropertiesTest::GetBool, &PropertiesTest::SetBool, 0);
  a.boolVal = val;

  ASSERT_EQ(val, (bool)b);
}

TEST(Properties_TBool, is_owner_setter_working_Bool)
{
  bool val = true;
  PropertiesTest a;
  TBool<PropertiesTest> b(!val);

  b.Init(&a, &PropertiesTest::GetBool, &PropertiesTest::SetBool, 0);

  a.boolVal = !val;

  b = val;

  ASSERT_EQ(val, a.boolVal);
}

TEST(Properties_TBool, is_CheckValue_working_Bool)
{
  bool val = true;
  PropertiesTest a;
  TBool<PropertiesTest> b;

  b.Init(&a, &PropertiesTest::GetBool, &PropertiesTest::SetBool,
    &PropertiesTest::CheckValue);

  b = !val;

  ASSERT_EQ(val, b);
}

TEST(Properties_TBool, is_GetAvailableData_working_Bool)
{
  bool val = true;
  PropertiesTest a;
  TBool<PropertiesTest> b(val);

  b.Init(&a, &PropertiesTest::GetBool, &PropertiesTest::SetBool, 0);

  b = !val;

  ASSERT_EQ(val, b.GetAvailableData());
}

TEST(Properties_TBool, is_GetGetter_and_SetGetter_working_Bool)
{
  TBool<PropertiesTest> b;

  ASSERT_NO_THROW(b.SetGetter(&PropertiesTest::GetBool));

  ASSERT_EQ(&PropertiesTest::GetBool, b.GetGetter());
}

TEST(Properties_TBool, is_GetIsHaveValue_working_Bool)
{
  bool val = true;
  TBool<PropertiesTest> b(val);

  ASSERT_EQ(true, b.GetIsHaveValue());

  ASSERT_NO_THROW(b.SetGetter(&PropertiesTest::GetBool));

  ASSERT_EQ(false, b.GetIsHaveValue());

  ASSERT_NO_THROW(b.SetGetter(0));

  ASSERT_EQ(true, b.GetIsHaveValue());

  ASSERT_NO_THROW(b.SetSetter(&PropertiesTest::SetBool));

  ASSERT_EQ(false, b.GetIsHaveValue());
}

TEST(Properties_TBool, is_SetIsHaveValue_working_Bool)
{
  bool val = true;
  PropertiesTest a;
  TBool<PropertiesTest> b(val);

  b.Init(&a, &PropertiesTest::GetBool, &PropertiesTest::SetBool,
    &PropertiesTest::CheckValue);

  a.boolVal = !val;

  ASSERT_EQ(false, b.GetIsHaveValue());

  ASSERT_EQ(!val, b);

  b.SetIsHaveValue(true);

  ASSERT_EQ(true, b.GetIsHaveValue());

  ASSERT_EQ(val, b);
}

TEST(Properties_TBool, is_GetSetter_and_SetSetter_working_Bool)
{
  TBool<PropertiesTest> b;

  ASSERT_NO_THROW(b.SetSetter(&PropertiesTest::SetBool));

  ASSERT_EQ(&PropertiesTest::SetBool, b.GetSetter());
}

TEST(Properties_TBool, is_SetCheckValue_and_GetCheckValue_working_Bool)
{
  TBool<PropertiesTest> b;

  ASSERT_NO_THROW(b.SetCheckValue(&PropertiesTest::CheckValue));

  ASSERT_EQ(&PropertiesTest::CheckValue, b.GetCheckValue());
}

TEST(Properties_TBool, is_CheckValue__working_Bool)
{
  PropertiesTest a;
  TBool<PropertiesTest> b;
  a.boolVal = false;

  b.Init(&a, &PropertiesTest::GetBool, &PropertiesTest::SetBool, &PropertiesTest::CheckValue);

  ASSERT_NO_THROW(b.CheckValue());

  ASSERT_EQ(true, a.boolVal);
}


TEST(Properties_TBool, is_Bool_ToString_working)
{
  string result = "true";
  bool val = true;
  TBool<PropertiesTest> b(val);

  ASSERT_EQ(result, b.ToString());
}

TEST(Properties_TBool, is_Bool_FromString_working)
{
  bool val = true;
  string sVal = "true";
  TBool<PropertiesTest> b(!val);

  ASSERT_NO_THROW(b.FromString(sVal));

  ASSERT_EQ(val, b);
}

TEST(Properties_TBool, is_Bool_operator_FromString_working)
{
  bool val = true;
  string sVal = "true";
  TBool<PropertiesTest> b(!val);

  b = sVal;

  ASSERT_EQ(val, b);
}

TEST(Properties_TBool, is_Bool_operator_ToString_working)
{
  string result = "true";
  bool val = true;
  TBool<PropertiesTest> b(val);

  ASSERT_EQ(result, (string)b);
}

/**
Проверка класса TFlag
*/

TEST(Properties_TFlag, can_create_default_Flag)
{
  ASSERT_NO_THROW(TFlag<PropertiesTest> a);
}

TEST(Properties_TFlag, can_create_Flag)
{
  bool val = true;
  ASSERT_NO_THROW(TFlag<PropertiesTest> a(val));
}

TEST(Properties_TFlag, is_init_Flag_value)
{
  bool val = true;
  TFlag<PropertiesTest> b(val);
  ASSERT_EQ(val, (bool)b);
}

TEST(Properties_TFlag, is_Flag_IsFlag_working)
{
  TFlag<PropertiesTest> b;

  ASSERT_EQ(true, b.IsFlag());
}

TEST(Properties_TFlag, is_Flag_operator_FromString_working)
{
  bool val = true;
  string sVal = "true";
  TFlag<PropertiesTest> b(!val);

  b = sVal;

  ASSERT_EQ(val, b);
}

TEST(Properties_TFlag, is_Flag_operator_ToString_working)
{
  string result = "true";
  bool val = true;
  TFlag<PropertiesTest> b(val);

  ASSERT_EQ(result, (string)b);
}

/**
Проверка класса TInt
*/

TEST(Properties_TInt, can_create_default_Int)
{
  ASSERT_NO_THROW(TInt<PropertiesTest> a);
}

TEST(Properties_TInt, can_create_Int)
{
  int val = 42;
  ASSERT_NO_THROW(TInt<PropertiesTest> a(val));
}

TEST(Properties_TInt, is_init_Int_value)
{
  int val = 42;
  TInt<PropertiesTest> b(val);
  ASSERT_EQ(val, (int)b);
}

TEST(Properties_TInt, is_Int_operator_FromString_working)
{
  int val = 17;
  string sVal = "17";
  TInt<PropertiesTest> b(val + 1);

  b = sVal;

  ASSERT_EQ(val, b);
}

TEST(Properties_TInt, is_Int_operator_ToString_working)
{
  int val = 17;
  string result = "17";
  TInt<PropertiesTest> b(val);

  ASSERT_EQ(result, (string)b);
}

/**
Проверка класса TInt
*/

TEST(Properties_TDouble, can_create_default_Double)
{
  ASSERT_NO_THROW(TDouble<PropertiesTest> a);
}

TEST(Properties_TDouble, can_create_Double)
{
  double val = 17.3;
  ASSERT_NO_THROW(TDouble<PropertiesTest> a(val));
}

TEST(Properties_TDouble, is_init_Double_value)
{
  double val = 17.3;
  TDouble<PropertiesTest> b(val);
  ASSERT_EQ(val, (double)b);
}

TEST(Properties_TDouble, is_Double_operator_FromString_working)
{
  double val = 17.3;
  string sVal = "17.3";
  TDouble<PropertiesTest> b(val + 1);

  b = sVal;

  ASSERT_EQ(val, (double)b);
}

TEST(Properties_TDouble, is_Double_operator_ToString_working)
{
  double val = 17.987654;
  string result = "17.987654";
  TDouble<PropertiesTest> b(val);

  ASSERT_EQ(result, (string)b);
}

/**
Проверка класса TString
*/

TEST(Properties_TString, can_create_default_String)
{
  ASSERT_NO_THROW(TString<PropertiesTest> a);
}

TEST(Properties_TString, can_create_String)
{
  string val = "abc";
  ASSERT_NO_THROW(TString<PropertiesTest> a(val));
}

TEST(Properties_TString, is_init_String_value)
{
  string val = "abc";
  TString<PropertiesTest> b(val);
  ASSERT_EQ(val, b.GetData());
}

/**
Проверка класса TStrings
*/

TEST(Properties_TStrings, can_create_default_Strings)
{
  ASSERT_NO_THROW(TStrings<PropertiesTest> a);
}

TEST(Properties_TStrings, is_init_Strings_value)
{
  string val[3] = { "abc", "def", "gih" };
  TStrings<PropertiesTest> b;
  b.SetSize(3);
  b = val;
  for (int i = 0; i < 3; i++)
    ASSERT_EQ(val[i], b.GetData()[i]);
}


TEST(Properties_TStrings, is_Strings_operator_FromString_working)
{
  string val = "a_b_c";
  string result[] = {"a", "b", "c"};
  TStrings<PropertiesTest> b;

  b = val;
  for (int i = 0; i < 3; i++)
    ASSERT_EQ(result[i], b.GetData()[i]);
  //ASSERT_EQ(val, b);
}

TEST(Properties_TStrings, is_Strings_operator_ToString_working)
{
  string sVal[] = { "a", "b", "c" };
  string result = "a_b_c";
  TStrings<PropertiesTest> b(sVal, 3);

  ASSERT_EQ(result, (string)b);
}

/**
Проверка класса TInts
*/

TEST(Properties_TInts, can_create_Ints)
{
  ASSERT_NO_THROW(TInts<PropertiesTest> a);
}

TEST(Properties_TInts, is_init_Ints_value)
{
  int val[3] = { 1, 2, 3 };
  TInts<PropertiesTest> b;
  b.SetSize(3);
  b = val;
  for (int i = 0; i < 3; i++)
    ASSERT_EQ(val[i], b.GetData()[i]);
}

TEST(Properties_TInts, is_Ints_operator_FromString_working)
{
  string val = "1_2_3";
  int result[] = { 1, 2, 3 };
  TInts<PropertiesTest> b;

  b = val;
  for (int i = 0; i < 3; i++)
    ASSERT_EQ(result[i], b.GetData()[i]);
}

TEST(Properties_TInts, is_Ints_operator_ToString_working)
{
  int sVal[] = { 1, 2, 3 };
  string result = "1_2_3";
  TInts<PropertiesTest> b(sVal, 3);

  ASSERT_EQ(result, (string)b);
}

/**
Проверка класса TDoubles
*/

TEST(Properties_TDoubles, can_create_Doubles)
{
  ASSERT_NO_THROW(TDoubles<PropertiesTest> a);
}

TEST(Properties_TDoubles, is_init_Doubles_value)
{
  double val[3] = { 1.1, 2.30, 3.54 };
  ASSERT_NO_THROW(TDoubles<PropertiesTest> a);
  TDoubles<PropertiesTest> b;
  b.SetSize(3);
  b = val;
  for (int i = 0; i < 3; i++)
    ASSERT_EQ(val[i], b.GetData()[i]);
}


TEST(Properties_TDoubles, is_TDoubles_operator_FromString_working)
{
  string val = "1.1_2.2_3.3";
  double result[] = { 1.1, 2.2, 3.3 };
  TDoubles<PropertiesTest> b;

  b = val;
  for (int i = 0; i < 3; i++)
    ASSERT_EQ(result[i], b.GetData()[i]);
}

TEST(Properties_TDoubles, is_TDoubles_operator_ToString_working)
{
  double sVal[] = { 1.1, 2.2, 3.3 };
  string result = "1.100000_2.200000_3.300000";
  TDoubles<PropertiesTest> b(sVal, 3);

  ASSERT_EQ(result, (string)b);
}

/**
Проверка класса TETypeMethod
*/

TEST(Properties_TETypeMethod, can_create_default_ETypeMethod)
{
  ASSERT_NO_THROW(TETypeMethod<PropertiesTest> a);
}

TEST(Properties_TETypeMethod, can_create_ETypeMethod)
{
  ETypeMethod val = ManyNumPointMethod;
  ASSERT_NO_THROW(TETypeMethod<PropertiesTest> a(val));
}

TEST(Properties_TETypeMethod, is_init_ETypeMethod_value)
{
  ETypeMethod val = ManyNumPointMethod;
  TETypeMethod<PropertiesTest> b(val);
  ASSERT_EQ(val, b);
}

/**
Проверка класса TESeparableMethodType
*/

TEST(Properties_TESeparableMethodType, can_create_default_ESeparableMethodType)
{
  ASSERT_NO_THROW(TESeparableMethodType<PropertiesTest> a);
}

TEST(Properties_TESeparableMethodType, can_create_ESeparableMethodType)
{
  ESeparableMethodType val = GridSearch;
  ASSERT_NO_THROW(TESeparableMethodType<PropertiesTest> a(val));
}

TEST(Properties_TESeparableMethodType, is_init_ESeparableMethodType_value)
{
  ESeparableMethodType val = GridSearch;
  TESeparableMethodType<PropertiesTest> b(val);
  ASSERT_EQ(val, b);
}

/**
Проверка класса TELocalMethodScheme
*/

TEST(Properties_TELocalMethodScheme, can_create_default_ELocalMethodScheme)
{
  ASSERT_NO_THROW(TELocalMethodScheme<PropertiesTest> a);
}


/**
Проверка класса TEStopCondition
*/

TEST(Properties_TEStopCondition, can_create_default_EStopCondition)
{
  ASSERT_NO_THROW(TEStopCondition<PropertiesTest> a);
}

TEST(Properties_TEStopCondition, can_create_EStopCondition)
{
  EStopCondition val = OptimumVicinity2;
  ASSERT_NO_THROW(TEStopCondition<PropertiesTest> a(val));
}

TEST(Properties_TEStopCondition, is_init_EStopCondition_value)
{
  EStopCondition val = OptimumVicinity2;
  TEStopCondition<PropertiesTest> b(val);
  ASSERT_EQ(val, b);
}

/**
Проверка класса TETypeCalculation
*/

TEST(Properties_TETypeCalculation, can_create_default_ETypeCalculation)
{
  ASSERT_NO_THROW(TETypeCalculation<PropertiesTest> a);
}

TEST(Properties_TETypeCalculation, can_create_ETypeCalculation)
{
  ETypeCalculation val = CUDA;
  ASSERT_NO_THROW(TETypeCalculation<PropertiesTest> a(val));
}

TEST(Properties_TETypeCalculation, is_init_ETypeCalculation_value)
{
  ETypeCalculation val = CUDA;
  TETypeCalculation<PropertiesTest> b(val);
  ASSERT_EQ(val, b);
}

/**
Проверка класса TETypeProcess
*/

TEST(Properties_TETypeProcess, can_create_default_ETypeProcess)
{
  ASSERT_NO_THROW(TETypeProcess<PropertiesTest> a);
}

TEST(Properties_TETypeProcess, can_create_ETypeProcess)
{
  ETypeProcess val = SynchronousProcessNew;
  ASSERT_NO_THROW(TETypeProcess<PropertiesTest> a(val));
}

TEST(Properties_TETypeProcess, is_init_ETypeProcess_value)
{
  ETypeProcess val = SynchronousProcessNew;
  TETypeProcess<PropertiesTest> b(val);
  ASSERT_EQ(val, b);
}

/**
Проверка класса TEMapType
*/

TEST(Properties_TEMapType, can_create_default_EMapType)
{
  ASSERT_NO_THROW(TEMapType<PropertiesTest> a);
}

TEST(Properties_TEMapType, can_create_EMapType)
{
  EMapType val = mpRotated;
  ASSERT_NO_THROW(TEMapType<PropertiesTest> a(val));
}

TEST(Properties_TEMapType, is_init_EMapType_value)
{
  EMapType val = mpRotated;
  TEMapType<PropertiesTest> b(val);
  ASSERT_EQ(val, b);
}