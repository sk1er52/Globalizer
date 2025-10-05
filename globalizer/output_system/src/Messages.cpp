#include "Messages.h"

#include <iostream>
#include <cstdlib>
#include <string>
#include <sstream>
#include <time.h>
#include <stdexcept>

bool OutputMessage::mIsInit = false;
bool OutputMessage::mIsMPIInit = false;
int OutputMessage::mProcessCount = 1;
int OutputMessage::mProcessNumber = 0;
std::vector<bool> OutputMessage::mIsPrintToFile = std::vector<bool>();
bool OutputMessage::mOldIsPrintToFile = false;
std::string* OutputMessage::mErrorsName = 0;
int* OutputMessage::mErrorsCode = 0;
int OutputMessage::mErrorsCount = 0;
std::vector<std::string> OutputMessage::mLogFileName = std::vector<std::string>();
int OutputMessage::rootNumber = 0;
std::vector<std::fstream*> OutputMessage::mLogOutputStream = std::vector<std::fstream*>();
int OutputMessage::streamCount = 0;
int OutputMessage::ImportantMessageNumber = -1;
bool OutputMessage::mIsPrintToConsole = true;

/// Глобальная переменная для вывода информации на экран, в debug печатает в лог-файл
OutputMessage printMessage;
/**
Глобальная переменная для вывода информации на экран только из первого процесса,
в debug печатает в лог-файл
*/
OutputMessage printMessageInFirstProcess(1);
/**
Глобальная переменная для вывода информации на экран в корне,
в debug печатает в лог-файл
*/
OutputMessage printMessageInRoot(true);

// ------------------------------------------------------------------------------------------------
void OutputMessage::SetDefaultErrors()
{
  mErrorsCount = 4;
  if (mErrorsName != 0)
    delete[] mErrorsName;
  if (mErrorsCode != 0)
    delete[] mErrorsCode;
  mErrorsName = new std::string[mErrorsCount];
  mErrorsName[0] = "";
  mErrorsName[1] = "SYSTEM CRASH";
  mErrorsName[2] = "INCORRECT ARGUMENT SIZE";
  mErrorsName[3] = "INCORRECT DIM IN TASK LEVEL";
  mErrorsCode = new int[mErrorsCount];
  mErrorsCode[0] = SYSTEM_OK;
  mErrorsCode[1] = SYSTEM_CRASH;
  mErrorsCode[2] = ERROR_ARGUMENT_SIZE;
  mErrorsCode[3] = ERROR_DIM_IN_TASK_LEVEL;
}

// ------------------------------------------------------------------------------------------------
OutputMessage::OutputMessage(int printPocessNumber)
{
  mPrintPocessNumber = printPocessNumber;
  output = "";
}

// ------------------------------------------------------------------------------------------------
OutputMessage::OutputMessage(bool printInRoot)
{
  mPrintPocessNumber = rootNumber;
  output = "";
}

// ------------------------------------------------------------------------------------------------
void OutputMessage::Init(bool isMPIInit, std::string logFileName, int processCount,
  int processNumber, bool isPrintToFile, std::string* errorsName,
  int* errorsCode, int errorsCount)
{
  streamCount = 2;
  mIsPrintToFile.resize(streamCount);
#ifdef _DEBUG
  mIsPrintToFile[0] = isPrintToFile;
  mIsPrintToFile[1] = false;
#else
  mIsPrintToFile[0] = false;
  mIsPrintToFile[1] = isPrintToFile;
#endif
  std::string log = logFileName + "_process_" + toString(processNumber);

  mLogFileName.resize(streamCount);
  mLogFileName[0] = log + "_Debug.txt";
  mLogFileName[1] = log + ".txt";
  Logger::init(mLogFileName[0]);

  //mLogOutputStream.resize(streamCount);

  for (int k = 0; k < streamCount; k++)
  {
    std::fstream* fs = new std::fstream(mLogFileName[streamCount - 1].c_str(), std::fstream::out);
    mLogOutputStream.push_back(fs);
  }
  if (mIsPrintToFile[1])
  {
    mLogOutputStream[1]->open(mLogFileName[1].c_str(), std::fstream::out);
    if (!mLogOutputStream[1]->is_open())
      throw std::runtime_error("Unable to create a log file\n");
    mLogOutputStream[1]->flush();
  }

  mProcessCount = processCount;
  mProcessNumber = processNumber;
  mErrorsCount = errorsCount;

  if ((mErrorsName == 0) || (mErrorsCount == 0) || (mErrorsCode == 0))
    SetDefaultErrors();
  else
  {
    if (mErrorsName != 0)
      delete[] mErrorsName;
    if (mErrorsCode != 0)
      delete[] mErrorsCode;

    mErrorsName = new std::string[mErrorsCount];
    mErrorsCode = new int[mErrorsCount];

    for (int i = 0; i < mErrorsCount; i++)
    {
      mErrorsName[i] = errorsName[i];
      mErrorsCode[i] = errorsCode[i];
    }
  }

  mIsMPIInit = isMPIInit;
  rootNumber = 0;
  mIsInit = true;

  printMessage << "";
  printMessageInFirstProcess << "";
  printMessageInRoot << "";
}

// ------------------------------------------------------------------------------------------------
std::string OutputMessage::ErrorName(int error)
{
  for (int i = 0; i < mErrorsCount; i++)
  {
    if (mErrorsCode[i] == error)
      return mErrorsName[i];
  }
  return toString(error);
}

// ------------------------------------------------------------------------------------------------
void OutputMessage::PrintError(int error, std::string file, int line, std::string expression)
{
  if (error != SYSTEM_OK)
  {
    if (!mIsInit)
    {
      int isMPIInit = 0;
      MPI_Initialized(&isMPIInit);
      int mpiC = 1, mpiR = 0;
      if (isMPIInit)
      {
        MPI_Comm_size(MPI_COMM_WORLD, &mpiC);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiR);
      }
      time_t t = time(0);
      Init(isMPIInit == 0 ? false : true, "Erorlog" + toString(t) + ".txt", mpiC, mpiR, true);
      ImportantMessageAppears();
    }
    else
    {
      time_t t = time(0);
      ImportantMessageAppears("Erorlog_process_" + toString(mProcessNumber) + "_" + toString(t) + ".txt");
    }
    printMessage << "Expression: " << expression << "\n" << file << "\t" << line << "\n";
    printMessage << "\n!!! Error " << ErrorName(error) << " !!!\n";
    ImportantMessageEnded();
    exit(error);
  }
}

// ------------------------------------------------------------------------------------------------
void OutputMessage::PrintImportantMessage(std::string expression,
  std::string file, int line, std::string fileName)
{
  if (!mIsInit)
  {
    int isMPIInit = 0;
    MPI_Initialized(&isMPIInit);
    int mpiC = 1, mpiR = 0;
    if (isMPIInit)
    {
      MPI_Comm_size(MPI_COMM_WORLD, &mpiC);
      MPI_Comm_rank(MPI_COMM_WORLD, &mpiR);
    }
    time_t t = time(0);
    Init(isMPIInit == 0 ? false : true, "Erorlog" + toString(t) + ".txt", mpiC, mpiR, true);
  }
  ImportantMessageAppears(fileName);
  if ((file != "") || (line != -1))
    printMessage << "Message: " << expression << "\n" << file << "\t" << line << "\n";
  else
    printMessage << "Message: " << expression << "\n";
  ImportantMessageEnded();
}

// ------------------------------------------------------------------------------------------------
void OutputMessage::PrintMessageToFile(std::string expression, std::string file,
  int line, std::string fileName)
{
  mIsPrintToConsole = false;
  ImportantMessageAppears(fileName);
  if ((file != "") || (line != -1))
    printMessage << "Message: " << expression << "\n" << file << "\t" << line << "\n";
  else
    printMessage << "Message: " << expression << "\n";

  mIsPrintToConsole = true;
  ImportantMessageEnded();
}

// ------------------------------------------------------------------------------------------------
void OutputMessage::ImportantMessageAppears(std::string importantMessageFile)
{
  if (importantMessageFile == "")
  {
#ifdef _DEBUG
    importantMessageFile = mLogFileName[0];
#else
    importantMessageFile = mLogFileName[1];
#endif
  }
  bool isImportantMessageFileOpen = false;
  for (int i = 0; i < streamCount; i++)
  {
    if (mLogFileName[i] == importantMessageFile)
    {
      isImportantMessageFileOpen = true;
      mOldIsPrintToFile = mIsPrintToFile[i];
      ImportantMessageNumber = i;
      mIsPrintToFile[i] = true;
      break;
    }
  }
  if (!isImportantMessageFileOpen)
  {
    streamCount++;
    mLogFileName.resize(streamCount);
    mLogFileName[streamCount - 1] = importantMessageFile;

    //mLogOutputStream.resize(streamCount);
    std::fstream* fs = new std::fstream(mLogFileName[streamCount - 1].c_str(), std::fstream::out);
    mLogOutputStream.push_back(fs);
    
    mLogOutputStream[mLogOutputStream.size() - 1]->open(mLogFileName[streamCount - 1].c_str(),
      std::fstream::out);
    if (!mLogOutputStream[mLogOutputStream.size() - 1]->is_open())
      throw std::runtime_error("Unable to create a log file\n");
    mLogOutputStream[mLogOutputStream.size() - 1]->flush();
    mIsPrintToFile.resize(streamCount);
    mIsPrintToFile[streamCount - 1] = true;
    ImportantMessageNumber = streamCount - 1;
    mOldIsPrintToFile = false;
  }
}

// ------------------------------------------------------------------------------------------------
void OutputMessage::OpenFileStream(std::string fileName)
{
  if (fileName == "")
  {
#ifdef _DEBUG
    fileName = mLogFileName[0];
#else
    fileName = mLogFileName[1];
#endif
  }

  bool isFileOpen = false;
  for (int i = 0; i < streamCount; i++)
  {
    if (mLogFileName[i] == fileName)
    {
      isFileOpen = true;
      mIsPrintToFile[i] = true;
      break;
    }
  }
  if (!isFileOpen)
  {
    streamCount++;
    mLogFileName.resize(streamCount);
    mLogFileName[streamCount - 1] = fileName;

    //mLogOutputStream.resize(streamCount);
    std::fstream* fs = new std::fstream(mLogFileName[streamCount - 1].c_str(), std::fstream::out);
    mLogOutputStream.push_back(fs);
    mLogOutputStream[mLogOutputStream.size() - 1]->open(mLogFileName[streamCount - 1].c_str(),
      std::fstream::out);
    if (!mLogOutputStream[mLogOutputStream.size() - 1]->is_open())
      throw std::runtime_error("Unable to create a log file\n");
    mLogOutputStream[mLogOutputStream.size() - 1]->flush();
    mIsPrintToFile.resize(streamCount);
    mIsPrintToFile[streamCount - 1] = true;
  }
}

// ------------------------------------------------------------------------------------------------
void OutputMessage::BreakFileStream(std::string fileName)
{
  if (fileName == "")
  {
#ifdef _DEBUG
    fileName = mLogFileName[0];
#else
    fileName = mLogFileName[1];
#endif
  }

  for (int i = 0; i < streamCount; i++)
  {
    if (mLogFileName[i] == fileName)
    {
      mIsPrintToFile[i] = false;
      break;
    }
  }
}

// ------------------------------------------------------------------------------------------------
void OutputMessage::ImportantMessageEnded()
{
  mIsPrintToFile[ImportantMessageNumber] = mOldIsPrintToFile;
  mOldIsPrintToFile = false;
}
