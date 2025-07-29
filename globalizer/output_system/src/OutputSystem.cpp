#include "OutputSystem.h"

#include <stdexcept>

// ------------------------------------------------------------------------------------------------
OutputDebugHelper print_dbg;
OutputHelper print;
LogDebugHelper print_dbg_file;
OutputLevel1Helper print_l1;
OutputLevel2Helper print_l2;

// ------------------------------------------------------------------------------------------------
bool Logger::mIsInitialized = false;
std::string Logger::mLogFileName = std::string();

// ------------------------------------------------------------------------------------------------
void Logger::init(const std::string& logFileName)
{
#ifdef _DEBUG
  if(logFileName.empty())
    throw std::runtime_error("Empty log file name in Logger");
  mIsInitialized = true;
  mLogFileName = logFileName;
  Logger& instance = Logger::instance();
  instance.flush();
#endif
}

// ------------------------------------------------------------------------------------------------
void Logger::flush()
{
  mLogOutputStream.flush();
}

// ------------------------------------------------------------------------------------------------
Logger& Logger::instance()
{
  static Logger instance;
  return instance;
}

// ------------------------------------------------------------------------------------------------
Logger::Logger()
{
#ifdef _DEBUG
  if(mIsInitialized)
  {
    if (!mLogOutputStream.is_open())
    {
      mLogOutputStream.open(mLogFileName.c_str(), std::fstream::out);
      if (!mLogOutputStream.is_open())
        throw std::runtime_error("Unable to create a debug log file\n");
    }
  }
  else
  {
    throw std::runtime_error("Logger singleton instance is not initialized\n");
  }
#endif
}