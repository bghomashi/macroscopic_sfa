#ifndef __LOGGER_H__
#define __LOGGER_H__

#include <string>

#define LOG_LEVEL_NONE        0
#define LOG_LEVEL_INFO        1
#define LOG_LEVEL_WARN        2
#define LOG_LEVEL_CRITICAL    3
#define LOG_LEVEL_DEBUG       4

#ifndef LOG_LEVEL
#define LOG_LEVEL LOG_LEVEL_DEBUG
#endif


namespace Log {
    void info(const std::string& text);
    void warn(const std::string& text);
    void critical(const std::string& text);
    void debug(const std::string& text);
    void set_logger(const std::string& log_file);
}

#if LOG_LEVEL >= LOG_LEVEL_INFO
#define LOG_INFO(x) Log::info(x)
#else
#define LOG_INFO(x)
#endif

#if LOG_LEVEL >= LOG_LEVEL_WARN
#define LOG_WARN(x) Log::warn(x)
#else
#define LOG_WARN(x)
#endif

#if LOG_LEVEL >= LOG_LEVEL_CRITICAL
#define LOG_CRITICAL(x) Log::critical(x)
#else
#define LOG_CRITICAL(x)
#endif

#if LOG_LEVEL >= LOG_LEVEL_DEBUG
#define LOG_DEBUG(x) Log::debug(x)
#else
#define LOG_DEBUG(x)
#endif

#endif