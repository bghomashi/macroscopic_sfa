#include <iostream>
#include <fstream>
#include "logger.h"

static std::ofstream s_logger;

namespace Log {
    void info(const std::string& text) {
        if (!s_logger.bad())
            s_logger << text << "\n";
        else
            std::cout << "INFO: " << text << "\n";
    }
    void warn(const std::string& text) {
        if (!s_logger.bad())
            s_logger << text << "\n";
        else
            std::cout << "WARN: " << text << "\n";
    }
    void critical(const std::string& text) {
        if (!s_logger.bad())
            s_logger << text << "\n";
        else
            std::cout << "CRITICAL: " << text << "\n";
    }
    void debug(const std::string& text) {
        if (!s_logger.bad())
            s_logger << text << "\n";
        else
            std::cout << "DEBUG: " << text << "\n";
    }
    void set_logger_file(const std::string& log_file) {
        s_logger = std::ofstream(log_file);
    }
    void flush() {
        s_logger << std::endl;
    }
}
