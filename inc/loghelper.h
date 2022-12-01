#ifndef LOGHELPER_H_
#define LOGHELPER_H_

#include <sstream>
#include <boost/format.hpp>
#include <iostream>
using namespace std;


enum log_level_t {
    LOG_CRITICAL = 0,
    LOG_ERROR = 1,
    LOG_WARNING = 2,
    LOG_INFO = 3,
    LOG_DEBUG = 4
};

extern log_level_t GLOBAL_LEVEL;

namespace log_impl {


    class formatted_log_t {
        public:
            formatted_log_t( log_level_t level, const wchar_t* msg ) : fmt(msg), level(level) {}
            ~formatted_log_t() {
                // GLOBAL_LEVEL is a global variable and could be changed at runtime
                // Any customization could be here
                if ( level <= GLOBAL_LEVEL ) wcout << level << L" " << fmt << endl;
            }        
            template <typename T> 
                formatted_log_t& operator %(T value) {
                    fmt % value;
                    return *this;
                }    
        protected:
            log_level_t     level;
            boost::wformat      fmt;
    };
}//namespace log_impl
// Helper function. Class formatted_log_t will not be used directly.
template <log_level_t level>
log_impl::formatted_log_t log(const wchar_t* msg) {
    return log_impl::formatted_log_t( level, msg );
}

#endif
