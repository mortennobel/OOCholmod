#include "ooc_exception.h"
#include "config_singleton.h"

namespace oocholmod {

#ifdef OOC_NO_EXCEPTION
    OOCException *OOCException::lastException = nullptr;
    
    void OOCException::createOOCException(std::string what){
        clearLastException();
        lastException = new OOCException(what);
    }
    
    OOCException* OOCException::getLastException(){
        return lastException;
    }
    
    void OOCException::clearLastException(){
        delete lastException;
        lastException = nullptr;
    }
#else
    inline void OOCException::createOOCException(std::string what){
        throw OOCException{what};
    }
#endif

    inline OOCException::OOCException(std::string what)
            :whatMsg(what)
    {}


    inline const char* OOCException::what() const throw() {
        return whatMsg.c_str();
    }
}