//
//  cholmod_exception.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 31/03/14.
//  Copyright (c) 2014 Morten Nobel-Joergensen. All rights reserved.
//

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
    void OOCException::createOOCException(std::string what){
        throw OOCException{what};
    }
#endif
    
    OOCException::OOCException(std::string what)
    :whatMsg{what}
    {}
    
    
    const char* OOCException::what() const throw() {
        return whatMsg.c_str();
    }
}