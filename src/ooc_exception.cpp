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
    
    OOCException *OOCException::lastException = nullptr;
    
    void OOCException::createOOCException(std::string what){
        clearLastException();
        lastException = new OOCException(what);
        if (ConfigSingleton::isUsingException()){
            throw OOCException{what};
        }
    }
    
    OOCException::OOCException(std::string what)
    :whatMsg{what}
    {}
    
    
    const char* OOCException::what() const throw() {
        return whatMsg.c_str();
    }
    
    OOCException* OOCException::getLastException(){
        return lastException;
    }
    
    void OOCException::clearLastException(){
        delete lastException;
        lastException = nullptr;
    }
}