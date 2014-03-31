//
//  cholmod_exception.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 31/03/14.
//  Copyright (c) 2014 Morten Nobel-Joergensen. All rights reserved.
//

#include "ooc_exception.h"

namespace oocholmod {
    OOCException::OOCException(std::string what)
    :whatMsg{what}
    {}
    
    
    const char* OOCException::what() const throw() {
        return whatMsg.c_str();
    }
}