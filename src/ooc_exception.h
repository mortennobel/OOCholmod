//
//  cholmod_exception.h
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen on 31/03/14.
//  Copyright (c) 2014 Morten Nobel-Joergensen. All rights reserved.
//
#pragma once
#include <iostream>
#include <exception>
#include <string>

namespace oocholmod {
    class OOCException : public std::exception {
    public:
        OOCException(std::string what);
        const char* what() const throw() override;
    private:
        std::string whatMsg;
    };
}