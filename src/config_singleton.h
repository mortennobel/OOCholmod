//
//  config_singleton.h
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//

#pragma once

#include <iostream>
#include <cholmod.h>

namespace oocholmod {
    
    class ConfigSingleton{
    public:
        static cholmod_common *getCommonPtr();
        static void destroy();
    private:
        static std::unique_ptr<cholmod_common> common;
    };
    
}