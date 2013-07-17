//
//  config_singleton.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//

#include "config_singleton.h"


namespace oocholmod {
    
    std::unique_ptr<cholmod_common> ConfigSingleton::common;
    cholmod_common *ConfigSingleton::getCommonPtr(){
        if (common.get() == nullptr){
            common = std::unique_ptr<cholmod_common>(new cholmod_common());
            cholmod_start(common.get());
        }
        return common.get();
    }
    
    void ConfigSingleton::destroy(){
        if (common.get() != nullptr){
            cholmod_finish(common.get()) ;
            common.reset(nullptr);
        }
    }
}

