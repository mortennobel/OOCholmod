//
//  config_singleton.cpp
//  OOCholmod
//
//  Created by Morten Nobel-Jørgensen / Asger Nyman Christiansen
//  Copyright (c) 2013 Morten Nobel-Joergensen. All rights reserved.
//

#include "config_singleton.h"

#include <memory>

namespace oocholmod {
    
    using namespace std;
    
    unique_ptr<cholmod_common> common;
    string lastError;
    
    bool ConfigSingleton::exceptionEnabled = true;
    
    void errorHandler(int status, const char *file,
                      int line, const char *message){
        lastError = string{message};
    }
    
    bool ConfigSingleton::isUsingException(){
        return exceptionEnabled;
    }
    
    void ConfigSingleton::setUseException(bool use){
        exceptionEnabled = use;
    }

    
    cholmod_common *ConfigSingleton::getCommonPtr(){
        if (!common.get()){
            common.reset(new cholmod_common());
            cholmod_start(common.get());
            common->error_handler = errorHandler;
            
        }
        return common.get();
    }
    
    void ConfigSingleton::config(cholmod_common *config){
        destroy();
        common.reset(config);
        cholmod_start(common.get());
    }
    
    void ConfigSingleton::destroy(){
        if (common.get()){
            cholmod_finish(common.get()) ;
            common.reset(nullptr);
        }
    }
    
    string ConfigSingleton::getLastError(){
        return lastError;
    }
}

