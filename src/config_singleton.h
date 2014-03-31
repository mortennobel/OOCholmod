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
#include <string>

namespace oocholmod {
    class DenseMatrix;
    class Factor;
    
	class ConfigSingleton{
    public:
		static void config(cholmod_common *);
		static cholmod_common *getCommonPtr();
		static void destroy();
    private:
        static std::string getLastError();
        ConfigSingleton(){};
        friend class DenseMatrix;
        friend class Factor;
    };
    
}