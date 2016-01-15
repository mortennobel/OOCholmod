#include "config_singleton.h" // used for correct syntax highlight

#include <memory>

namespace oocholmod {

    using namespace std;


    string& lastError(){
        static string lastError;
        return lastError;
    }

    unique_ptr<cholmod_common>& common(){
        static unique_ptr<cholmod_common> common;
        return common;
    }



    inline void errorHandler(int status, const char *file,
            int line, const char *message){
        lastError() = string{message};
    }

    inline cholmod_common *ConfigSingleton::getCommonPtr(){
        if (!common().get()){
            common().reset(new cholmod_common());
            cholmod_start(common().get());
            common()->error_handler = errorHandler;
        }
        return common().get();
    }

    inline void ConfigSingleton::config(cholmod_common *config){
        destroy();
        common().reset(config);
        cholmod_start(common().get());
    }

    inline void ConfigSingleton::destroy(){
        if (common().get()){
            cholmod_finish(common().get()) ;
            common().reset(nullptr);
        }
    }

    inline string ConfigSingleton::getLastError(){
        return lastError();
    }
}

