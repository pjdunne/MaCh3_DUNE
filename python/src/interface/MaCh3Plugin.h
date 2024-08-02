#pragma once
/*
Derived classes are interfaced with via plugins
*/
#include <memory>
#include <iostream>

#include "manager/manager.h"

template <typename T>
class MaCh3Plugin{
   public:
    MaCh3Plugin(){
        override_default_constructor = false;
    }

    virtual void process(std::unique_ptr<T> obj)=0;

    bool get_override_default_constructor(){return override_default_constructor;}

    template<class... Args>
    std::unique_ptr<T> construct_obj(const Args&... args){
        std::cerr<<"Not Implemented Error: Constructor not implemented here"<<std::endl;
        std::cerr<<__FILE__<<":"<<__LINE__<<std::endl;
        throw;
    }

   protected:
    bool override_default_constructor;

};