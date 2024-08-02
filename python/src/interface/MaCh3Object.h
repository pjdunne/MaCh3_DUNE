#pragma once

/*
HW : A Generic MaCh3 object
*/

// C++ Includes
#include <string>
#include <vector>
#include <functional>
#include <iostream>
#include <memory>


#include "MaCh3Plugin.h"
#include "manager.h"

template <typename T>
class MaCh3Object{
   public:
    /// @brief MaCh3Object constructor
    MaCh3Object(std::string manger_name){
        __has_been_processed=false;
        _stored_object = nullptr;

        // Not neat, should be passed as a YAML node
        fit_manager =  new manager(yaml_config);

    }

    /// @brief destructor
    virtual ~MaCh3Object(){}

    /// @brief Pure virtual call around get_likelihood method
    virtual double get_likelihood()=0;

    /// @brief adds processing method to object
    void add_processor(MaCh3Plugin func){
        bool default_constructor_overriden=false;

        if(func.get_override_default_constructor()){
            if(default_constructor_overriden){
                std::cerr<<"ERROR::Tried to overwrite already overwriten constructor, this will not work"<<std::endl;
                std:cerr<<__FILE__<<":"<<__LINE__<<std::endl;
                throw;
            }
            
            // This is stupid but means I can override the default processor
            _processor_vec.insert(_processor_vec.begin(), func);
            default_constructor_overriden = true;
        }
        _processor_vec.push_back(func);
    }

    /// @brief applies processing to object
    void apply_processing(){
        if(__has_been_processed){
            std::cout<<"INFO: Object has been processed, skipping"<<std::endl;
        }

        if(_processor_vec[0]->get_override_default_constructor()){
            _stored_object = _processor_vec[0].construct_obj(fit_manager);
        }
        else{
            construct_default();
        }

        for(auto func : _processor_vec){
            else{

                func->process(*this);
            }
        }
        __has_been_processed = true;
    }


    std::unique_ptr<T> get_stored_object(){
        if(_stored_object==nullptr){
            std::cerr<<"Cannot access object before it has been processed\n";
            std::cerr<<__FILE__<<":"<<__LINE__<<std::endl;
            throw;
        }
        return _stored_object;
    }

    std::unique_ptr<manager> get_manager(){
        return fit_manager;
    }

   protected:

    /// Pure Virtual Config reader method
    virtual void read_config(std::string config)=0;

    /// map of arguments
    std::unordered_map<std::string, std::string> args_map;

    /// Actual MaCh3 object
    std::unique_ptr<T> _stored_object;

    /// Pure Virtual, default constructor of _stored _object, can be overridden with a plugin
    virtual std::unique_ptr<T> construct_default()=0
    
    std::unique_ptr<manager> fit_manager;

   private:
    /// List of processors to be used on _stored object
    std::vector<MaCh3Plugin> _processor_vec;

    bool __has_been_processed;
};