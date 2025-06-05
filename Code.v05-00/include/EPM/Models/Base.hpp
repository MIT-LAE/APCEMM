#ifndef EPM_MODELS_BASE_HPP_INCLUDED
#define EPM_MODELS_BASE_HPP_INCLUDED

#include "Core/Input_Mod.hpp"


namespace EPM::Models {

class Base {
public:
    Base() = delete;
    Base(const OptInput &optInput);
    virtual ~Base() = default;

    // Disable copy and move semantics: EPM model instances will only ever be
    // accessed via std::unique_ptr values.
    Base(const Base&) = delete;
    Base& operator=(const Base&) = delete;
    Base(Base&&) = delete;
    Base& operator=(Base&&) = delete;

    // Virtual function to be overridden by derived classes
    virtual void run() = 0;
};

}

#endif
