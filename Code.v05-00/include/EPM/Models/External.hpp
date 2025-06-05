#ifndef EPM_MODELS_EXTERNAL_HPP_INCLUDED
#define EPM_MODELS_EXTERNAL_HPP_INCLUDED

#include "EPM/Models/Base.hpp"


namespace EPM::Models {

class External : public Base {
public:
    External(const OptInput &optInput);
    void run() override {}
};

}

#endif
