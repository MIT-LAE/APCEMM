#ifndef EPM_MODELS_NEW_PHYSICS_HPP_INCLUDED
#define EPM_MODELS_NEW_PHYSICS_HPP_INCLUDED

#include <stdexcept>

#include "EPM/Models/Base.hpp"


namespace EPM::Models {

class NewPhysics : public Base {
public:
    NewPhysics(const OptInput &optInput) : Base(optInput) {
        throw std::runtime_error("New EPM model is not implemented yet.");
    }

    void run() override {}
};

} // namespace EPM::Models

#endif
