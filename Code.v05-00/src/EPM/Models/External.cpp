#include <stdexcept>

#include "EPM/Models/External.hpp"


EPM::Models::External::External(const OptInput &optInput) : Base(optInput) {
    throw std::runtime_error("External EPM model is not implemented yet.");
}
