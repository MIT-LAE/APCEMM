#ifndef EPM_MODELS_ORIGINAL_HPP_INCLUDED
#define EPM_MODELS_ORIGINAL_HPP_INCLUDED

#include "EPM/Models/Base.hpp"


namespace EPM::Models {

class Original : public Base {
public:
    Original(const OptInput &optInput);

    void run() override;
};

} // namespace EPM::Models

#endif
