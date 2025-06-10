#ifndef EPM_MODELS_ORIGINALIMPL_INDEXES_H_INCLUDED
#define EPM_MODELS_ORIGINALIMPL_INDEXES_H_INCLUDED

namespace EPM::Models::OriginalImpl {

    // Indices for the gas-aerosol system
    enum EPM_ind {
        EPM_ind_Trac = 0, // Tracer
        EPM_ind_T,       // Temperature
        EPM_ind_P,       // Pressure
        EPM_ind_H2O,     // Water vapor mixing ratio
        EPM_ind_SO4,     // Sulfate gas mixing ratio
        EPM_ind_SO4l,    // Sulfate liquid mixing ratio
        EPM_ind_SO4g,    // Sulfate gas mixing ratio
        EPM_ind_SO4s,    // Sulfate solid mixing ratio
        EPM_ind_HNO3,    // Nitric acid mixing ratio
        EPM_ind_Part,    // Particle mixing ratio
        EPM_ind_ParR,    // Particle radius
        EPM_ind_the1,    // Surface coverage 1
        EPM_ind_the2     // Surface coverage 2
    };

} // namespace EPM::Models::OriginalImpl

#endif
