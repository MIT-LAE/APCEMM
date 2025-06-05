/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/*                      Early Plume Microphysics                    */
/*                              (EPM)                               */
/*                                                                  */
/* Right-hand side Header File                                      */
/*                                                                  */
/* File                 : RHS.hpp                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#ifndef EPM_MODELS_ORIGINAL_RHS_H_INCLUDED
#define EPM_MODELS_ORIGINAL_RHS_H_INCLUDED

#include "Util/ForwardDecl.hpp"
#include "AIM/Aerosol.hpp"


namespace EPM::Models::Original {

struct gas_aerosol_rhs {
    const double m_temperature_K;
    const double m_pressure_Pa;
    const double m_delta_T;
    const double m_H2O_mixingratio;
    const double m_SO4_mixingratio;
    const double m_SO4l_mixingratio;
    const double m_SO4g_mixingratio;
    const double m_HNO3_mixingratio;
    const double m_part_mixingratio;
    const double m_part_r0;

    const double sticking_SO4;
    const double sigma_SO4;

    const Vector_1D KernelSO4Soot;

    AIM::Aerosol &nPDF_SO4;


    gas_aerosol_rhs(double temperature_K, double pressure_Pa, double delta_T,
                    double H2O_mixingratio, double SO4_mixingratio, double SO4l_mixingratio,
                    double SO4g_mixingratio, double HNO3_mixingratio, double part_mixingratio,
                    double part_r0, Vector_1D Kernel_, AIM::Aerosol &nPDF_SO4_) :
        m_temperature_K(temperature_K),
        m_pressure_Pa(pressure_Pa),
        m_delta_T(delta_T),
        m_H2O_mixingratio (H2O_mixingratio),
        m_SO4_mixingratio (SO4_mixingratio),
        m_SO4l_mixingratio (SO4l_mixingratio),
        m_SO4g_mixingratio (SO4g_mixingratio),
        m_HNO3_mixingratio (HNO3_mixingratio),
        m_part_mixingratio(part_mixingratio),
        m_part_r0(part_r0),
        sticking_SO4(1.0),
        sigma_SO4(5.0E+14),
        KernelSO4Soot(Kernel_),
        nPDF_SO4(nPDF_SO4_) {}

    void operator()(const Vector_1D &x, Vector_1D &dxdt, const double t = 0) const;
};

double entrainmentRate(const double time);
double depositionRate(const double r, const double T, const double P, const double H2O,
                      const double r_0,  const double theta );
double dT_Vortex(const double time, const double delta_T, bool deriv = 0);
bool isFreezable(const double r, const double T, const double H2O,
                 const double r0);
double condensationRate(const double r, const double T, const double P,
                        const double H2O, const double theta);

} // namespace EPM::Models::Original

#endif
