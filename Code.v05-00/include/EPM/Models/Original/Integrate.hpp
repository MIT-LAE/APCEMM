#ifndef EPM_MODELS_ORIGINALIMPL_INTEGRATE_H_INCLUDED
#define EPM_MODELS_ORIGINALIMPL_INTEGRATE_H_INCLUDED

#include <cstddef>

#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>

#include "Util/ForwardDecl.hpp"
#include "EPM/Models/Original/StateObserver.hpp"

namespace EPM::Models::OriginalImpl
{
    /* Advance the EPM state x over one outer interval [t0, t1], then record
     * the end state x(t1) at time t1. The observer sees every adaptive step
     * but does not touch the history. Returns the number of adaptive steps. */
    template <class Stepper, class System, class Observer>
    std::size_t advanceInterval(
        Stepper stepper, System system, Vector_1D &x,
        double t0, double t1, double dt0,
        Observer observer, StateHistory &history)
    {
        const std::size_t nSteps = boost::numeric::odeint::integrate_adaptive(
            stepper, system, x, t0, t1, dt0, observer);
        history.record(x, t1);
        return nSteps;
    }
}

#endif
