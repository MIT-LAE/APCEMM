#ifndef EPM_MODELS_ORIGINALIMPL_STATEOBSERVER_H_INCLUDED
#define EPM_MODELS_ORIGINALIMPL_STATEOBSERVER_H_INCLUDED

#include <string>
#include <vector>

#include "Util/ForwardDecl.hpp"

namespace EPM::Models::OriginalImpl
{
    /* History of the EPM state on the outer time grid */
    class StateHistory
    {
    public:
        StateHistory(std::vector<UInt> indices, std::string fileName);

        void record(const Vector_1D &x, double t);
        const Vector_2D& states() const { return m_states; }
        const Vector_1D& times() const { return m_times; }

        void print2File() const;

    private:
        Vector_2D m_states;
        Vector_1D m_times;
        std::string m_fileName;
        const std::vector<UInt> m_indices;
    };

    /* Boost odeint observer that sets a flag once the plume reaches water
     * saturation. Called every adaptive step by the integrator
     * so that it can catch a saturation peak between two outer grid points.
     * Contains only a reference to an external boolean flag because the boost
     * integrator copies observers by value. */
    class WaterSaturationObserver
    {
    public:
        explicit WaterSaturationObserver(bool &reached): m_reached(reached) {}
        void operator()(const Vector_1D &x, double t) const;

    private:
        bool &m_reached;
    };
}

#endif
