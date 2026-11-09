#ifndef EPM_MODELS_ORIGINALIMPL_STATEOBSERVER_H_INCLUDED
#define EPM_MODELS_ORIGINALIMPL_STATEOBSERVER_H_INCLUDED

#include <string>
#include <vector>

#include "Util/ForwardDecl.hpp"

namespace EPM::Models::OriginalImpl
{
    /* History of the EPM state on the outer time grid. The caller records
     * each state explicitly; this class is not an odeint observer. odeint
     * copies observers by value, which made the old counter-based observer
     * record the wrong states (TASK-3). */
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

    /* odeint observer that sets a flag once the plume reaches water
     * saturation (RH_w >= 1). odeint calls it at every adaptive step, so it
     * also catches a saturation peak between two outer grid points. It holds
     * only a reference to an external flag and no state of its own: odeint
     * copies observers by value, and every copy writes to the same flag. */
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
