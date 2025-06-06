#ifndef EPM_MODELS_ORIGINALIMPL_STATEOBSERVER_H_INCLUDED
#define EPM_MODELS_ORIGINALIMPL_STATEOBSERVER_H_INCLUDED

#include <string>
#include <vector>

#include "Util/ForwardDecl.hpp"

namespace EPM::Models::OriginalImpl
{
    class StateObserver
    {
    public:
        StateObserver(
            Vector_2D &states, Vector_1D &times,
            std::vector<UInt> indices,
            std::string fileName, UInt write_every = 100);
        StateObserver& operator=(const StateObserver &obs);
        void operator()(const Vector_1D &x, double t);
        double getLastElement() const;
        void print2File() const;
        bool checkwatersat() const;

        UInt m_write_every;
        UInt m_count = 0;
        std::string fileName;

        Vector_2D &m_states;
        Vector_1D &m_times;

    private:
        const std::vector<UInt> m_indices;
    };
}

#endif
