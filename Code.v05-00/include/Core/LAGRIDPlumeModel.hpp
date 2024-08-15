#ifndef LAGRIDPLUMEMODEL_H
#define LAGRIDPLUMEMODEL_H
#include "AIM/Aerosol.hpp"
#include "AIM/Settling.hpp"
#include "LAGRID/RemappingFunctions.hpp"
#include "FVM_ANDS/FVM_Solver.hpp"
#include "EPM/Integrate.hpp"
#include "Core/Diag_Mod.hpp"
#include "Core/MPMSimVarsWrapper.hpp"
#include "Core/TimestepVarsWrapper.hpp"
#include "Core/Meteorology.hpp"
#include "Util/VectorUtils.hpp"
#include "Util/PlumeModelUtils.hpp"
#include <filesystem>
#include "Core/Status.hpp"
class LAGRIDPlumeModel {
    public:
        static constexpr bool COCIP_MIXING = 0; // Results in less accurate mixing representation, only meant for comparisions vs. the CoCiP model.
        static constexpr double NUM_FILTER_RATIO = 1e-5;
        static constexpr double TOP_BUFFER_SCALING = 1.5;
        static constexpr double BOT_BUFFER_SCALING = 1.1;
        static constexpr double LEFT_BUFFER_SCALING = 1.5;
        static constexpr double RIGHT_BUFFER_SCALING = 1.5;

        LAGRIDPlumeModel() = delete;
        LAGRIDPlumeModel(const OptInput &Input_Opt, const Input &input);
        SimStatus runFullModel();
        SimStatus runEPM();
        struct BufferInfo {
            double leftBuffer;
            double rightBuffer;
            double topBuffer;
            double botBuffer;
        };
    private:
        const OptInput& optInput_;
        const Input& input_;
        int numThreads_;
        SZA sun_;
        Aircraft aircraft_;
        Fuel jetA_;
        Emission EI_;
        MPMSimVarsWrapper simVars_;
        TimestepVarsWrapper timestepVars_;
        AIM::Grid_Aerosol iceAerosol_;
        std::pair<EPM::EPMOutput, SimStatus> EPM_result_;
        Meteorology met_;
        Vector_2D diffCoeffX_;
        Vector_2D diffCoeffY_;
        double survivalFrac_;
        Vector_1D yCoords_;
        Vector_1D yEdges_;
        Vector_1D xCoords_;
        Vector_1D xEdges_;
        Vector_2D H2O_;
        Vector_1D vFall_;
        double initNumParts_;
        double simTime_h_;
        double solarTime_h_;
        double shear_rep_;

        typedef std::pair<std::vector<std::vector<int>>, VectorUtils::MaskInfo> MaskType;
        inline MaskType iceNumberMask(double cutoff_ratio = NUM_FILTER_RATIO) {
            Vector_2D iceTotalNum = iceAerosol_.TotalNumber();
            double maxNum = VectorUtils::VecMax2D(iceTotalNum);
            auto iceNumMaskFunc = [maxNum](double val) {
                return val > maxNum * NUM_FILTER_RATIO;
            };
            return VectorUtils::Vec2DMask(iceTotalNum, xEdges_, yEdges_, iceNumMaskFunc);
        }

        inline MaskType H2OMask() {
            Vector_2D metH2O = met_.H2O_field();
            Vector_2D diffH2O = Vector_2D(yCoords_.size(), Vector_1D(xCoords_.size()));
            for (std::size_t j = 0; j < yCoords_.size(); j++) {
                for (std::size_t i = 0; i < xCoords_.size(); i++) {
                    diffH2O[j][i] = std::abs(H2O_[j][i] - metH2O[j][i]);
                }
            }
            double maxDiff = 1.0e6;
            auto maskFunc = [maxDiff](double val) {
                return val > maxDiff;
            };
            return VectorUtils::Vec2DMask(diffH2O, xEdges_, yEdges_, maskFunc);
        }

        void createOutputDirectories();
        void initializeGrid();
        void saveTSAerosol();
        void initH2O();
        void updateDiffVecs();
        void runTransport(double timestep);
        void remapAllVars(double remapTimestep);
        void trimH2OBoundary();
        std::pair<LAGRID::twoDGridVariable,LAGRID::twoDGridVariable> remapVariable(const VectorUtils::MaskInfo& maskInfo, const BufferInfo& buffers, const Vector_2D& phi, const std::vector<std::vector<int>>& mask);
        double totalAirMass();
        void runCocipH2OMixing(const Vector_2D& h2o_old, const Vector_2D& h2o_amb_new, MaskType& mask_old, MaskType& mask_new);


};

#endif
