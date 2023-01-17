# Source files that need to be compiled
set(SRCS
    KPP_Function.cpp
    KPP_Hessian.cpp
    KPP_HessianSP.cpp
    KPP_HetRates.cpp
    KPP_Integrator_ADJ.cpp
    KPP_Integrator.cpp
    KPP_Jacobian.cpp
    KPP_JacobianSP.cpp
    KPP_LinearAlgebra.cpp
    KPP_Main_ADJ.cpp
    KPP_Monitor.cpp
    KPP_Photol.cpp
    KPP_Rates.cpp
    KPP_Stoichiom.cpp
    KPP_StoichiomSP.cpp
    KPP_Util.cpp
    )

# This command ensures the static library gets build
add_library(KPP STATIC ${SRCS})

# This command defines the dependencies of libKPP.a
target_link_libraries(KPP Core ${NETCDF_LIBRARIES}
                        ${netCDFCxx_LIBRARIES} OpenMP::OpenMP_CXX)