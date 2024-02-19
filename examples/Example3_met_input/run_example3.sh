cd ../../
mkdir build
cd build
cmake ../Code.v05-00 && cmake --build . || exit 1
cd ../examples/Example3_met_input
export APCEMM_runDir="."
./../../build/APCEMM input.yaml
