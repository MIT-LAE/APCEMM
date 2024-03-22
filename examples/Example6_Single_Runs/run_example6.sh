cd ../../Code.v05-00/
cmake . && cmake --build . || exit 1
cd ../examples/Example4_NIPC
export APCEMM_runDir="."
./../../Code.v05-00/APCEMM input.yaml
