# # Uncomment the following lines to compile APCEMM
# git submodule update --init --recursive
# mkdir ../../build
# cd ../../build/
# cmake ../Code.v05-00 && cmake --build . || exit 1
# cd ../examples/Example2_Impose_Depth

# Run Example 2
export APCEMM_runDir="."
./../../build/APCEMM input.yaml
