echo "-----------------------Deleting build dir-----------------------"
rm -rf build
echo "-----------------------Generating build setup-----------------------"
cmake -B build -S .
echo "-----------------------Compiling and linking-----------------------"
cmake --build build
echo "-----------------------Running-----------------------"
./build/ntt_implementation