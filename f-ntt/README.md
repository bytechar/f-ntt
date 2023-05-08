# How to build

- ***Make sure you have `cmake` installed. You can install it from brew using `brew install cmake`.***

```bash
cd /path/to/f-ntt
# to build cmake artifacts
cmake -B build -S .
# to compile and link the project
cmake --build build
# to run the executable
./build/ntt_implementation
```
