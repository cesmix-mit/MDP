!g++ -fPIC -O3 -c cpuGolff.cpp
!g++ --shared cpuGolff.o -o libcpuGolff.so
!ar rvs cpuGolff.a cpuGolff.o        

