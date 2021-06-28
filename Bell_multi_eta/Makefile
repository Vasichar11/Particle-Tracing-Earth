BIN=tracer   
CC=g++ 
LNK=g++ -fopenmp  -g
OPT=-O3
OMP=-fopenmp
CXXFLAGS=-std=c++17 -lstdc++fs -Wall


#Include paths and flags for libraries used
INC_HIVIFE= -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I/home/family/HighFive/include/
#Add paths that correspond to the running system (e.g -I/home/"username"/HighFive/include/ ).


#Objects to build 
OBJ= Bell_multi_eta.o   functions.o   time_rates.o    struct_Particles.o   struct_Species.o


#Rules

defaul:tracer

clean: 
	rm  *.o *.gch

allclean:
	rm *.o ${BIN} *h5 *.gch

struct_Particles.o: struct_Particles.cc
	${CC} ${CXXFLAGS} ${OPT} -c struct_Particles.cc   

struct_Species.o: struct_Species.cc
	${CC} ${CXXFLAGS} ${OPT} -c struct_Species.cc 

functions.o: functions.cc
	${CC} ${CXXFLAGS} ${OPT} -c functions.cc

time_rates.o: time_rates.cc
	${CC} ${CXXFLAGS} ${OPT} -c time_rates.cc

Bell_multi_eta.o: Bell_multi_eta.cc
	${CC} ${CXXFLAGS} ${OPT} -c Bell_multi_eta.cc    ${INC_HIVIFE}

tracer: $(OBJ)
	${LNK} -o ${BIN} ${OBJ} ${INC_HIVIFE} ${INC_READER} ${CXXFLAGS} ${OMP} ${INC_INI}