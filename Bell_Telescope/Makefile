BIN=tracer   
CC=g++ 
LNK=g++ -fopenmp  -pg
OPT=-O3
CXXFLAGS=-std=c++17 -lstdc++fs -Wall -pg


INC_HIVIFE= -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I/home/vchariton/HighFive/include/


INC_FAST= -I/usr/include/hdf5/serial -I/home/vchariton/HighFive/include/ -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5


OBJ=  functions.o   time_rates.o   struct_Particles.o   struct_Species.o   struct_Telescope.o


#Rules

default:performant

clean: 
	rm  *.o 

allclean:
	rm *.o ${BIN} *.out *.png

struct_Particles.o: struct_Particles.cc
	${CC} ${CXXFLAGS} ${OPT} -c struct_Particles.cc   

struct_Species.o: struct_Species.cc
	${CC} ${CXXFLAGS} ${OPT} -c struct_Species.cc 

struct_Telescope.o: struct_Telescope.cc
	${CC} ${CXXFLAGS} ${OPT} -c struct_Telescope.cc    ${INC_FAST}

motion.o: motion.cc
	${CC} ${CXXFLAGS} ${OPT} -c motion.cc

functions.o: functions.cc
	${CC} ${CXXFLAGS} ${OPT} -c functions.cc

time_rates.o: time_rates.cc
	${CC} ${CXXFLAGS} ${OPT} -c time_rates.cc

Bell_Telescope_performant.o: Bell_Telescope_performant.cc
	${CC} ${CXXFLAGS} ${OPT} -c Bell_Telescope_performant.cc    ${INC_FAST}

Bell_Telescope_concise.o: Bell_Telescope_concise.cc
	${CC} ${CXXFLAGS} ${OPT} -c Bell_Telescope_concise.cc    ${INC_FAST}

performant: $(OBJ) Bell_Telescope_performant.o
	${LNK} -o ${BIN} ${OBJ}  Bell_Telescope_performant.o ${INC_FAST} ${CXXFLAGS} 

concise: $(OBJ) Bell_Telescope_concise.o motion.o #add motion function
	${LNK} -o ${BIN} ${OBJ} Bell_Telescope_concise.o motion.o ${INC_FAST} ${CXXFLAGS} 

