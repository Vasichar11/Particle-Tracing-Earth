BIN=tracer   
CC=g++ 
LNK=g++ 
OPT=-O3
CXXFLAGS=-std=c++17 -lstdc++fs -Wall -g  -fopenmp


INC_HIVIFE= -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I/home/vchariton/HighFive/include/


INC_FAST= -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial  -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I/home/vchariton/HighFive/include/


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
	${CC} ${CXXFLAGS} ${OPT} -c struct_Telescope.cc    

adiabatic_motion.o: adiabatic_motion.cc
	${CC} ${CXXFLAGS} ${OPT} -c adiabatic_motion.cc

wpi.o: wpi.cc
	${CC} ${CXXFLAGS} ${OPT} -c wpi.cc

functions.o: functions.cc
	${CC} ${CXXFLAGS} ${OPT} -c functions.cc

time_rates.o: time_rates.cc
	${CC} ${CXXFLAGS} ${OPT} -c time_rates.cc



Bell_Telescope_performant.o: Bell_Telescope_performant.cc
	${CC} ${CXXFLAGS} ${OPT} -c Bell_Telescope_performant.cc    ${INC_FAST}
	
Bell_Telescope_adiabatic.o: Bell_Telescope_adiabatic.cc
	${CC} ${CXXFLAGS} ${OPT} -c Bell_Telescope_adiabatic.cc    ${INC_FAST}

Bell_Telescope_wpi.o: Bell_Telescope_wpi.cc
	${CC} ${CXXFLAGS} ${OPT} -c Bell_Telescope_wpi.cc    ${INC_FAST}



performant: $(OBJ) Bell_Telescope_performant.o
	${LNK} -o ${BIN} ${OBJ}  Bell_Telescope_performant.o ${INC_HIVIFE} ${CXXFLAGS} 

adiabatic: $(OBJ) Bell_Telescope_adiabatic.o adiabatic_motion.o #add adiabatic_motion function
	${LNK} -o ${BIN} ${OBJ} Bell_Telescope_adiabatic.o adiabatic_motion.o ${INC_HIVIFE} ${CXXFLAGS} 

wpi: $(OBJ) Bell_Telescope_wpi.o wpi.o #add wpi function
	${LNK} -o ${BIN} ${OBJ} Bell_Telescope_wpi.o wpi.o ${INC_HIVIFE} ${CXXFLAGS} 

