BIN = tracer_bell_nowpi  tracer_bell_wpi    ray_int  tracer_li_wpi 
CC  = g++ 
LNK = g++ 
OPT = -O3
CXXFLAGS = -std=c++17 -lstdc++fs -Wall -g -fopenmp


INC_HIFIVE = -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I/home/billy/HighFive/include/

INC_H5LESS = -I/usr/include/hdf5/serial -Wdate-time -D_FORTIFY_SOURCE=2 -g  -fdebug-prefix-map=/build/hdf5-X9JKIg/hdf5-1.10.0-patch1+docs=. -fstack-protector-strong -Wformat -Werror=format-security -L/usr/lib/x86_64-linux-gnu/hdf5/serial  -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial -I/home/billy/HighFive/include/


OBJ_BELL =  functions.o   time_rates.o   struct_Particles.o   struct_Species.o   struct_Telescope.o  

OBJ_LI =  Li_Ray_Telescope.o functions.o   time_rates.o   struct_Particles.o    struct_Telescope.o   is_in_packet.o   li_wpi.o   read_hdf5.o

OBJ_INT = Read_Ray_Write.o   read_csv.o   interpolate.o



#RULES
default:tracer_bell_nowpi
li_ray_wpi: ray_int   tracer_li_wpi
all:${BIN}
clean: 
	rm  *.o 
allclean:
	rm *.o ${BIN} *.out *.png


#For all code
struct_Particles.o: struct_Particles.cc
	${CC} ${CXXFLAGS} ${OPT} -c struct_Particles.cc   
struct_Species.o: struct_Species.cc 	#Only for Bell code
	${CC} ${CXXFLAGS} ${OPT} -c struct_Species.cc 
struct_Telescope.o: struct_Telescope.cc
	${CC} ${CXXFLAGS} ${OPT} -c struct_Telescope.cc    
functions.o: functions.cc
	${CC} ${CXXFLAGS} ${OPT} -c functions.cc
time_rates.o: time_rates.cc
	${CC} ${CXXFLAGS} ${OPT} -c time_rates.cc


#For Ray Interpolation
read_hdf5.o: read_hdf5.cc
	${CC} ${CXXFLAGS} ${OPT} -c read_hdf5.cc  ${INC_H5LESS}
read_csv.o: read_csv.cc
	${CC} ${CXXFLAGS} ${OPT} -c read_csv.cc
interpolate.o: interpolate.cc
	${CC} ${CXXFLAGS} ${OPT} -c interpolate.cc
Read_Ray_Write.o: Read_Ray_Write.cc
	${CC} ${CXXFLAGS} ${OPT} -c Read_Ray_Write.cc  ${INC_H5LESS}



bell_nowpi.o: bell_nowpi.cc
	${CC} ${CXXFLAGS} ${OPT} -c bell_nowpi.cc
bell_wpi.o: bell_wpi.cc
	${CC} ${CXXFLAGS} ${OPT} -c bell_wpi.cc
li_wpi.o: li_wpi.cc
	${CC} ${CXXFLAGS} ${OPT} -c li_wpi.cc ${INC_H5LESS}

Li_Ray_Telescope.o: Li_Ray_Telescope.cc
	${CC} ${CXXFLAGS} ${OPT} -c Li_Ray_Telescope.cc  ${INC_H5LESS}
Bell_Telescope_adiabatic.o: Bell_Telescope_adiabatic.cc
	${CC} ${CXXFLAGS} ${OPT} -c Bell_Telescope_adiabatic.cc    ${INC_H5LESS}
Bell_Telescope_wpi.o: Bell_Telescope_wpi.cc
	${CC} ${CXXFLAGS} ${OPT} -c Bell_Telescope_wpi.cc    ${INC_H5LESS}








#BELL ADIABATIC OR WPI
tracer_bell_nowpi: $(OBJ_BELL) Bell_Telescope_adiabatic.o bell_nowpi.o #add bell_nowpi function
	${LNK} -o tracer_bell_nowpi ${OBJ_BELL} Bell_Telescope_adiabatic.o bell_nowpi.o ${INC_HIFIVE} ${CXXFLAGS} 

tracer_bell_wpi: $(OBJ_BELL) Bell_Telescope_wpi.o bell_wpi.o #add wpi function
	${LNK} -o tracer_bell_wpi ${OBJ_BELL} Bell_Telescope_wpi.o bell_wpi.o ${INC_HIFIVE} ${CXXFLAGS} 


#LI RAY TRACING AND WPI
ray_int:  ${OBJ_INT}
	${LNK} -o ray_int ${OBJ_INT} ${INC_HIFIVE} ${CXXFLAGS}

tracer_li_wpi: ${OBJ_LI}
	${LNK} -o tracer_li_wpi ${OBJ_LI} ${INC_HIFIVE} ${CXXFLAGS} ${OMP}
