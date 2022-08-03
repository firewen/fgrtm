FC=ifort
#FC=gfortran
vpath %.f src
vpath %.F src
vpath %.f90 src
vpath %.o OBJ

OBJS = MainProgram.o bessel.o comvar.o dwim.o fft.o Gbasic.o Ginput.o Grtcoefs.o G-s-r.o\
	 Gsubs.o integrang_calc.o MTX_SUB.o iptam.o SOURCE.o UKO.o Ydumtx.o write2sac.o \
	 raytracing.o directarrival.o

OBJ_PATH = OBJ
SRC_PATH = src
BIN_PATH = bin

Flags = -L${SACHOME}/lib -lsacio -I${SACHOME}/include -lsac

$(BIN_PATH)/MainProgram : $(OBJS)
	$(FC) $^ -o $@ ${Flags}

$(OBJ_PATH)/MainProgram.o : MainProgram.f
	$(FC) -c $< -o $@

$(OBJ_PATH)/bessel.o : bessel.f
	$(FC) -c $< -o $@

$(OBJ_PATH)/comvar.o : comvar.f90
	$(FC) -c $< -o $@

$(OBJ_PATH)/dwim.o : dwim.f90
	$(FC) -c $< -o $@

$(OBJ_PATH)/fft.o : fft.f
	$(FC) -c $< -o $@

$(OBJ_PATH)/Gbasic.o : Gbasic.F
	$(FC) -c $< -o $@

$(OBJ_PATH)/Ginput.o : Ginput.F
	$(FC) -c $< -o $@

$(OBJ_PATH)/Grtcoefs.o : grtcoefs.f90
	$(FC) -c $< -o $@

$(OBJ_PATH)/G-s-r.o : G-s-r.f
	$(FC) -c $< -o $@

$(OBJ_PATH)/Gsubs.o : Gsubs.f
	$(FC) -c $< -o $@

$(OBJ_PATH)/integrang_calc.o : integrang_calc.f90
	$(FC) -c $< -o $@

$(OBJ_PATH)/MTX_SUB.o : MTX_SUB.F
	$(FC) -c $< -o $@

$(OBJ_PATH)/iptam.o : iptam.f90
	$(FC) -c $< -o $@

$(OBJ_PATH)/SOURCE.o : source.f90
	$(FC) -c $< -o $@

$(OBJ_PATH)/UKO.o : UKO.F
	$(FC) -c $< -o $@

$(OBJ_PATH)/Ydumtx.o  : ydumtx.f90
	$(FC) -c $< -o $@

$(OBJ_PATH)/write2sac.o : write2sac.f90
	$(FC) -c $< -o $@ ${Flags}

$(OBJ_PATH)/raytracing.o : raytracing.f90
	$(FC) -c $< -o $@

$(OBJ_PATH)/directarrival.o : directarrival.f90
	$(FC) -c $< -o $@

clean:
	rm OBJ/*.o bin/MainProgram
