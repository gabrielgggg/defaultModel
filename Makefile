
tools = intel

ifeq ($(tools),intel)
  compiler = ifort
  flags = -Ofast -parallel -align -xHost -qopenmp -ipo -std18 -warn all,noexternal # -CB -g 
  libs = -mkl
else ifeq ($(tools),intelLLVM)
  compiler = ifx
  flags = -Ofast -align -xHost -qopenmp -ipo -std18 -g -warn all,noexternal 
  libs = -qmkl
else ifeq ($(tools),gnu)
  compiler = gfortran
  flags = -Ofast -march=native -mtune=native -pipe -fopenmp -flto -Wall -pedantic -g
  libs = -llapack -lblas -lm
else ifeq ($(tools),arm)
  compiler = armflang
  flags = -Ofast -fopenmp -mcpu=a64fx -flto -Wall
  libs = -armpl
else ifeq ($(tools),nvidia)
  compiler = nvfortran
  flags = -fast -Mconcur -Minfo -Mvect -mp=gpu -gpu=ccall 
  libs = 
else ifeq ($(tools),fujitsu)
  compiler = mpifrt
  flags = -Kparallel -Kopenmp -Kreduction -KSVE -Kfast -flto -fi 
  libs = -SSL2 
else ifeq ($(tools),cray)
  compiler = mpifort 
  flags = -O3 -h omp,scalar3,vector3,ipa3,contiguous,nobounds,wp,pl=./_pl
  libs =
else
  $(error Unknown toolchain?)
endif

.PHONY: all clean

all:	bin/defaultModel

bin/defaultModel:	src/defaultModel.f90 defMod.o
	$(compiler) $(flags) -o bin/defaultModel *.o src/defaultModel.f90 ${libs}

sim.o:	src/sim.f90 
	$(compiler) $(flags) -c src/sim.f90

NL.o:	src/NL.f90
	$(compiler) $(flags) -c src/NL.f90

defMod.o:	src/defMod.f90 NL.o sim.o
	$(compiler) $(flags) -c src/defMod.f90

clean:
	rm -f bin/* ./*.mod ./*.o

cleanresults:
	rm -f results/*.tab results/*.bin

cleanmlab:
	rm -f results/*.mat results/*.fig

ccc:	clean cleanresults

