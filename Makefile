objects = beamforming.o sacio.o

beamforming : $(objects)
	cc -o beamforming $(objects) -lm

$(objects) : sacio.h

clean :
	rm -f beamforming $(objects)
