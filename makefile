srt: Hybrid.cpp Hybrid.h
	mpic++ -fopenmp -c Hybrid.cpp Hybrid.h -o hy.o -I.
	ar -rsc srt.a hy.o