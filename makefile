all: Disk_Math.h Polynomial.h test.cpp Monca.h
	g++ -g -Wall Disk_Math.h Polynomial.h Monca.h LA1.h test.cpp -I.
clean:
	
