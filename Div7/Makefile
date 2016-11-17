all: serial enumspec1 enumspec2 purespec1 purespec2 purespec15 enumconv_mimd enumspec_mimd

serial: div_serial.cpp
	icc -mmic -std=c++11 -O3 -no-vec div_serial.cpp -o serial

enumspec1: div_enumspec1.cpp
	icc -mmic -std=c++11 -O3 div_enumspec1.cpp -o enumspec1

purespec1: div_purespec1.cpp
	icc -mmic -std=c++11 -O3 div_purespec1.cpp -o purespec1

enumspec2: div_enumspec2.cpp
	icc -mmic -std=c++11 -O3 div_enumspec2.cpp -o enumspec2

purespec2: div_purespec2.cpp
	icc -mmic -std=c++11 -O3 div_purespec2.cpp -o purespec2

purespec15: div_purespec15.cpp
	icc -mmic -std=c++11 -O3 div_purespec15.cpp -o purespec15

enumconv_mimd: div_enumconv_mimd.cpp
	icc -mmic -std=c++11 -O3 div_enumconv_mimd.cpp -openmp -o enumconv_mimd

enumspec_mimd: div_enumspec_mimd.cpp
	icc -mmic -std=c++11 -O3 div_enumspec_mimd.cpp -openmp -o enumspec_mimd

clean:
	rm serial enumspec1 enumspec2 purespec1 purespec2 purespec15 enumconv_mimd enumspec_mimd
