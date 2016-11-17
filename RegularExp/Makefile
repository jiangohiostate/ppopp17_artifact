all: serial enumspec1 purespec1 enumspec3 purespec3 purespec15 enumconv_mimd enumspec_mimd 

serial: rex_serial.cpp
	icc -mmic -std=c++11 -O3 -no-vec rex_serial.cpp -o serial

enumspec1: rex_enumspec1.cpp
	icc -mmic -std=c++11 -O3 rex_enumspec1.cpp -o enumspec1

purespec1: rex_purespec1.cpp
	icc -mmic -std=c++11 -O3 rex_purespec1.cpp -o purespec1

enumspec3: rex_enumspec3.cpp
	icc -mmic -std=c++11 -O3 rex_enumspec3.cpp -o enumspec3

purespec3: rex_purespec3.cpp
	icc -mmic -std=c++11 -O3 rex_purespec3.cpp -o purespec3

purespec15: rex_purespec15.cpp
	icc -mmic -std=c++11 -O3 rex_purespec15.cpp -o purespec15

enumconv_mimd: rex_enumconv_mimd.cpp
	icc -mmic -std=c++11 -O3 rex_enumconv_mimd.cpp -openmp -o enumconv_mimd

enumspec_mimd: rex_enumspec_mimd.cpp
	icc -mmic -std=c++11 -O3 rex_enumspec_mimd.cpp -openmp -o enumspec_mimd

clean:
	rm serial enumspec1 purespec1 enumspec3 purespec3 purespec15 enumconv_mimd enumspec_mimd 
