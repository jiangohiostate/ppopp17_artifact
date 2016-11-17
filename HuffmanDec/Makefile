all: serial enumspec1 purespec1 purespec15 enumspec_mimd enumconv_mimd

serial: huffman_serial.cpp
	icc -mmic -std=c++11 -O3 -no-vec huffman_serial.cpp -o serial

enumspec1: huffman_enumspec1.cpp
	icc -mmic -std=c++11 -O3 huffman_enumspec1.cpp -o enumspec1

purespec1: huffman_purespec1.cpp
	icc -mmic -std=c++11 -O3 huffman_purespec1.cpp -o purespec1

purespec15: huffman_purespec15.cpp
	icc -mmic -std=c++11 -O3 huffman_purespec15.cpp -o purespec15

enumspec_mimd: huffman_enumspec_mimd.cpp
	icc -mmic -std=c++11 -O3 huffman_enumspec_mimd.cpp -openmp -o enumspec_mimd

enumconv_mimd: huffman_enumconv_mimd.cpp
	icc -mmic -std=c++11 -O3 huffman_enumconv_mimd.cpp -openmp -o enumconv_mimd

clean:
	rm serial enumspec1 purespec1 purespec15 enumspec_mimd enumconv_mimd
