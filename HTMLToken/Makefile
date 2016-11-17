all: serial enumspec1 enumspec3 enumspec7 purespec1 purespec3 purespec7 purespec15

serial: html_serial.cpp
	icc -mmic -std=c++11 -O3 -no-vec html_serial.cpp -o serial

enumspec1: html_enumspec1.cpp
	icc -mmic -std=c++11 -O3 html_enumspec1.cpp -o enumspec1

purespec1: html_purespec1.cpp
	icc -mmic -std=c++11 -O3 html_purespec1.cpp -o purespec1

enumspec3: html_enumspec3.cpp
	icc -mmic -std=c++11 -O3 html_enumspec3.cpp -o enumspec3

purespec3: html_purespec3.cpp
	icc -mmic -std=c++11 -O3 html_purespec3.cpp -o purespec3

enumspec7: html_enumspec7.cpp
	icc -mmic -std=c++11 -O3 html_enumspec7.cpp -o enumspec7

purespec7: html_purespec7.cpp
	icc -mmic -std=c++11 -O3 html_purespec7.cpp -o purespec7

purespec15: html_purespec15.cpp
	icc -mmic -std=c++11 -O3 html_purespec15.cpp -o purespec15

clean:
	rm serial enumspec1 enumspec3 enumspec7 purespec1 purespec3 purespec7 purespec15
