#!/bin/bash

echo "===============Test Huffman Decoding==========="
echo "Build the programs....."
make
echo "Single-Core Test"
echo "######Dataset: 76.txt.utf-8#########"
echo "----------Run Serial-----------"
./serial datasets/76.txt.utf-8
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/76.txt.utf-8
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/76.txt.utf-8
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/76.txt.utf-8
echo ""
echo "######Dataset: 50247.txt.utf-8#########"
echo "----------Run Serial-----------"
./serial datasets/50247.txt.utf-8
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/50247.txt.utf-8
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/50247.txt.utf-8
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/50247.txt.utf-8
echo ""
echo "######Dataset: 98.txt.utf-8#########"
echo "----------Run Serial-----------"
./serial datasets/98.txt.utf-8
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/98.txt.utf-8
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/98.txt.utf-8
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/98.txt.utf-8
echo ""
echo "######Dataset: 74.txt.utf-8#########"
echo "----------Run Serial-----------"
./serial datasets/74.txt.utf-8
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/74.txt.utf-8
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/74.txt.utf-8
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/74.txt.utf-8

echo ""
echo "Multicore Test"
echo "Generate a 100MB textfile randomly"
cd datasets
python gen_data.py
cd ..
echo "######Dataset: big.txt#########"
echo "----------Run Serial-----------"
./serial datasets/big.txt
echo "----------Run EnumSpec-----------"
for i in 1 2 4 8 16 32 60; do
  echo "%%%%%number of threads = $i"
  ./enumspec_mimd datasets/big.txt $i
done
echo "----------Run EnumConv-----------"
for i in 1 2 4 8 16 32 60; do
  echo "%%%%%number of threads = $i"
  ./enumconv_mimd datasets/big.txt $i
done

