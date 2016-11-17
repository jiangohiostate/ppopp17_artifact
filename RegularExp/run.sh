#!/bin/bash

echo "===============Test Regular Expression Matching==========="
echo "Build the programs....."
make
echo "Single-Core Test"
echo "@@@@@@Dataset: the first regular expression @@@@@@@@"
echo "######Dataset: 76.txt.utf-8#########"
echo "----------Run Serial-----------"
./serial datasets/regex1 datasets/76.txt.utf-8
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/regex1 datasets/76.txt.utf-8
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/regex1 datasets/76.txt.utf-8
echo "----------Run 3-way EnumSpec-----------"
./enumspec3 datasets/regex1 datasets/76.txt.utf-8
echo "----------Run 3-way PureSpec-----------"
./purespec3 datasets/regex1 datasets/76.txt.utf-8
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/regex1 datasets/76.txt.utf-8
echo ""
echo "######Dataset: 50247.txt.utf-8#########"
echo "----------Run Serial-----------"
./serial datasets/regex1 datasets/50247.txt.utf-8
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/regex1 datasets/50247.txt.utf-8
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/regex1 datasets/50247.txt.utf-8
echo "----------Run 3-way EnumSpec-----------"
./enumspec3 datasets/regex1 datasets/50247.txt.utf-8
echo "----------Run 3-way PureSpec-----------"
./purespec3 datasets/regex1 datasets/50247.txt.utf-8
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/regex1 datasets/50247.txt.utf-8
echo ""
echo "######Dataset: 98.txt.utf-8#########"
echo "----------Run Serial-----------"
./serial datasets/regex1 datasets/98.txt.utf-8
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/regex1 datasets/98.txt.utf-8
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/regex1 datasets/98.txt.utf-8
echo "----------Run 3-way EnumSpec-----------"
./enumspec3 datasets/regex1 datasets/98.txt.utf-8
echo "----------Run 3-way PureSpec-----------"
./purespec3 datasets/regex1 datasets/98.txt.utf-8
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/regex1 datasets/98.txt.utf-8
echo ""
echo "######Dataset: 74.txt.utf-8#########"
echo "----------Run Serial-----------"
./serial datasets/regex1 datasets/74.txt.utf-8
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/regex1 datasets/74.txt.utf-8
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/regex1 datasets/74.txt.utf-8
echo "----------Run 3-way EnumSpec-----------"
./enumspec3 datasets/regex1 datasets/74.txt.utf-8
echo "----------Run 3-way PureSpec-----------"
./purespec3 datasets/regex1 datasets/74.txt.utf-8
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/regex1 datasets/74.txt.utf-8

echo ""
echo ""
echo "@@@@@@Dataset: the second regular expression @@@@@@@@"
echo "######Dataset: 76.txt.utf-8#########"
echo "----------Run Serial-----------"
./serial datasets/regex2 datasets/76.txt.utf-8
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/regex2 datasets/76.txt.utf-8
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/regex2 datasets/76.txt.utf-8
echo "----------Run 3-way EnumSpec-----------"
./enumspec3 datasets/regex2 datasets/76.txt.utf-8
echo "----------Run 3-way PureSpec-----------"
./purespec3 datasets/regex2 datasets/76.txt.utf-8
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/regex2 datasets/76.txt.utf-8
echo ""
echo "######Dataset: 50247.txt.utf-8#########"
echo "----------Run Serial-----------"
./serial datasets/regex2 datasets/50247.txt.utf-8
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/regex2 datasets/50247.txt.utf-8
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/regex2 datasets/50247.txt.utf-8
echo "----------Run 3-way EnumSpec-----------"
./enumspec3 datasets/regex2 datasets/50247.txt.utf-8
echo "----------Run 3-way PureSpec-----------"
./purespec3 datasets/regex2 datasets/50247.txt.utf-8
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/regex2 datasets/50247.txt.utf-8
echo ""
echo "######Dataset: 98.txt.utf-8#########"
echo "----------Run Serial-----------"
./serial datasets/regex2 datasets/98.txt.utf-8
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/regex2 datasets/98.txt.utf-8
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/regex2 datasets/98.txt.utf-8
echo "----------Run 3-way EnumSpec-----------"
./enumspec3 datasets/regex2 datasets/98.txt.utf-8
echo "----------Run 3-way PureSpec-----------"
./purespec3 datasets/regex2 datasets/98.txt.utf-8
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/regex2 datasets/98.txt.utf-8
echo ""
echo "######Dataset: 74.txt.utf-8#########"
echo "----------Run Serial-----------"
./serial datasets/regex2 datasets/74.txt.utf-8
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/regex2 datasets/74.txt.utf-8
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/regex2 datasets/74.txt.utf-8
echo "----------Run 3-way EnumSpec-----------"
./enumspec3 datasets/regex2 datasets/74.txt.utf-8
echo "----------Run 3-way PureSpec-----------"
./purespec3 datasets/regex2 datasets/74.txt.utf-8
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/regex2 datasets/74.txt.utf-8

echo ""
echo ""
echo ""
echo "Multicore Test"
echo "Generate a 1GB textfile randomly"
cd datasets
python gen_data.py
cd ..
echo "######Dataset: big.txt#########"
echo "@@@@@@Dataset: the first regular expression @@@@@@@@"
echo "----------Run Serial-----------"
./serial datasets/regex1 datasets/big.txt
echo "----------Run EnumSpec-----------"
for i in 1 2 4 8 16 32 60; do
  echo "%%%%%number of threads = $i"
  ./enumspec_mimd datasets/regex1 datasets/big.txt $i
done
echo "----------Run EnumConv-----------"
for i in 1 2 4 8 16 32 60; do
  echo "%%%%%number of threads = $i"
  ./enumconv_mimd datasets/regex1 datasets/big.txt $i
done
echo ""
echo "######Dataset: big.txt#########"
echo "@@@@@@Dataset: the second regular expression @@@@@@@@"
echo "----------Run Serial-----------"
./serial datasets/regex2 datasets/big.txt
echo "----------Run EnumSpec-----------"
for i in 1 2 4 8 16 32 60; do
  echo "%%%%%number of threads = $i"
  ./enumspec_mimd datasets/regex2 datasets/big.txt $i
done
echo "----------Run EnumConv-----------"
for i in 1 2 4 8 16 32 60; do
  echo "%%%%%number of threads = $i"
  ./enumconv_mimd datasets/regex2 datasets/big.txt $i
done

