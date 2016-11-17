#!/bin/bash

echo "===============Test Div7==========="
echo "Build the programs....."
make
echo "Single-Core Test"
echo "######Dataset: number1 #########"
echo "----------Run Serial-----------"
./serial datasets/number1
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/number1
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/number1
echo "----------Run 2-way EnumSpec-----------"
./enumspec2 datasets/number1
echo "----------Run 2-way PureSpec-----------"
./purespec2 datasets/number1
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/number1
echo ""
echo "######Dataset: number2 #########"
echo "----------Run Serial-----------"
./serial datasets/number2
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/number2
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/number2
echo "----------Run 2-way EnumSpec-----------"
./enumspec2 datasets/number2
echo "----------Run 2-way PureSpec-----------"
./purespec2 datasets/number2
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/number2
echo ""
echo "######Dataset: number3 #########"
echo "----------Run Serial-----------"
./serial datasets/number3
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/number3
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/number3
echo "----------Run 2-way EnumSpec-----------"
./enumspec2 datasets/number3
echo "----------Run 2-way PureSpec-----------"
./purespec2 datasets/number3
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/number3
echo ""
echo "######Dataset: number4 #########"
echo "----------Run Serial-----------"
./serial datasets/number4
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/number4
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/number4
echo "----------Run 2-way EnumSpec-----------"
./enumspec2 datasets/number4
echo "----------Run 2-way PureSpec-----------"
./purespec2 datasets/number4
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/number4
echo ""


echo ""
echo "Multicore Test"
echo "Generate a longer binary number randomly"
cd datasets
python gen_data.py
cd ..
echo "######Dataset: number_longer #########"
echo "----------Run Serial-----------"
./serial datasets/number_longer
echo "----------Run EnumSpec-----------"
for i in 1 2 4 8 16 32 60; do
  echo "%%%%%number of threads = $i"
  ./enumspec_mimd datasets/number_longer $i
done
echo "----------Run EnumConv-----------"
for i in 1 2 4 8 16 32 60; do
  echo "%%%%%number of threads = $i"
  ./enumconv_mimd datasets/number_longer $i
done

