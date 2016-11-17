#!/bin/bash

echo "===============Test HTML Tokenization==========="
echo "Build the programs....."
make
echo "Single-Core Test"
echo "######Dataset: webpage1 #########"
echo "----------Run Serial-----------"
./serial datasets/nytime1.html
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/nytime1.html
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/nytime1.html
echo "----------Run 3-way EnumSpec-----------"
./enumspec3 datasets/nytime1.html
echo "----------Run 3-way PureSpec-----------"
./purespec3 datasets/nytime1.html
echo "----------Run 7-way EnumSpec-----------"
./enumspec7 datasets/nytime1.html
echo "----------Run 7-way PureSpec-----------"
./purespec7 datasets/nytime1.html
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/nytime1.html
echo ""
echo "######Dataset: webpage2 #########"
echo "----------Run Serial-----------"
./serial datasets/nytime2.html
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/nytime2.html
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/nytime2.html
echo "----------Run 3-way EnumSpec-----------"
./enumspec3 datasets/nytime2.html
echo "----------Run 3-way PureSpec-----------"
./purespec3 datasets/nytime2.html
echo "----------Run 7-way EnumSpec-----------"
./enumspec7 datasets/nytime2.html
echo "----------Run 7-way PureSpec-----------"
./purespec7 datasets/nytime2.html
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/nytime2.html
echo ""
echo "######Dataset: webpage3 #########"
echo "----------Run Serial-----------"
./serial datasets/nytime3.html
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/nytime3.html
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/nytime3.html
echo "----------Run 3-way EnumSpec-----------"
./enumspec3 datasets/nytime3.html
echo "----------Run 3-way PureSpec-----------"
./purespec3 datasets/nytime3.html
echo "----------Run 7-way EnumSpec-----------"
./enumspec7 datasets/nytime3.html
echo "----------Run 7-way PureSpec-----------"
./purespec7 datasets/nytime3.html
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/nytime3.html
echo ""
echo "######Dataset: webpage4 #########"
echo "----------Run Serial-----------"
./serial datasets/nytime4.html
echo "----------Run 1-way EnumSpec-----------"
./enumspec1 datasets/nytime4.html
echo "----------Run 1-way PureSpec-----------"
./purespec1 datasets/nytime4.html
echo "----------Run 3-way EnumSpec-----------"
./enumspec3 datasets/nytime4.html
echo "----------Run 3-way PureSpec-----------"
./purespec3 datasets/nytime4.html
echo "----------Run 7-way EnumSpec-----------"
./enumspec7 datasets/nytime4.html
echo "----------Run 7-way PureSpec-----------"
./purespec7 datasets/nytime4.html
echo "----------Run 15-way PureSpec-----------"
./purespec15 datasets/nytime4.html

