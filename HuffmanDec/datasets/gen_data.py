import random
import sys

outfile = open('big.txt', 'w')
for i in range(100*1024*1024):
    outfile.write(random.choice('abcdefghijklmnopqrstubwxyz '))
    if(random.random() < 0.01):
        outfile.write('\n')
