import random
import sys

outfile = open('big.txt', 'w')
for j in range(1024):
    str = ''
    for i in range(1024*1024):
        str += (random.choice('abcdefghijklmnopqrstubwxyz '))
        if(random.random() < 0.01):
            str += ('\n')
    outfile.write(str);
