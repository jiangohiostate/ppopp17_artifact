import random
from sys import argv

outfile = open('number_longer', 'w')
for j in range(1024):
    str = ''
    for i in range(1024*1024):
      if random.random() < 0.5:
        str += '0'
      else:
        str += '1'
    outfile.write(str)

outfile.close()

