import fileinput
import sys

for line in fileinput.input(['./merged.gtf'], inplace=True):
    sys.stdout.write('chr{l}'.format(l=line))
