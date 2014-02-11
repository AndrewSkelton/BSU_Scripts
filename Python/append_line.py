import fileinput
import sys

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

for line in fileinput.input(['./merged.gtf'], inplace=True):
	info = line.split()
	var1 = info[0]
	if is_number(var1) == True:
		sys.stdout.write('chr{l}'.format(l=line))
	else:
		sys.stdout.write('{l}'.format(l=line))