import sys

for line in open(sys.argv[1]):
    dat = line.rstrip().split(',')

    if dat[1].endswith('_STD'):
        dat[1] = dat[1][:-4]
    print(','.join(dat))

