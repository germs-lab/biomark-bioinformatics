import sys

#  python combine-biomark-csv.py *csv  > combined.csv

for n, f in enumerate(sys.argv[1:]):
    for m, line in enumerate(open(f)):
        dat = line.rstrip().split(',')
        dat.append(f.split()[0])
        print(','.join(dat))

    
