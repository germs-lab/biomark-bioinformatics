import sys

#  python combine-biomark-csv.py *csv  > combined.csv

for n, f in enumerate(sys.argv[1:]):
    if n == 1:
        for m, line in enumerate(open(f)):
            if m > 10:
                print(line.rstrip())
    else:
        for m, line in enumerate(open(f)):
            if m > 11:
                print(line.rstrip())
