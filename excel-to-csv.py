import sys
import pandas

#Takes a list of files and converts it to csv.  python excel-to-csv.py *xlsx

for x in sys.argv[1:]:
    read_file = pandas.read_excel(x)
    new_csv_name = x[:-5] + ".csv"
    read_file.to_csv(new_csv_name, index=False)

