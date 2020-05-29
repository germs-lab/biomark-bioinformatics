# biomark-private

1.  Make sure all your data files are the only thing in the folder (e.g., no metadata)
2.  Run excel-to-csv.py:  python excel-to-csv.py *xlsx to convert all the excel files to csv
3.  Check to make sure each file has a corresponding csv file
4.  Combine these into one data file by running:  python combine-biomark-csv.py *csv (This is why you cannot have any other csvs besides what you want combined in this location)
