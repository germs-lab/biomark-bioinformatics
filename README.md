# biomark-private

1.  Make sure all your data files are the only thing in the folder (e.g., no metadata)
2.  Run excel-to-csv.py:  python excel-to-csv.py *xlsx to convert all the excel files to csv (See Step 4 also)
3.  Check to make sure each file has a corresponding csv file
4.  Open each csv file and remove all the header files (Alternately, you can do this before step 2 as excel files).  Note, I think biomark creates these headers differently for each run, so I just remove them manually
5.  Combine these into one data file by running:  python combine-biomark-csv.py *csv (This is why you cannot have any other csvs besides what you want combined in this location)
6.  Metadata:  Each row should be unique in the metadata, an example is shown in the metadata folder.  Your sample name header should be "Name".  You should have as many rows in your metadata (besides the header) as unique Name in your data table.
7.  Standards:  Concentrations of standards should be in the biomark output, or at least that is how this program works, it assumes your rConc is accurate in your data table

