import os
import sys
import pandas as pd
from sortedcontainers import SortedDict
from pybedtools import BedTool

# Set the path of your directory
path = "/mnt/home3/reid/av638/ENCODE/"
binsize = 1000
# Peak file should have chr at col1 , start at col2, end at col3, MACS score at col7. This format you generally obtain from MACS2 callpeak.
sample_list1 = open(sys.argv[1], "r")
sample_list1 = list(sample_list1)
# print(sample_list1)
fulldict = dict()
# Read the blacklisted regions
blacklist = ''.join(path + sys.argv[2]).strip()
blacklist_name = "blacklist"
# print(blacklist_name)
blacklist_bed = BedTool(blacklist)
# print(blacklist_bed)

if os.path.exists(sys.argv[3]):
    os.remove(sys.argv[3])
for samples in sample_list1:
    samples = samples.strip('').split('\n')
    # print(samples[0])
    peakfile = open(''.join('' + samples[0]), 'r')
    peakfile_name = os.path.basename(peakfile.name).split("_")[0]
    peaks = BedTool(''.join('' + samples[0]))
    # print(peaks)
    # filter the object by column 6 (macs score) > 100
    # get the number of intervals in the filtered object
    peaks_filtered = peaks.filter(lambda x: float(x[6]) >= 100)
    peaks_free_of_blacklist = peaks_filtered.intersect(blacklist_bed, v=True)
    peaks_filtered_count = peaks_free_of_blacklist.count()
    # print(peaks_free_of_blacklist)
    # print(peaks_filtered_count)
    # Create and empty dictionary to store bins
    bindict = dict()
    bindict[peakfile_name] = {}
    # bindict[peakfile_name][peakfile_name] = 1
    # print(bindict)
    # Apply condition if the file is not empty after intersection
    if peaks_filtered_count > 0:
        # print("Filtered Bedtools object is not empty")
        # Intersect the filtered peaks with the blacklisted regions
        # Create Bins and store in dictionary with key value 1
        for peak_lines in peaks_free_of_blacklist:
            # print(peak_lines)
            chrom = peak_lines[0]
            # Divide by binsize
            newstart = (str(int(peak_lines[1]) / binsize)).split('.')[0]
            # print(newstart)
            newend = (str(int(peak_lines[2]) / binsize)).split('.')[0]
            # print(chrom, newstart, newend)
            # If bin is short that is bin/1kb Start = End
            if int(newstart) == int(newend):
                # print(chrom,"Start:",peak_lines[1], "End:",peak_lines[2], "BinStart:",newstart, "BinEnd: ",newend, "No")
                # Store Start
                binkey = ''.join(chrom + '_' + newstart)
                bindict[peakfile_name][binkey] = 1
                # Store End
                binkey = ''.join(chrom + '_' + newend)
                bindict[peakfile_name][binkey] = 1
                # print(bindict)
            # If bin is short that is bin/1kb Start < End
            elif (int(newend) - int(newstart)) >= 1:
                # print(peakfile_name, "Start:",lines[1], "End:",lines[2], "BinStart:",newstart, "BinEnd: ",newend, "Yes")
                # print(int(newend) - int(newstart))
                for hiddenbins in range(int(newstart), int(newend), 1):
                    # print(hiddenbins)
                    # Store Start
                    binkey = ''.join(chrom + '_' + newstart)
                    bindict[peakfile_name][binkey] = 1
                    # Store Bins in-between start and end
                    binkey = ''.join(chrom + '_' + str(hiddenbins))
                    bindict[peakfile_name][binkey] = 1
                    # Store End
                    binkey = ''.join(chrom + '_' + newend)
                    bindict[peakfile_name][binkey] = 1
                    # print(bindict)

    # print(bindict)
    # print(type(fulldict))
    if peakfile_name not in fulldict:
        fulldict.update(bindict)
    elif peakfile_name in fulldict:
        # print(bindict[peakfile_name])
        # If key exist add values to same key
        fulldict[peakfile_name].update(bindict[peakfile_name])
        # fulldict.update(bindict)
    # Put all data in new dictionary
    # Iterate through lines and write them in the output file
    for i in SortedDict(bindict):
        for j in SortedDict(bindict[i]):
            # print(i, j)
            print(''.join(i + '\t' + j + '\t' + str(bindict[i][j])))
            with open(sys.argv[3], 'a') as myfile:
                myfile.write(''.join(i + '\t' + j + '\t' + str(bindict[i][j])) + '\n')
                myfile.close()
    else:
        pass

# print(fulldict)
# Make a dataframe from full dictionary
fulldictdf = pd.DataFrame(fulldict)
print(fulldictdf.to_string(na_rep="0"))

if os.path.exists(sys.argv[4]):
    os.remove(sys.argv[4])
with open(sys.argv[4], 'a') as yourfile:
    yourfile.write(fulldictdf.to_string(na_rep="0"))
    yourfile.close()
