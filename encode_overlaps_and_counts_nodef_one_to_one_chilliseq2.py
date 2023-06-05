from pybedtools import BedTool
import os,sys
path = "./"
sample_list1 = open(''.join(path + sys.argv[1]), "r")
sample_list2 = open(''.join(path + sys.argv[2]), "r")
#sample_list = open(''.join(path + sys.argv[1]))
sample_list1 = list(sample_list1)
sample_list2 = list(sample_list2)
#print(sample_list)
if os.path.exists('One_to_one_encode_overlaps_counts_out_chilli2.txt'):
    os.remove('One_to_one_encode_overlaps_counts_out_chilli2.txt')
for list1 in sample_list1:
    print(list1)
    for list2 in sample_list2:
        print(list2)
        overlaps = {}
        counts = {}
        mybed_peak = ''.join(path + list1).strip()
        mybed_peak_id = mybed_peak.split('/')[-1:][0].replace('_peaks_id.bed','')
        #print(mybed_peak_id)
        bedfile1 = BedTool(mybed_peak)
        counts[mybed_peak_id] = len(bedfile1)
        #print(type(bedfile1))
        mybed1_peak = ''.join(path + list2).strip()
        mybed1_peak_id = mybed1_peak.split('/')[-1:][0].replace('_peaks_id.bed','')
        #print(mybed1_peak_id)
        bedfile2 = BedTool(mybed1_peak)
        counts[mybed1_peak_id] = len(bedfile2)
        #print(type(bedfile2))
        intersect_file = bedfile1.intersect(bedfile2, wa=True, wb=True)
        for intcont in intersect_file:
            # first peaks ids dictionaries
            if intcont[11] not in overlaps:
                overlaps[intcont[11]] = {}
            # dictionaries of second peaks ids intersected with first peaks
            if intcont[23] not in overlaps[intcont[11]]:
                overlaps[intcont[11]][intcont[23]] = set()
            overlaps[intcont[11]][intcont[23]].add(intcont[3])
    
            # second peaks ids dictionaries
            if intcont[23] not in overlaps:
                overlaps[intcont[23]] = {}
            # dictionaries of first peaks ids intersected with second peaks
            if intcont[11] not in overlaps[intcont[23]]:
                overlaps[intcont[23]][intcont[11]] = set()
            overlaps[intcont[23]][intcont[11]].add(intcont[15])
        print(len(overlaps))
        if len(overlaps) == 0:
            print(''.join(mybed_peak_id + "\t" + mybed1_peak_id + "\t" + "0" + "\t" + str(len(bedfile1)) + "\t" + "0" + "\n" + mybed1_peak_id + "\t" + mybed_peak_id + "\t" + "0" + "\t" + str(len(bedfile2)) + "\t" + "0"))
            with open("One_to_one_encode_overlaps_counts_out_chilli2.txt" , "a") as myfinalout:
                myfinalout.write(''.join(mybed_peak_id + "\t" + mybed1_peak_id + "\t" + "0" + "\t" + str(len(bedfile1)) + "\t" + "0" + "\n" + mybed1_peak_id + "\t" + mybed_peak_id + "\t" + "0" + "\t" + str(len(bedfile2)) + "\t" + "0" + "\n"))
                myfinalout.close()
        else:
            for i in overlaps:
                for j in overlaps[i]:
                        print(''.join(i + "\t" + j + "\t" + str(len(overlaps[i][j])) + "\t" + str(counts[i]) + "\t" + str(float(int(len(overlaps[i][j]))) / float(int(counts[i])))))
                        with open("One_to_one_encode_overlaps_counts_out_chilli2.txt" , "a") as myfinalout:
                            myfinalout.write(''.join(i + "\t" + j + "\t" + str(len(overlaps[i][j])) + "\t" + str(counts[i]) + "\t" + str(float(int(len(overlaps[i][j]))) / float(int(counts[i]))) + "\n"))
                            myfinalout.close()


#awk '{print $1"%"$2"%"$3"%"$4"%"$5}' One_to_one_encode_overlaps_counts_out_chilli2.txt | sort -k1,1 -u | awk -F'%' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' > One_to_one_encode_overlaps_counts_out_chilli2_dedup.txt

print("Intersection and Proportion Analysis Finished !")

print("Must run deduplication command of awk as in this script")
