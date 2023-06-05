import os
import sys
from operator import itemgetter
from collections import Counter
from collections import defaultdict
import argparse
import re
print("_________________________________________________________________________________________\n")
print('''
################################
Complex Mapper:version 5.0
################################

--> For only Protein eg. EZH2
Usage: python complex_mapper_v5.py --summitsc summit_mydata_p1_test_re.txt --complex humanComplexes.txt --chip SUZ12 --top 1

--> For Protein and Identifier eg. EZH2_GSM1427070_EZH2-F_rep1
Usage:  python complex_mapper_v5.py --summitsc summit_mydata_p1_test_re.txt --complex humanComplexes.txt --chip EZH2_GSM1427070_EZH2-F_rep1 --top 1
''')
print("_________________________________________________________________________________________\n")
#grep CRAMP1_NTC_IgG_NTC summit_mydata_p1_re.txt | awk '{if($3 >= 5) print $0}' | sort -g -r -k3,3 | head -n 10 | awk -F'_' '{print $1}' | sort -k1,1 -u
parser = argparse.ArgumentParser(prog='complex_mapper_v5.py', description='Complex Mapper:version 5.0: Map coassociation score to protein complexes')
parser.add_argument('--summitsc', help='summit output with co-association score from R script', type=str)
parser.add_argument('--complex', help='file containing human complexes in specific format', type=str)
parser.add_argument('--chip', help='Your choice of protein', default=None, type=str)
parser.add_argument('--top', help='Range', default=None, type=int)
# Print help message if no arguments
if len(sys.argv)==1:
   parser.print_help(sys.stderr)
   sys.exit(1)

args = parser.parse_args()

'''Parse summit coassociation score file'''
summito = open(args.summitsc, "r")
summito = summito.read().strip().split('\n')
overlapslist = []
for overlaps in summito:
    overlaps = overlaps.strip().split("\t")
    #print(overlaps)
    #print(len(list(args.chip.split('_'))))
    if len(list(args.chip.split('_'))) > 1:
        myinputTF = args.chip
        #print(myinputTF)
        '''Remove extra encode IDs / identifier after "_" so that I can map all replicates of that protein with data eg. SUZ12 will map to both SUZ12_ENXCCGXCC and SUZ12_ENCVAJCBC now'''
        inTFleft = overlaps[0]
        inTFright = overlaps[1]
        #print(inTFleft, inTFright)
        '''Keep the IDs attached to protein when I prepared the list'''
        if inTFleft == args.chip:
            if float(overlaps[2]) >= float(5):
                overlapslist.append(overlaps)
        elif inTFright == args.chip:
            if float(overlaps[2]) >= float(5):
                overlapslist.append(overlaps)
    else:
        myinputTF = args.chip
        #print(myinputTF)
        '''Remove extra encode IDs / identifier after "_" so that I can map all replicates of that protein with data eg. SUZ12 will map to both SUZ12_ENXCCGXCC and SUZ12_ENCVAJCBC now'''
        inTFleft = overlaps[0].split('_')[0]
        inTFright = overlaps[1].split('_')[0]
        '''Keep the IDs attached to protein when I prepared the list'''
        if inTFleft == args.chip:
            if float(overlaps[2]) >= float(5):
                overlapslist.append(overlaps)
        elif inTFright == args.chip:
            if float(overlaps[2]) >= float(5):
                overlapslist.append(overlaps)

#print(overlapslist)
'''Now loop over the list content (which is itself a list) and convert column 3 to float'''
for i in range(len(overlapslist)):
    overlapslist[i][2] = float(overlapslist[i][2])

'''Sort the list high to to low coassociation scores'''
overlapslist_sorted = []
for line in sorted(overlapslist, key=itemgetter(2), reverse=True):
    overlapslist_sorted.append(line)
#print(overlapslist_sorted)
'''Create interaction list'''
interacts_list = []
if len(list(args.chip.split('_'))) > 1:
    for j in overlapslist_sorted[1:int(args.top)]:
        if j[0].split('_')[0] not in interacts_list:
            interacts_list.append(j[0].split('_')[0])
        if j[1].split('_')[0] not in interacts_list:
            interacts_list.append(j[1].split('_')[0])
else:
    for j in overlapslist_sorted[1:int(args.top)]:
        if j[0].split('_')[0] not in interacts_list:
            interacts_list.append(j[0].split('_')[0])
        if j[1].split('_')[0] not in interacts_list:
            interacts_list.append(j[1].split('_')[0])
print("\n")
print("The identified interactors are: ", ','.join(interacts_list))
print("\n")
'''Find interactors in complex'''
#complex_count_list = [] #only for complement check by Counter fun
complex_dict = {}
complex_list = []
file = open(args.complex, "r")
complexesfile = file.read().strip().split('\n')

for interactors in interacts_list:
    for i in complexesfile:
        i = i.split("\t")
        complex = ''.join(i[1])
        components1 = i[18].split(';')
        components2 = re.split(', | |;', i[19])
        components = [*components1, *components2]
        #print(interactors)
        if interactors in components:
            #print(complex, interactors, components)
            #complex_count_list.append(complex)
            complex_list.append(''.join(complex + '\t' + interactors + '\t' + str(components)))
            list_components = list(components)
            while 'None' in list_components:
                list_components.remove('None')
            else:
                pass
            # print(list_components)
            complex_dict[complex] = list_components


#complex_count = Counter(complex_count_list)
# print(complex_count)
'''defaultdict will allow to append repeated TFs for same complex. The output of defaultdict is a dictionary'''
complex_tf_dict = defaultdict(list)
for j in complex_list:
    j = j.split('\t')
    complex_tf_dict[j[0]].append(j[1])
# #
# print(complex_dict)
# print(complex_tf_dict)
'''Use two dictionaries to get the required data'''
if os.path.exists(args.chip + '_' +'complex_mapper_out.txt'):
    os.remove(args.chip + '_' +'complex_mapper_out.txt')
else:
    pass
for items in complex_dict:
    print('_'.join(items.split(' ')) + '\t' + str(','.join(complex_tf_dict[items])) + '\t' +  str(','.join(complex_dict[items])) + '\t'+  str(len(complex_tf_dict[items])) + '\t'+  str(len(complex_dict[items]))  + '\t'+  str(len(complex_tf_dict[items])/len(complex_dict[items])))
    with open(args.chip + '_' + 'complex_mapper_out.txt', 'a') as mycomplex:
        mycomplex.write('_'.join(items.split(' ')) + '\t' + str(','.join(complex_tf_dict[items])) + '\t' +  str(','.join(complex_dict[items])) + '\t'+  str(len(complex_tf_dict[items])) + '\t'+  str(len(complex_dict[items]))  + '\t'+  str(len(complex_tf_dict[items])/len(complex_dict[items])) + '\n')
        mycomplex.close()


print('''
---------->>>> Now you can run for example below command to see complexes where the input protein and its interactors fall
''')

print("awk " + "'{print " + "$1" + '"\\t"' + "$2"   + '"\\t"'+ "$4" + '"\\t"' + "$5"  + '"\\t"'+ "$6" + "}' " + args.chip + "_complex_mapper_out.txt"  + " | sort -k3,3 -g -r | grep " + args.chip.split('_')[0] + " | head -n 10")


print('''
-------->>> or if your chip_protein is unknown, can check where interactors match with complex
''')

print("awk " + "'{print " + "$1" + '"\\t"' + "$2"  + '"\\t"'+ "$4" + '"\\t"' + "$5"  + '"\\t"'+ "$6" + "}' " + args.chip + "_complex_mapper_out.txt"  + " | sort -k3,3 -g -r | head -n 10 \n")