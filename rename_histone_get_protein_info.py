import os, re


def read_accession(ID):
    ids = ID.read().strip().split('\n')
    del ids[0]
    accession_dict = {}
    for links in ids:
        links = links.split('/')
        accession = links[6].split('.bam')[0]
        accession_dict[accession] = accession
    return accession_dict

def read_exp_sheet(sheet):
    experiment_sheet = sheet.read().strip().split('\n')
    del experiment_sheet[0]
    content_dict = {}
    for content in experiment_sheet:
        content = content.split('\t')
        newcontent = content[15].replace("/files/", "")
        content_dict[content[6]] = newcontent
    del content_dict["Target of assay"]
    return content_dict


def find_accession_in_sheet(acc, cont):
    accession_in_sheet_dict = {}
    count = 0
    # print(cont)
    for code in acc:
        # print(code)
        for values in cont:
            for matchpos in re.finditer(code, cont[values]):
                count += 1
                accession_in_sheet_dict[code] = ''.join(str(count) + '\t' + code + '\t' + values + "\t" + str(matchpos))
    return accession_in_sheet_dict

# if os.path.exists('new_script_to_add_id.sh'):
#     os.remove('new_script_to_add_id.sh')
def rename_and_adjust(code_and_sheet):
        print('Renaming the BAM files')
        for i in code_and_sheet:
            components = code_and_sheet[i].split("\t")[2]
            print(i, components)
            os.system('cp ' + i + '.bam ' + components + '_' + i + '.bam')

ID = open("files_filter_download_histone.txt", 'r')
sheet = open('experiment_report_2023_6_19_11h_28m.tsv', 'r')


myacc = read_accession(ID)
# print(myacc)
mysheet = read_exp_sheet(sheet)
# print(mysheet)
match_acc_sheet = find_accession_in_sheet(myacc, mysheet)
# print(match_acc_sheet)
# for i in match_acc_sheet:
#     print(i, match_acc_sheet[i])
# print(match_acc_sheet)
rename_and_adjust(match_acc_sheet)

