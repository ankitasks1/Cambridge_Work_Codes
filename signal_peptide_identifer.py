import re
import os
def read_fasta(in_fa):
    out_fa = []
    fasta = open(in_fa)
    fasta = fasta.read()
    for x in fasta.strip().split('>'):
        content = x.split('\n')
        #Since after splitting by '>' first element will be blank so len(element) >1 is required
        if len(content) > 1:
            identifier = content[0]
            identifier = identifier.split(' ')
            header = identifier[0]
            geneid = identifier[3].split(':')[1]
            symbol = identifier[7].split(':')[1]
            sequence = str(''.join(content[1:]))
            sequence = sequence.upper()
            data = ''.join(header + '\t' + geneid + '\t' + symbol + '\t' + sequence)
            out_fa.append(data)
    return out_fa

def read_tm(in_tm):
    out_tm = []
    tmfile = open(in_tm)
    tmfile = tmfile.read().strip().split('\n')
    for tmcontent in tmfile:
        tmcontent = tmcontent.split('\t')
        #print(tmcontent)
        tmid = tmcontent[0]
        tmval = str('\t'.join(tmcontent[0:7]))
        tmdata = ''.join(tmid + '\t' +tmval)
        out_tm.append(tmdata)
    return out_tm

def read_exp(in_exp):
    out_exp = []
    expfile = open(in_exp)
    expfile = expfile.read().strip().split('\n')
    for expcontent in expfile:
        expcontent = expcontent.split('\t')
        expid = expcontent[0]
        expval = str('\t'.join(expcontent[0:11]))
        expdata = ''.join(expid + '\t' + expval)
        out_exp.append(expdata)
    return out_exp

def signal_over_fasta(in_sig, in_fa):
    signal_over_fasta_list = []
    count = 0
    for signal in in_sig:
        for fas in in_fa:
            fas = fas.split('\t')
            #print(fas[3])
            for matchpos in re.finditer(signal, fas[3]):
                count += 1
                fasta_sig = ''.join(str(count) + '\t' + str(''.join(matchpos.group())) + '\t' +  str(matchpos.start() + 1) + '\t' + str(matchpos.end()) + '\t' + '\t'.join(fas[0:4]))
                signal_over_fasta_list.append(fasta_sig)
    return signal_over_fasta_list

def signal_over_fasta_tm(signal_fasta_file, tm_file):
    count = 0
    signal_over_fasta_tm_list = []
    for sig_fa in signal_fasta_file:
        sig_fa = sig_fa.split("\t")
        #print(sig_fa[4])
        for tm_info in tm_file:
            tm_info = tm_info.split("\t")
            if sig_fa[4] == tm_info[0]:
                count += 1
                sig_fa_tm_info = ''.join(str(count) + '\t' + str('\t'.join(sig_fa)) + '\t' + str('\t'.join(tm_info)))
                #print(sig_fa_tm_info)
                signal_over_fasta_tm_list.append(sig_fa_tm_info)
    return signal_over_fasta_tm_list

def signal_over_fasta_tm_expr(signal_fasta_tm_file, exp_file):
    count = 0
    signal_over_fasta_tm_expr_list = []
    for sig_fa_tm in signal_fasta_tm_file:
        sig_fa_tm = sig_fa_tm.split("\t")
        #print(sig_fa_tm)
        for exp_info in exp_file:
            exp_info = exp_info.split("\t")
            #print(exp_info)
            if sig_fa_tm[6] == exp_info[0]:
                count += 1
                sig_fa_tm_exp_info = ''.join(str(count) + '\t' + str('\t'.join(sig_fa_tm)) + '\t' + str('\t'.join(exp_info)))
                #print(sig_fa_tm_exp_info)
                signal_over_fasta_tm_expr_list.append(sig_fa_tm_exp_info)
    return signal_over_fasta_tm_expr_list


in_sig = ["[DE][A-Z][A-Z][A-Z]L[LI]", "Y[A-Z][A-Z][VILFWYM]", "NP[A-Z]Y", "GD[A-Z]Y"]
myfasta = read_fasta('/Users/ankitverma/Documents/peptide_data/withpython/Drosophila_melanogaster.BDGP6.32.pep.all.fa')
mytm = read_tm('/Users/ankitverma/Documents/peptide_data/withpython/transmembrane_me_count.txt')
myexp = read_exp('/Users/ankitverma/Documents/peptide_data/withpython/sc_expression_data.txt')

# count = 0
# for i in myfasta:
#     i = i.split('\t')
#     count += 1
#     print(count, i[3])
print('..............................')
print('Matching signal and fasta file')
print('..............................\n')
out_signal_over_fasta = signal_over_fasta(in_sig, myfasta)
count = 0
for i in out_signal_over_fasta:
    count += 1
print(f'{count} hits were identified\n')
#print(out_signal_over_fasta)
print('..............................................................')
print('Intersecting signal matched fasta file with Transmembrane file')
print('..............................................................\n')
out_signal_over_fasta_tm = signal_over_fasta_tm(out_signal_over_fasta, mytm)
count = 0
for i in out_signal_over_fasta_tm:
    count += 1
print(f'{count} hits were identified\n')

print('........................................................................')
print('Intersecting signal and TM helix matched file with expression data')
print('........................................................................\n')
out_signal_over_fasta_tm_exp = signal_over_fasta_tm_expr(out_signal_over_fasta_tm, myexp)

colname = "".join("num_exp"+"\t"+ "num_tm"+"\t"+ "num_fasta"+"\t"+ "signal"+"\t"+ "signal_start"+"\t"+ "signal_end"+"\t"+ "Protein_Stable_ID"+"\t"+ "Gene_Stable_ID"+"\t"+ "Gene_Symbol"+"\t"+ "Sequence"+"\t"+ "FlyBase_Protein_ID"+"\t"+ "FlyBase_Protein_ID_repeat"+"\t"+"Gene_Symbol"+"\t"+ "Flybase_Gene_ID"+"\t"+ "TM"+"\t"+ "TM_start"+"\t"+ "TM_end"+"\t"+ "gene.id1"+"\t"+ "gene.id2"+"\t"+ "gene.id3"+"\t"+ "aEC1gene_exp"+"\t"+ "aEC2gene_exp"+"\t"+ "aEC3gene_exp"+"\t"+ "aEC4gene_exp"+"\t"+ "pEC1gene_exp" +"\t"+ "pEC2gene_exp"+"\t"+ "pEC3gene_exp"+"\t"+ "dECgene_exp"+"\t"+ "mECgene_exp")

if os.path.exists("/Users/ankitverma/Documents/peptide_data/withpython/signal_over_fasta_tm_exp_output.txt"):
    os.remove("/Users/ankitverma/Documents/peptide_data/withpython/signal_over_fasta_tm_exp_output.txt")

with open('/Users/ankitverma/Documents/peptide_data/withpython/signal_over_fasta_tm_exp_output.txt', 'a') as myfile:
        myfile.write(colname + '\n')
        myfile.close()
#print(out_signal_over_fasta_tm_exp)
#print(out_signal_over_fasta_tm)

count = 0
for i in out_signal_over_fasta_tm_exp:
    #print(i)
    count += 1
    with open('/Users/ankitverma/Documents/peptide_data/withpython/signal_over_fasta_tm_exp_output.txt', 'a') as myfile:
        myfile.write(i + '\n')
        myfile.close()

print(f'{count} hits were identified\nThe results are in signal_over_fasta_tm_exp_output.txt')
