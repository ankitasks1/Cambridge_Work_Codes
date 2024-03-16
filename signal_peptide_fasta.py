import os
import re
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

def signal_over_fasta(in_sig, in_fa):
    signal_adjusted_fasta = []
    count = 0
    for signal in in_sig:
        for fas in in_fa:
            fas = fas.split('\t')
            #print(fas[3])
            for matchpos in re.finditer(signal, fas[3]):
                count += 1
                fasta_sig = ''.join(str(''.join(matchpos.group())) + '\t' + str(count) + '\t' + str(matchpos.start() + 1) + '\t' + str(matchpos.end()) + '\t' + '\t'.join(fas[0:4]))
                signal_adjusted_fasta.append(fasta_sig)
    return signal_adjusted_fasta


in_sig = ["[DE][A-Z][A-Z][A-Z]L[LI]", "Y[A-Z][A-Z][VILFWYM]", "NP[A-Z]Y", "GD[A-Z]Y"]
myfasta = read_fasta('/Users/ankitverma/Documents/peptide_data/Drosophila_melanogaster.BDGP6.32.pep.all.fa')
#print(myfasta)
# count = 0
# for i in myfasta:
#     i = i.split('\t')
#     count += 1
#     print(count, i[3])
out_signal_over_fasta = signal_over_fasta(in_sig, myfasta)
if os.path.exists("/Users/ankitverma/Documents/peptide_data/signal_over_fasta_output.txt"):
    os.remove("/Users/ankitverma/Documents/peptide_data/signal_over_fasta_output.txt")
count = 0
for i in out_signal_over_fasta:
    print(i)
    count += 1
    with open('/Users/ankitverma/Documents/peptide_data/signal_over_fasta_output.txt', 'a') as myfile:
        myfile.write(i + '\n')
        myfile.close()
print("\n\n\nAnalysis finished \n\n")

print(f'{count} hits were identified\n\nThe results are in signal_over_fasta_output.txt')

