with open('genesR.txt') as f:
    content = f.readlines()

wanted = [line.rstrip('\n') for line in open('genesR.txt')]


from Bio import SeqIO
genbank_file = "Gmerged.gbff" 
for record in SeqIO.parse("Gmerged.gbff", "genbank"):
    for f in record.features:
    if f.type == "mRNA" and "locus_tag" in f.qualifiers:
        locus_tag = f.qualifiers["locus_tag"][0]
        if locus_tag in wanted:
        try:
            out1 = f.qualifiers["locus_tag"][0] + "," + f.qualifiers["gene"][0] + "\n"
            file.write(out1)
        except:
            out2 = f.qualifiers["locus_tag"][0] + ", no gene name" + "\n"
            file.write(out2)
