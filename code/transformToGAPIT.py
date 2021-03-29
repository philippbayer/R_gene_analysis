import gzip
with gzip.open('./data/NBS_PAV.txt.gz', mode='rt') as fh:
    header = fh.readline().split()

    # ['Individual', 'AB-01', 'AB-02', 'BR-01', 'BR-02', 'BR-03', 'BR-04', #
    ind_dict = {ind:list() for ind in header[1:]}

    genes = []
    for line in fh:
        ll = line.split()
        genes.append(ll[0])
        for ind, allele in zip(header[1:], ll[1:]):
            if allele == '1':
                # 1 is heterozygous in GAPIT
                allele = '2'
            ind_dict[ind].append(allele)

# We also need the genes' positions, for now I'll just take the middle of the gene model
gene_set = set(genes)
# GlymaLee.18G225400.1.p
pos_dict = {}

with open('data/Lee.pan.v1.renamed.gff') as gff:
    for line in gff:
        ll = line.split()
        # ['Gm01', 'phytozomev13', 'gene', '37775', '37993', '.', '+', '.', 'ID=GlymaLee.01G000100.1.v1.1;Name=GlymaLee.01G000100']
        if ll[2] != 'gene':
            continue
        names = ll[-1].split(';')
        thisid = names[0]
        assert 'ID' in thisid
        thisid = thisid.replace('ID=','').replace('.v1.1', '.p')
        if thisid not in gene_set:
            # this an NLR?
            continue
        chrom, start, end = ll[0], ll[3], ll[4]
        start, end = int(start), int(end)
        if start > end:
            start, end = end, start
        pos = start + (end-start)/2
        pos = int(pos) # python2 has an int here, python3 has a float here
        pos = str(pos)
        # GAPIT expects numeric chromosomes, which we don't have. let's transform
        if 'Gm' in chrom:
            # Gm01 becomes 1
            chrom = str(int(chrom.replace('Gm',''))) # the str(int) thing is a hack to make 0 go away from Gm01
        elif chrom.startswith('sc'):
            # scaffolds become 21
            chrom = '21'
        elif chrom.startswith('Uwa'):
            # pangenome becomes 22
            chrom = '22'
        else:
            assert False, chrom
        pos_dict[thisid] = [chrom, pos]

# now let's write the numeric table GAPIT wants

newheader = ['taxa'] + genes
with open('data/NLR_PAV_GD.txt', 'w') as geno_out, open('data/NLR_PAV_GM.txt', 'w') as map_out:
    geno_out.write('\t'.join(newheader) + '\n')
    map_out.write('Name\tChromosome\tPosition\n')
    for i in header[1:]:
        this_line = [i] + ind_dict[i]
        geno_out.write('\t'.join(this_line) + '\n')
    for g in genes:
        this_m = [g] + pos_dict[g]
        map_out.write('\t'.join(this_m) + '\n')

