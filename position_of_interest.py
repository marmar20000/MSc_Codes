##### Prediction in  Genes with Above a Gene Length Threshold and above a First Intron Length threshold ####
import numpy as np
def genes_above_threshold(file_with_genes, gene_of_interest_positions, intron_exons_junctions):
    with open(gene_of_interest_positions, 'w+') as interest_positions:
        with open(intron_exons_junctions, 'w+') as intron_exon_junctions:
            with open(file_with_genes, 'r+') as genes:
                gene_list = genes.readlines()
                for gene in gene_list:
                    name = gene.strip('\n').split('\t')[0]
                    chr =  gene.strip('\n').split('\t')[1]
                    start = gene.strip('\n').split('\t')[2]
                    end = gene.strip('\n').split('\t')[3]
                    length =gene.strip('\n').split('\t')[4]
                    list_with_exons = eval(gene.strip('\n').split('\t')[5])
                    starts =  list(np.arange(int(start), int(end), 20))
                    ends =   list(np.arange((int(start)+1), (int(end)+1),20))
                    number_of_exons = len(list_with_exons)
                    list_with_introns =  [(list_with_exons[i][1], list_with_exons[i + 1][0]) for i in range(len(list_with_exons) - 1)]
                    exon_rank = 0 
                    intron_rank = 0 
                    #print(list_with_exons)
                    #print(list_with_introns)
                    # report if exon as 1 and if it is intron as a 0
                    for exon in list_with_exons:
                        exon_rank += 1 
                        # generate random positions in all exons and report their rank and that they are exons (1)
                        random_numbers_list = np.random.randint(exon[0], exon[1], size=(int((exon[1]-exon[0])/20)))
                        for i in random_numbers_list:
                            interest_positions.write(str(chr)+'\t'+str(i)+'\t'+str(i+1)+'\t'+str(exon_rank)+'\t'+str(1)+'\n')
                    for intron in list_with_introns:
                        intron_rank += 1
                        #if intron[1]>intron[0]:
                        random_numbers_list = np.random.randint(intron[0], intron[1], size=(int((intron[1]-intron[0])/20)))
                        intron_exon_junctions.write(str(chr)+'\t'+str(intron[1])+'\t'+str(intron[1]+1)+'\t'+str(intron_rank)+'\t'+str(1)+'\n')
                        intron_exon_junctions.write(str(chr)+'\t'+str(intron[0]-1)+'\t'+str(intron[0])+'\t'+str(intron_rank)+'\t'+str(1)+'\n')
                        for i in random_numbers_list:
                            interest_positions.write(str(chr)+'\t'+str(i)+'\t'+str(i+1)+'\t'+str(intron_rank)+'\t'+str(0)+'\n')
                        #else:
                        #    print(intron, intron_rank, name)
                    #for i in range(int(int(length)/20)):
                    #    interest_positions.write(str(chr)+'\t'+str(starts[i])+'\t'+(str(ends[i]))+'\n')


genes_above_threshold('output_file_of_genes_greater_5Kbp_first_intron.tsv', 'positions_of_interest_fix.tsv', 'intron_exon_juctions_fix.tsv')

