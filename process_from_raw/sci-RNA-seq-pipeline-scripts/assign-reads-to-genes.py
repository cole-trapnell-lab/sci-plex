#!/usr/bin/env python

import sys

empty_gene_set = set(["."])

gene_start = dict()
gene_end = dict()

with open(sys.argv[1]) as genes_bed:
    for line in genes_bed:
        fields = line.strip().split("\t")
        gene_start[fields[3]] = int(fields[1])
        gene_end[fields[3]] = int(fields[2])

for line in sys.stdin:
    fields = line.strip().split("\t")
    
    read_info = fields[0]
    read_chr = fields[1]
    read_start = int(fields[2])
    read_end = int(fields[3])
    read_strand = fields[5]
    
    alignment_type = "intergenic";

    if fields[6] == ".":
        # no read fragment is contained within an exon on the proper strand
        if fields[7] == ".":
            # no read fragment is contained within a gene on the proper strand
            # gapped reads will not follow this code path because they will have fields[7]
            # be like ".,."
            # we'll handle these later
            sys.stdout.write("%s\t%s\t%s\n" % (read_info, "NA", alignment_type))
            continue
        else:
            alignment_type = "intronic"
            compatible_gene_sets = [set(x.split("|")) - empty_gene_set for x in fields[7].split(",")]
    else:
        alignment_type = "exonic";
        compatible_gene_sets = [set(x.split("|")) - empty_gene_set for x in fields[6].split(",")]

    compatible_genes_intersection = compatible_gene_sets[0]
    compatible_genes_union = compatible_gene_sets[0]

    for gene_set in compatible_gene_sets[1:]:
        compatible_genes_intersection &= gene_set
        compatible_genes_union |= gene_set

    if len(compatible_genes_union) == 0:
        # gapped read is not contained within any gene on the proper strand, so skip it
        # (ungapped reads that do not overlap any gene already have been handled above)
        alignment_type = "intergenic"
        sys.stdout.write("%s\tNA\t%s\n" % (read_info, alignment_type))
        continue
    elif len(compatible_genes_union) == 1:
        # good: read is consistent with sequence from exactly one gene
        gene = next(iter(compatible_genes_union))
        sys.stdout.write("%s\t%s\t%s\n" % (read_info, gene, alignment_type))
        continue
    else:
        if len(compatible_genes_intersection) == 0:
            # bad: read has a gapped alignment and fragments of the alignment
            #      map to different non-overlapping genes on the same strand,
            #      shouldn't be possible

            sys.stdout.write("%s\tNA\tinconsistent_%s\n" % (read_info, alignment_type))
            continue
        elif len(compatible_genes_intersection) == 1:
            # good: fragments of gapped read alignment are compatible with multiple genes
            #       but only one gene is consistent with all alignment fragments
            gene = next(iter(compatible_genes_intersection))
            sys.stdout.write("%s\t%s\t%s\n" % (read_info, gene, alignment_type))
        else:
            # read is ambiguous, being compatible with more than one gene
            # e.g. contained within intronic sequence that's shared by more than one gene,
            # or contained within an exon that's the last exon for one gene,
            # but also used as an internal exon by another gene (ugh)

            nearby_genes = []
            for gene in compatible_genes_intersection:
                if read_strand == "+":
                    dist = gene_end[gene] - read_end
                elif read_strand == "-": 
                    dist = read_start - gene_start[gene]
                if dist >= 0 and dist <= 100:
                    nearby_genes.append(gene)

            if len(nearby_genes) == 1:
                sys.stdout.write("%s\t%s\t%s\n" % (read_info, nearby_genes[0], alignment_type))
                continue          
            else:
                sys.stdout.write("%s\tNA\tambiguous_%s\n" % (read_info, alignment_type))
                continue;

