grep $'\tgene\t' genome.gff | awk -F '\t' -v OFS='\t' '{print $1, $4 -1, $5, '.', '.', $7}' > bedfile.bed