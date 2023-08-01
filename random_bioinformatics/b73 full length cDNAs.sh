# Get b73 full length cDNAs: http://www.maizecdna.org/download/
# Combine them
cat maize_flcdna_1.txt maize_flcdna_group2.txt maize_flcdna_group3.txt maize_flcdna_group4.txt maize_flcdna_group5.txt > all_b73_flcdna.fasta
sed "/^>/s/$/_b73/" all_b73_flcdna.fasta > all_b73_flcdna_named.fasta
