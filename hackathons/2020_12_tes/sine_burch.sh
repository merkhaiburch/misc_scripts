# Script to create SINE annotaitons for EP1


# get code from Michelle's directory on blfs1
scp -r mbb262@cbsublfs1.tc.cornell.edu:/data1/users/mcs368/annotaTE_TEs/sine /workdir/mbb262

# Get EP1 genome from maizegdb
wget https://download.maizegdb.org/Zm-EP1-REFERENCE-TUM-1.0/Zm-EP1-REFERENCE-TUM-1.0.fa.gz /workdir/mbb262/sine
gunzip /workdir/mbb262/sine/Zm-EP1-REFERENCE-TUM-1.0.fa.gz

# Need to change permissions of sine directory
chmod 777

# Run script
./annotate_sine.sh EP1 Zm00010a /workdir/mbb262/sine/Zm-EP1-REFERENCE-TUM-1.0.fa /workdir/mbb262/sine/maize.SINE.post-B73Ab10

# Run second version of script
./annotate_sine_2.sh \
    EP1 Zm00010a \
    /workdir/mbb262/sine/Zm-EP1-REFERENCE-TUM-1.0.fa \
    /workdir/mbb262/sine/maize.SINE.post-B73Ab10 \
    /workdir/mbb262/temp


# Run code while supplying tempdir
# $OUTDIR = /workdir/mbb262/sine
## now create a db for the new genome (in a directory because there's lots of output files)
mkdir -p mmseqsDBs
mmseqs createdb /workdir/mbb262/sine/EP1/sine/EP1-matches.noTSD.fa /workdir/mbb262/sine/mmseqsDBs/EP1.SINE 

## concatenate the old (from cmd line) to the new genome entries
OLDDB=/workdir/mbb262/sine/mmseqsDBs/maize.SINE.post-B73Ab10
mmseqs concatdbs $OLDDB \
    /workdir/mbb262/sine/mmseqsDBs/EP1.SINE \
    /workdir/mbb262/sine/mmseqsDBs/maize.SINE.EP1.combined

mmseqs concatdbs \
    ${OLDDB}_h \
    /workdir/mbb262/sine/mmseqsDBs/EP1.SINE_h \
    /workdir/mbb262/sine/mmseqsDBs/maize.SINE.EP1.combined_h

#### then need to cluster using the updated clusters
## kinda confusing list of arguments, shown here for clarity
# mmseqs clusterupdate <oldDB> <newDB> <oldDB_clustering> <outDB> <tmpDir> [opts] ## was getting confused with oldDB_clustering and rewriting names!

# New code
OLDDB=/workdir/mbb262/sine/mmseqsDBs/maize.SINE.post-B73Ab10
mmseqs clusterupdate \
    /workdir/mbb262/sine/mmseqsDBs/maize.SINE.post-B73Ab10 \
    /workdir/mbb262/sine/mmseqsDBs/maize.SINE.EP1.combined \
    /workdir/mbb262/sine/mmseqsDBs/maize.SINE.post-B73Ab10.8080 \
    /workdir/mbb262/sine/mmseqsDBs/maize.SINE.post-EP1 \
    /workdir/mbb262/sine/mmseqsDBs/maize.SINE.post-EP1.8080 \
    /workdir/mbb262/temp \
    --min-seq-id 0.8 -c 0.8 --cov-mode 1 --cluster-mode 2 --search-type 3

mmseqs createtsv \
    /workdir/mbb262/sine/mmseqsDBs/maize.SINE.post-EP1 \
    /workdir/mbb262/sine/mmseqsDBs/maize.SINE.post-EP1 \
    /workdir/mbb262/sine/mmseqsDBs/maize.SINE.post-EP1.8080 \
    maize.SINE.post-EP1.8080.tsv
 

## then generate a gff3
Rscript generate_gff_SINE.R $GENOME $SHORTID $OUTDIR
