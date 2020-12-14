
## required software
# python2
# mmseqs2 https://github.com/soedinglab/MMseqs2
## on biohpc, `export PATH=/programs/mmseqs-avx2-11/bin:$PATH`
export PATH=/programs/mmseqs-avx2-11/bin:$PATH
## on farm, using precompiled `wget https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz; tar xvfz mmseqs-linux-sse41.tar.gz; export PATH=$(pwd)/mmseqs/bin/:$PATH
#export PATH=/group/jrigrp10/mstitzer/multiple/mmseqs/bin/:$PATH

### for nam, be able to access genome info on command line if needed
if [ ! -z "$1" ]
then
        GENOME=$1
        SHORTID=$2
        GENOMEFASTA=$3
        OUTDIR="${GENOME}/sine" ## need to remake because different from CONFIG
	OLDDB=$4
	TEMPDIR=$5
fi

CPU=1

mkdir -p $OUTDIR
mkdir -p $TEMPDIR

# #### get the SINE-Finder program from Wenke et al. 2011
# if [ ! -f sine_finder.py ]; then
# 	wget http://www.plantcell.org/content/suppl/2011/08/29/tpc.111.088682.DC1/Supplemental_Data_Set_1-sine_finder.txt
# 	#### change name
# 	mv Supplemental_Data_Set_1-sine_finder.txt sine_finder.py
# 	chmod 755 sine_finder.py
# fi

# #### run sinefinder
# #### I haven't been able to get sine_finder to work with reverse sequences, as it seems to report TSDs wrong on the reverse strand.
# ####   so I'm only reporting on the forward strand.
# ### -f both : outputs csv and fasta
# #python2 sine_finder.py -T chunkwise -V1 -f both -o F ${GENOMEFASTA}
# python2 sine_finder.py -T chunkwise -V1 -f both ${GENOMEFASTA}

# ## sine finder has gross hard coded output that throws output in same directory as the genome fasta, so move them here.
# SF_OUT="$(dirname -- $GENOMEFASTA)"
# mv ${SF_OUT}/*-matches.fasta $OUTDIR/$GENOME-matches.fasta
# mv ${SF_OUT}/*-matches.csv $OUTDIR/$GENOME-matches.csv

# ## output includes TSD, so remove
# python2 remove_tsd_sinefinder.py $OUTDIR/${GENOME}-matches.fasta $OUTDIR/${GENOME}-matches.noTSD.fa


## cluster families - mmseqs version, still single linkage clustering
#wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz; tar xvfz mmseqs-linux-avx2.tar.gz;export PATH=$(pwd)/mmseqs/bin/:$PATH

#### these are what I did to generate the first db, using old results from v4, Mo17, W22, and Ph207
#mmseqs createdb post-Mo17.existingRST.fa maize.SINE ## this generates the maize.SINE db for the first time (and these IDs are flipped to old fam names)
#mmseqs cluster maize.SINE maize.SINE.8080 /scratch/  --min-seq-id 0.8 -c 0.8 --cov-mode 1 --cluster-mode 2

## now create a db for the new genome (in a directory because there's lots of output files)
mkdir -p mmseqsDBs
mmseqs createdb $OUTDIR/${GENOME}-matches.noTSD.fa mmseqsDBs/${GENOME}.SINE 

## concatenate the old (from cmd line) to the new genome entries
mmseqs concatdbs mmseqsDBs/$OLDDB mmseqsDBs/${GENOME}.SINE  mmseqsDBs/maize.SINE.${GENOME}.combined
mmseqs concatdbs mmseqsDBs/${OLDDB}_h mmseqsDBs/${GENOME}.SINE_h mmseqsDBs/maize.SINE.${GENOME}.combined_h
#### then need to cluster using the updated clusters
## kinda confusing list of arguments, shown here for clarity
# mmseqs clusterupdate <oldDB> <newDB> <oldDB_clustering> <outDB> <tmpDir> [opts] ## was getting confused with oldDB_clustering and rewriting names!
mmseqs clusterupdate mmseqsDBs/$OLDDB mmseqsDBs/maize.SINE.${GENOME}.combined mmseqsDBs/$OLDDB.8080 mmseqsDBs/maize.SINE.post-${GENOME} mmseqsDBs/maize.SINE.post-${GENOME}.8080 $TEMPDIR --min-seq-id 0.8 -c 0.8 --cov-mode 1 --cluster-mode 2 --search-type 3
mmseqs createtsv mmseqsDBs/maize.SINE.post-${GENOME} mmseqsDBs/maize.SINE.post-${GENOME} mmseqsDBs/maize.SINE.post-${GENOME}.8080 maize.SINE.post-${GENOME}.8080.tsv
 
#mmseqs clusterupdate maize.SINE ${GENOME}.SINE maize.SINE.8080 maize.SINE.post-${GENOME} maize.SINE.post-${GENOME}.8080 /scratch/ --min-seq-id 0.8 -c 0.8 --cov-mode 1 --cluster-mode 2

#mmseqs createtsv maize.SINE maize.SINE maize.SINE.8080 maize.SINE.8080.tsv
 
#mmseqs createdb DB.fasta DB_new
#mmseqs clusterupdate DB_trimmed DB_new DB_trimmed_clu DB_new_updated DB_update_clu tmp

## then generate a gff3
Rscript generate_gff_SINE.R $GENOME $SHORTID $OUTDIR

