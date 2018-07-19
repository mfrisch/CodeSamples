#========================================================
# Contributors: Michael Frisch, Kranthi Vavikolanu
# Read Lab @ Emory University
# Based on Mothur MiSeq SOP: http://www.mothur.org/wiki/MiSeq_SOP
#========================================================

#========================================================
# Variable assignment
#========================================================
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )" # Find shell file's directory
DATE_MADE=$(date +"%Y-%h-%d") # Set date for labeling biom files
PROC=8 #------------ Number of processors to use when running analysis
MINLEN=1 #---------- The minimum size of contigs
MAXLEN=385 #-------- The maximum size of contigs
SILVASTART=13862 #-- Beginning of 16s region of interest
SILVAEND=23444 #---- End of 16s region of interest
SCREENSTART=1046 #-- After the alignment step, position where aligned sequences should start
SCREENEND=6424 #---- After the alignment step, position where aligned sequences should end
CLASSCUTOFF=80 #---- Minimum sequence identity between during classify.seqs
TAXLEVEL=6 #-------- Depth of taxonomic rank desired (max = 6 if using SILVA file): Kingdom = 1, Phylum = 2, ..., Species = 7
TAXCUTOFF=0.05 #---- Cutoff value for clustering 
METADATA=false #---- Boolean value to indicate if metadata should be included in biom file

# Check if 'stability.files' exists in directory
if [ ! -e $DIR/stability.files ]
then
	echo Missing stability.files in the directory $DIR
	exit 1
fi
# Check if 'metadata.txt' exists in directory
if [ $METADATA ] && [ ! -e $DIR/metadata.txt ]
then
	echo Missing metadata in the directory $DIR
	exit 1
fi

#========================================================
# Quality Control for biom Files
#========================================================
mothur "#make.contigs(file=stability.files, processors=$PROC);\
	summary.seqs(fasta=stability.trim.contigs.fasta, processors=$PROC);\
	screen.seqs(fasta=stability.trim.contigs.fasta, group=stability.contigs.groups, maxambig=0, maxlength=$MAXLEN);\
	summary.seqs(fasta=stability.trim.contigs.good.fasta);\
	unique.seqs(fasta=stability.trim.contigs.good.fasta);\
	count.seqs(name=stability.trim.contigs.good.names, group=stability.contigs.good.groups);\
	pcr.seqs(fasta=$DIR/silva.nr_v119.align, start=$SILVASTART, end=$SILVAEND, keepdots=F, processors=$PROC);\
	system(mv $DIR/silva.nr_v119.pcr.align silva.v4.fasta);\
	summary.seqs(fasta=stability.trim.contigs.fasta, processors=$PROC);\
	align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.nr_v119.align);\
	summary.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table);\
	screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, summary=stability.trim.contigs.good.unique.summary, start=$SCREENSTART, end=$SCREENEND, maxhomop=8);\
	filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=.);\
	unique.seqs(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table);\
	pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.good.unique.good.filter.count_table, diffs=2);\
	chimera.uchime(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t);\
	remove.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.accnos);\
	classify.seqs(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table, reference=$DIR/silva.nr_v119.align, taxonomy=$DIR/silva.nr_v119.tax, cutoff=$CLASSCUTOFF);\
	remove.lineage(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v119.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)"

#========================================================
# Analysis for Silva biom File
#========================================================
mothur "#cluster.split(fasta=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v119.wang.pick.taxonomy, splitmethod=classify, taxlevel=$TAXLEVEL, cutoff=$TAXCUTOFF, processors=$PROC);\
	make.shared(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table);\
	classify.otu(list=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v119.wang.pick.taxonomy);\
	make.biom(shared=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared, constaxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.unique.cons.taxonomy, metadata=metadata.txt)"

mv stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.unique.biom stability.$DATE_MADE.biom

#========================================================
# Production of GreenGenes biom File for PICRUst Analysis
#========================================================
mothur "#system(cp stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta stability_gg.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta);\
	cluster.split(fasta=stability_gg.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v119.wang.pick.taxonomy, splitmethod=classify, taxlevel=$TAXLEVEL, cutoff=$TAXCUTOFF, processors=$PROC);\
	make.shared(list=stability_gg.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table);\
	classify.otu(list=stability_gg.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=stability.trim.contigs.good.unique.good.filter.unique.precluster.uchime.pick.pick.count_table, taxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v119.wang.pick.taxonomy);\
	make.biom(shared=stability_gg.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared,  reftaxonomy=$DIR/gg_13_5_99.gg.tax,  constaxonomy=stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.unique.cons.taxonomy, metadata=metadata.txt, picrust=$DIR/99_otu_map.txt)"
mv stability_gg.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.unique.biom stability_gg.$DATE_MADE.biom

#========================================================
# Analysis Using PICRUst
#========================================================
#normalize_by_copy_number.py -i stability_gg.$DATE_MADE.biom -o normalized_otus.biom
#predict_metagenomes.py -i normalized_otus.biom -g 13_5 -o metagenome_predictions.biom --with_confidence
#categorize_by_function.py -i metagenome_predictions.biom -o predicted_metagenomes.L3.biom -c KEGG_Pathways -l 1
#metagenome_contributions.py -i normalized_otus.biom -o ko_metagenome_contributions.tab

#========================================================
# Linear Discriminate Analysis Effect Size (LEfSe)
#========================================================
#mothur "#summary.seqs(stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta,processors=$PROC); make.shared(biom=stability_gg.$DATE_MADE.biom);  lefse(shared=stability_gg.$DATE_MADE.shared, design=sampledata, class=Extraction_Method, subclass=BodySite)"

#========================================================
# Clean Up
#========================================================
#rm *.filter stability.*.align *.sum *.train *.shared *.logfile *.groups *.summary *.count_table *.taxonomy *.rabund *.list *.map *.fasta *.accnos *.names *.qual *.chimeras *.8mer *.prob *.numNonZero *.report *.dist *.temp

#========================================================
# Make Report
#========================================================
#Rscript Phyloseq_mfrisch.R $DIR

