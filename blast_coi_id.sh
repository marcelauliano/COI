#! /bin/bash
#Bash script to COI if fasta sequences. Under development

export PATH=/software/ensembl/compara/ncbi-blast-2.2.30+/bin:${PATH}

software="/software/team311/mu2/coi" #location of all scripts to run the COI-id pipeline
database="/lustre/scratch116/vr/projects/vgp/user/mu2/coi-db/coi.updated.fa" #latest database for now
FASTA="" #path to the fasta file just arriving from sequencing that you want to COI id.
fasta_prefix=${FASTA#*2020/}

seqtk sample $FASTA 40000 > $fasta_prefix.40000 #I'm subsampling the fasta file to have 40000 sequences. This is arbritary, can be changed. 
blastn -query $fasta_prefix.40000 -db $database -evalue 1e-05 -out $fasta_prefix.40000.blastn
wait
perl $software/blastcov.pl $fasta_prefix.40000.blastn >  $fasta_prefix.40000.blastn.cov #parsing blast output. This is the same as blast output format 8, with two extra colums. They represent, respectively, the % of the query and subject sequence in the blast match.
wait
cat $fasta_prefix.40000.blastn.cov | awk '{ if ($3>=97 && $14>=97) print $0}' >  $fasta_prefix.40000.blastn.cov.97

echo "All done. Please check file <output>.97 to confirm the COI id for your fasta"
