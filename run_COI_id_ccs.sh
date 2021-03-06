#!bin/bash
software="/software/team311/mu2/coi"
database="/lustre/scratch116/vr/projects/vgp/user/mu2/coi-db/coi.updated.fa" #latest database for now
fasta="$1"
fasta_prefix=${fasta#*2020/}
database_prefix=${database#*coi-db/}

#Step 1: replace variable new_species with the scientific name of the species you want to include COI-5P sequences for in the blast database.
new_species="$2" #Replace this with species name underscored. Ex: homo_sapiens.


if [ "$1" == "-h" ]; then
        echo "Usage: <path/species.ccs.fasta>  <genus_species> "
        echo -e "<path/species.ccs.fasta>  fasta file to be analyzed. Usually a gzip fasta file CCS PacBio"     
        echo -e "<genus_species> New species to be added to the blast database separated by a underscore. Ex: homo_sapiens"     
        exit 0
fi


if [ -z $1 ]; then
        echo "fasta sequence not provided"
        exit 1
else
        echo "fasta file is $fasta"
fi


if [ -z $2 ]; then
        echo "No new species to include in the blast db? This must be wrong"
        exit 1
else
        echo "New species to include on the blast db is  $new_species"
fi

echo "First let's subsample our fasta file to contain 30000 sequences"

seqtk sample $fasta 20000 > $fasta_prefix.20000

echo "Seqk done. Subsampled fasta file is called $fasta_prefix.20000"

echo "Downloading and adding COI-5P sequences for the new species to the blast db"

python $software/get_bold.py -s $new_species -o out #python script that generates a wget line for coi-5p of the species present in the variable new_specie 
sh out
sleep 2

if [[ -f newseqs ]]
then
        echo "File 'newseqs' created. Editing and adding it to blast db"
else
        echo "no newseqs created"
        exit 1
fi
cat newseqs $database  > $database_prefix.newseqs
sed s'/ /_/g' $database_prefix.newseqs > $database_prefix.new_seqs
sleep 2
echo "Adding done. Making blastdb and running blast"

#Step 1.1: let's create a new blastdb including COI-5P sequences from out new species
export PATH=/software/ensembl/compara/ncbi-blast-2.2.30+/bin:${PATH}
makeblastdb -in $database_prefix.new_seqs -dbtype nucl

blastn -query $fasta_prefix.20000 -db $database_prefix.new_seqs -evalue 1e-05 -out $fasta_prefix.20000.blastn
sleep 2
perl $software/blastcov.pl $fasta_prefix.20000.blastn >  $fasta_prefix.20000.blastn.cov #parsing blast output. This is the same as blast output format 8, with two extra colums. They represent, respectively, the % of the query and subject sequence in the blast match.
sleep 2
cat $fasta_prefix.20000.blastn.cov | awk '{ if ($3>=97 && $14>=97) print $0}' >  $fasta_prefix.20000.blastn.cov.97

rm out newseqs $fasta_prefix.20000.blastn $fasta_prefix.20000.blastn.cov $database_prefix.new_seqs* $fasta_prefix.20000 $database_prefix.newseqs
echo "Intermediate files removed"
echo "All done. Please check file <output>.97 to confirm the COI id for your fasta"
