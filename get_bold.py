import argparse
parser= argparse.ArgumentParser(add_help=False)
parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help= "This script generates a line that can be executed to download coi-5p fasta sequences for a given species. Usage is: python bold_path.py -s 'genus_species' -o <output>") 
parser.add_argument("-s", help= "-s is the species name to be downloaded. In this fashion: 'genus_species'. It can be only 'genus'. Ex 'homo_sapiens'", required = "True")
parser.add_argument("-o", help= "output file where a line to download coi-5p will be written", required = "True")

args = parser.parse_args()
output=open(args.o, "w")

name = args.s.split("_")
genus=name[0]
if len(name) ==2:
    species=name[1]
    output.write("wget -q -O newseqs" + " " + "http://v3.boldsystems.org/index.php/API_Public/sequence?taxon=" + genus +"%20" + species + "&marker=COI-5P")
else:
    species = None
    output.write("wget -q -O newseqs" + " " + "http://v3.boldsystems.org/index.php/API_Public/sequence?taxon=" + genus +"%20" + "&marker=COI-5P")
