# Repository for Blastocrithidia project
## What you can find here

* Useful utils in `/bin` folder
* Sources for annotator and helpers for utils in `/lib` folder

## Utilities
### Translated blast hits .csv -> non-translated .gff converter

It can translate `.csv` files with amino acid coordinates like this
```
qseqid,qlen,sseqid,slen,length,evalue,pident,bitscore,mismatch,gaps,qstart,qend,sstart,send
0,26,NODE_1096_length_1211_cov_913.715_1,403,22,1.94e-10,100.000,55.5,0,0,5,26,1,22
10000,23,NODE_328_length_19730_cov_141.815_2,6576,23,8.45e-09,100.000,51.6,0,0,1,23,5176,5198
```
to `.gff` file with nucleotide coordinates like this:
```
##gff-version 3
NODE_1096_length_1211_cov_913.715	blast	gene	1	66	.	+	1	ID=0_6739010979747404_1
NODE_328_length_19730_cov_141.815	blast	gene	15527	15595	.	+	2	ID=10000_8207093410319402_2
```

#### Conventions
To determine nucleotide indicies the script shold know about amioacid *nucleotide length* of the sequence and a *frame*. We agreed to put these values in seqid in this format: 

`SEQID_length_<nucleotide length>[_optional_string]_<frame>`.

So, the length can be everywhere in seqid after the word `length`. The frame is the last number of the seqid.

#### Usage
The script has two required parameters
* -t (--target) determines what you want to get in the `.gff` annotation: query or subject
* -m (--mode) determines how the script format the seqid:
  * genome - simply cuts frame from seqid (`NODE_1096_length_1211_cov_913.715_1` becomes `NODE_1096_length_1211_cov_913.715`)
  * transcriptome - cuts everything after the length (`NODE_1096_length_1211_cov_913.715_1` becomes `NODE_1096`)

If `-in` (or `--input_file`) option is not provided the script will open the window and you can select the `.csv` file manually.
If `-out` (or `--output_file`) option is not provided the script will save the output file in the input file folder.

Examples:
```
bin/translated_hits_to_gff.rb -in hits.csv -out result.gff -t query -m genome
bin/translated_hits_to_gff.rb -t query -m genome
```
You can see full instructions by executing script with `-h` (or `--help`) option.
