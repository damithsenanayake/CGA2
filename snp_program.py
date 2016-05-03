import sys
import pysam as ps

'''
Assignment 2 : Author - Damith Senanayake
Task 2.

a. output summary for rs12913832:

sample 1 :
There are 36 reads overlapping the coordinate at
chr15:28365618
The reads support the following nucleotides (total/with minimum base quality of 20):
A: 22/20 (of these 36%/40% map to the forward strand)
C: 0/0
T: 0/0
G: 14/14 (of these 42%/42% map to the forward strand)
N: 0/0
C/T haplotype
sample 2 :
There are 33 reads overlapping the coordinate at
chr15:28365618
The reads support the following nucleotides (total/with minimum base quality of 20):
A: 33/33 (of these 48%/48% map to the forward strand)
C: 0/0
T: 0/0
G: 0/0
N: 0/0
T/T haplotype
sample 3 :
There are 13 reads overlapping the coordinate at
chr15:28365618
The reads support the following nucleotides (total/with minimum base quality of 20):
A: 0/0
C: 0/0
T: 1/1 (of these 0%/0% map to the forward strand)
G: 12/12 (of these 25%/25% map to the forward strand)
N: 0/0
C/C haplotype
------------------------------------------------------------------------------------------------------------------------

output summary for rs1129038

sample 1 :
chr15:28356859
There are 38 reads overlapping the coordinate at
chr15:28356859
The reads support the following nucleotides (total/with minimum base quality of 20):
A: 2/1 (of these 50%/0% map to the forward strand)
C: 14/13 (of these 42%/46% map to the forward strand)
T: 20/17 (of these 50%/52% map to the forward strand)
G: 2/2 (of these 0%/0% map to the forward strand)
N: 0/0
Roughly equal C and T values : - > in the HERC2 gene, this would have a A/G  haplotype

sample 2 :
There are 26 reads overlapping the coordinate at
chr15:28356859
The reads support the following nucleotides (total/with minimum base quality of 20):
A: 0/0
C: 24/22 (of these 75%/72% map to the forward strand)
T: 0/0
G: 2/1 (of these 50%/0% map to the forward strand)
N: 0/0
(G/G haplotype)

sample 3:
There are 9 reads overlapping the coordinate at
chr15:28356859
The reads support the following nucleotides (total/with minimum base quality of 20):
A: 0/0
C: 1/1 (of these 0%/0% map to the forward strand)
T: 7/7 (of these 28%/28% map to the forward strand)
G: 1/1 (of these 0%/0% map to the forward strand)
N: 0/0
(: A/A haplotype)

b : In each of the sample the following information can be inferred.

for the SNP rs1129038:
Sample 1 has A/G haplotype and is heterozygous.
Sample 2 has G/G haplotype and is homozygous
Sample 3 has A/A haplotype and is homozygous

for the SNP rs12913832:
sample 1 has C/T haplotype and is heterozygous
sample 2 has T/T haplotype and is homozygous
sample 3 has C/C haplotype and is homozygous


c : According to the table 3 in the paper, the SNP, the samples are as follows

 sample         rs1129038       rs12913832
 1              A/G             C/T
 2              G/G             T/T
 3              A/A             C/C

 Therefore the probable colors of each sample are (confidence: based on table 3 of the paper. )

 1. Brown (81.9%)
 2. Brown (100%)
 3. Blue (97.7%)

d : Each of these may have their own confidence due to multiple reasons such as :

- other SNPs that may have significant correlation with the phenotypes
- Differences in expressing the genotypes
- other genes which may contribute to the phenotypes

'''

files = sys.argv

if not len(files) == 4:
    print "Error, run as \n" \
          "snp_program.py <file.bam> <seq>:<position> <quality_score>"

samfile = ps.AlignmentFile(files[1], "rb")

chromosome = files[2].split(":")[0]

position = int(files[2].split(":")[1]) - 1

fullcount = 0
forward_count = 0

counts_qual_for = {"T": 0, "A": 0, "G": 0, "C": 0, "N": 0}
counts_tot = {"T": 0, "A": 0, "G": 0, "C": 0, "N": 0}
counts_tot_for = {"T": 0, "A": 0, "G": 0, "C": 0, "N": 0}
counts_qual = {"T": 0, "A": 0, "G": 0, "C": 0, "N": 0}

total_count = 0
for read in samfile.fetch(chromosome, position, position + 1):
    total_count += 1
    counts_tot[read.seq[position - read.pos]] += 1
    if read.query_qualities[position - read.pos] >= int(files[3]):
        counts_qual[read.seq[position - read.pos]] += 1
        if not read.is_reverse:
            counts_qual_for[read.seq[position - read.pos]] += 1
    if not read.is_reverse:
        counts_tot_for[read.seq[position - read.pos]] += 1

print "There are " + str(total_count) + " reads overlapping the coordinate at"
print files[2]

print "The reads support the following nucleotides (total/with minimum base quality of " + str(files[3]) + "):"

out_strings = []
for k in counts_tot.keys():
    count_str = "0/0"
    if counts_tot[k] > 0:
        count_str = str(counts_tot[k]) + "/" + str(counts_qual[k]) + " (of these " + str(
            counts_tot_for[k] * 100 / counts_tot[k]) + "%/" + (str(0) if counts_qual[k] == 0 else  str(
            counts_qual_for[k] * 100 / counts_qual[k])) + "% map to the forward strand)"
    out_strings.append(k + ": " + count_str)

for string in out_strings:
    print string


samfile.close()
