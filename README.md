# wrapper-dev
Wrappers for alignment tools.

Example usage:  
Download and install a wrapper:  
```
./wrapper_daligner.py install
```
  
Run alignment (reads in the data folder are simulated reads):  
```
/wrapper_daligner.py run data/reads/reads-pacbio.fa data/reference/escherichia_coli.fa pacbio data/alignment/
```
Output of DALIGNER is then given in the file:
```
data/alignment/escherichia_coli.fa.reads-pacbio-daligner.fasta.las
```
Extracted overlaps (in a text form) can be found in:
```
data/alignment/escherichia_coli.fa.reads-pacbio-daligner.fasta.las.txt
```
