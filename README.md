## wrapper-dev - Wrappers for alignment tools.  
Wouldn't it be great if you could simply run any aligner/mapper with the same command line, and get a SAM file as an output, while also measuring time/memory usage statistics?  
And setting up the tool would simply be done by saying 'install'?  

Well, that's what this project is all about.  

This is an early development version, and currently only DALIGNER is implemented.  
Also note that the command line interfaces might change from version to version, until it iterates out to a stable form.  

Stay tuned, more comming soon!  


***Example usage:***  
Download and install a wrapper:  
```
./wrapper_daligner.py install  
```  
  
Run alignment (reads in the data folder are simulated reads):  
```
/wrapper_daligner.py run data/reads/reads-pacbio.fa data/reference/escherichia_coli.fa pacbio data/alignment/ suffix  
```  
SAM file is then given in the file:  
```
data/alignment/DALIGNER-suffix.sam  
```  

To specify the full name of the output SAM file, simply add a '.sam' to the last parameter.  
```
/wrapper_daligner.py run data/reads/reads-pacbio.fa data/reference/escherichia_coli.fa pacbio data/alignment/ output.sam  
```  
SAM file is then given in the file:  
```
data/alignment/output.sam  
```  

To run DALIGNER as an overlapper:  
```
/wrapper_daligner.py overlap data/reads/reads-pacbio.fa data/reference/escherichia_coli.fa pacbio data/alignment/  
```  
Output of DALIGNER is then given in the file:  
```
data/alignment/escherichia_coli.fa.reads-pacbio-daligner.fasta.las  
```  
Extracted overlaps (in a text form) can be found in:  
```
data/alignment/escherichia_coli.fa.reads-pacbio-daligner.fasta.las.txt  
```  
