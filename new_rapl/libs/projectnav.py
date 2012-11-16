content = """<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" 
  "http://www.w3.org/TR/html4/loose.dtd">
<html><head>
<meta http-equiv="content-type" content="text/html; charset=ISO-8859-1">
<title>RAPL output</title>
</head>
<body>

generated: 2012-06-25 05:50

<h1>RAPL output</h1>

<h2><a name="toc">Table of content</a></h2>

<ul>
  <li><a href="#introduction">Introduction</a></li>
  <li><a href="#input_files">Input files</a></li>
  <li><a href="#read_mappings">Read mappings</a></li>
  <li><a href="#coverages">Read coverages</a></li>
  <li><a href="#coverages">Coverager files</a></li>
  <li><a href="#read_quantification">Feature based read quantification</a></li>
  <li><a href="#report">Reports and statistics</a></li>
 </ul>

<h2><a name="introduction">Introduction</a></h2>
<p><a href="#toc">top</a></p>

<p> PERMANENET URL RAPL is a pipeline to perform the common steps of
an RNA-Seq analysis. This overview intends to guide you trough the
input and ouput</p>

<p>At first: Relax - advanced skills in bioinformatics/computation
biology are not required. Most of the files are plain text files which
can be opened with simple text viewer/ editor programs like less,
gedit, emacs, vi in GNU/Linux, Mac OS and other Unix-like operating
systems. In Windows notepad is a sufficient tool but the files can
also be loaded in more sophisticated office applications like
OpenOffice Doc or Microsoft Word. Some files are so called CSV files
(character separeted files) in which columns are separated by
tabs. These files can also be loaded in spreadsheet application like
LibreOffice Calc or Microsoft Excel and then saved in the native
format of the application. The graphics are stored as PDF files and
are generated using the programming language R. The names of the R
scripts end with .R and can be modified with the above mentioned text
editors.</p>

<p>To make the navigation to required files and the generated output
easy and intuitive there are two folder in the root of the project
folder. The folder $input_folder contains the files required
for the run. In the subfolder $RNA_lib_folder the read files
are located, the subfolder $genome_file_folder contains the
genome files (i.e. the FASTA files of all chromosom and plasmids)
while the subfolder $annotation_file_folder harbours
annotation files.

The created files are stored in subfolder of the folder
$output_folder. Detailed description of theses subfolder and
their content can be found in the following sections.</p>

<h2><a name="input_files">Input files</a></h2>
<p><a href="#toc">top</a></p>

<p>
<a href="./RAPL_analysis/input">input</a>
<a href="./RAPL_analysis/input/annotations">annotations</a>
<a href="./RAPL_analysis/input/genomes">genome</a>
<a href="./RAPL_analysis/input/reads">reads</a>
</p>

<p>The first step of the process is the actual sequencing of the cDNA
libraries. The result of this sequencing is a set of FASTA files which
contain the raw reads. The barcode sequenced that are used for
multiplexed runs (more than one lib per sequencing lane) are already
removed from these. For your project we used the following read files
stored in the folder $RNA_lib_folder:</p>

<h2><a name="read_mappings">Read mappings</a></h2>
<p><a href="#toc">top</a></p>

<p>The mapping against the reference genome was done using the program
segemehl (Hoffmann et. al, PLoS Computational Biology, 2009) and took
place in two steps. In the first step it was tried to align all raw
reads to the reference genome. Usually this is not successful for all
of them e.g due to sequencing errors. Some read cannot be mapped as
they contain a poly-A tail. To tackle this problem all reads for which
the mapping failed in the first run were undergoing a procedure in
which a potential poly-A tail was searched and clipped. After the
clipping the resulting sequences were filtered by sequence size: All
sequences that were shorter than $min_seq_length nucleotides were
removed. The other remaining reads were aligned to the genome in a
second round. After this the mapping results of the two runs were
combined. These reads of these combined files were then filtered by
their A-content. All entries with more than $max_a_content\\% of A's
were removed.</p>

<p>If you look at the mapping results you might noticed that some reads
are listed more than once. We allowed only one best segemehl
hit per read but for certain reads there were some hits sharing the
same best value. In such case two or more placements are accepted.</p>

<a href="./RAPL_analysis/output/read_mappings-clipped_reads">read_mappings-clipped_reads</a>
<a href="./RAPL_analysis/output/read_mappings-index">read_mappings-index</a>
<a href="./RAPL_analysis/output/read_mappings-mappings">read_mappings-mappings</a>
<a href="./RAPL_analysis/output/read_mappings-read_tracing">read_mappings-read_tracing</a>
<a href="./RAPL_analysis/output/read_mappings-unmapped_reads">read_mappings-unmapped_reads</a>

<h2><a name="coverages">Read coverages</a></h2>
<p><a href="#toc">top</a></p>

<a href="./RAPL_analysis/output/gr-coverages_raw">wiggle</a>
<a href="./RAPL_analysis/output/gr-coverages_read_normalized">wiggle norm</a>


<p>Based on the combined and A-content filtered mapping results GR file
that can be viewed in the IGB were generated. These GR files describe
the amount of mapped reads per nucleotide position for a each genome
file. Above it was mentioned that some reads are equally well mapped
to different locations on a genome. In such a case we divide the raw
counting (i.e. 1) by the number of mappings. E.g. if a read is mapped
to 4 different places each position only get 1/4 point counted for
this read.</p>

<p>For each genome file (chromosome/plasmid) two different kinds of sets
of GR files were generated. One set represents the raw countings of
reads, the other one contains normalized values to make a comparative
view on the data possible. The normalization was done by dividing the
value of each position by the number of successfully mapped
reads. Unfortunately the IGB has problems with small number and will
round them down to zero. Due to this the normalized values were
multiplied with a constant value for each set. We decided to take the
smallest number of successfully mapped reads in this set. As an
example: We have the two libraries Foo+ and Foo- mapped to the genomes
file Bar007.fa. If Bar007.fa has 2000,000 mappings of reads from Foo+
and 1000,000 from Foo- the values for Foo+ are divided by 2000,000 and
then multiplied by 1000,000. The values for lib Foo- are divided by
1000,000 and then multiply by 1000,000 (yes, for the libraries with
the lowest number of mapped reads the whole normalization and
multiplication process results in the raw countings).</p>

<p>For each normalized/not-normalized library-genome-file-pair you have
two different GR files. One for the plus strand, one for the minus
strand. You can recognize the type by the file name. The plus strand
files end with "plus.gr", the minus strand files with
"minus.gr". While the plus strand file should contain only positive
intensities, the minus strand file should contain only negative
intensities.</p>

<p>The Integrated Genome Browser is an open-source application that can be
downloaded from http://www.bioviz.org/igb/. It requires Java.</p>

<p>To make sense of the coordinates of the GR files data of the reference
genome has to be supplied. This data is taken from a FASTA file (.fna)
and a gff file of the analysed chromosome/plasmid. These file can be
downloaded from ftp://ftp.ncbi.nih.gov/genomes/Bacteria/ or
ftp://ftp.ncbi.nih.gov/genbank/genomes/Bacteria/</p>

<p>For example: To analyse mapped date of E. coli K12 MG1655
the files NC_000913.fna and NC_000913.gff are
downloaded.  Unfortunately some manual work is needed to modify the
FASTA so that IGB can handle it. At first the ending has to be changed
to ".fa". Additionally the head of the Fasta file has to be
shorted. In the example we change it to "&gt;NC_000913.2". Once the
file changes are done we can start the IGB and load the ggf and the
FASTA by clicking "File" then "Open File", then selecting the two
files followed by clicking "Open". After the reference files are load
the GR files can be added the same way. It makes sense to only load
normalized and or not-normalized GR files together. The next step
should be the adjustment of the displayed intesities. If you select a
library display by clicking on it (you can select more than one by
holding Shift) and then select the "Graph Adjuster" tab you can set
the minimal and maximal values shown. For plus strand GR files a
minimum of 0 and a maximum of 100 is good starting point. For minus
strand GR files a minimum of -100 and a maximum of 0 should be
set. Try different value to see which serves your purpose best. You
can rearrange the order of the libs by shifting them up or down. Once
you have set up the libraries as described you can start exploring
your data.</p>

<h2><a name="read_quantification">Feature based read quantification</a></h2>
<p><a href="#toc">top</a></p>

Could be gene expression.

<h2><a name="report">Reports and statistics</a></h2>
<p><a href="#toc">top</a></p>

<a href="./RAPL_analysis/output/reports_and_stats/">stats</a>

</body></html>
"""
