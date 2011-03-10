from string import Template
from datetime import datetime
class RaplReporter(object):

    def __init__(self, rapl_instance):
        """Create instance

        Arguments:
        - `rapl_instance`: The rapl instance that is reported about
                           
        """
        self.rapl_instance = rapl_instance

    def report(self):
        """
        """
        template = Template(self._template_text())
        return(template.substitute(
              date = datetime.today().strftime("%d.%m.%Y"),
              read_lib_listing = self._read_lib_listing(),
              genome_file_listing = self._genome_file_listing(),
              min_seq_length = self.rapl_instance.min_seq_length,
              max_a_content = self.rapl_instance.max_a_content,
              read_mapping_folder = self._latex_safe(
                    self.rapl_instance.read_mapping_folder),
              unmapped_read_first_run_folder = self._latex_safe(
                self.rapl_instance.umapped_reads_of_first_mapping_folder)))

    def _read_lib_listing(self):
        """ """
        listing = "\\begin{itemize}\n"
        for read_file in self.rapl_instance.read_files:
            listing += "\\item %s with %s reads " % (
                self._latex_safe(read_file),
                self.rapl_instance._number_of_fasta_entries(
                    self.rapl_instance._read_file_path(read_file)))
        listing += "\\end{itemize}\n"
        return(listing)

    def _genome_file_listing(self):
        listing = "\\begin{itemize}\n"
        listing += "\n".join(
            ["\\item %s\n" % self._latex_safe(genome_file)
             for genome_file in self.rapl_instance.genome_files])
        listing += "\\end{itemize}\n"
        return(listing)

    def _latex_safe(self, input_string):
        return(input_string.replace("_", "\_"))

    def _template_text(self):
        return("""\\documentclass[12pt,a4paper,oneside]{article}
\\usepackage[latin1]{inputenc}
\\usepackage[pdfborder={0 0 0}]{hyperref}
\\setlength{\\parindent}{0ex}

\\begin{document}
\\begin{titlepage}
  \\title{\\huge Report}
  \\date{$date}
  \\author{}
  \\maketitle
\\end{titlepage}
\\pagestyle{plain}

\\section{Introduction}

This report intends to help you to understand the results of a Solexa
or 454 sequncing and analyzing process.

\\section{Used libraries}

The first step of the process is the actual sequencing of the cDNA
libraries. The result of this sequencing is a set of FASTA files
which contain the raw reads. For your project we used the following read
files:

$read_lib_listing

\\section{Read mapping}

\\subsection{Description of the mapping process}

The mapping against the reference genome was done using the program
\\texttt{segemehl} (Hoffmann \\textit{et. al}, \\textit{PLoS
Computational Biology}, 2009) and in two steps. In the first step it
was tried to align all raw reads to the reference genome. Usually this
is not successful for all of them e.g due to sequencing errors. Some
read cannot be mapped as they contain a poly-A tail. To tackle this
problem all reads for which the mapping failed in the first run were
undergoing a procedure in which a potential poly-A tail was searched
and clipped. After the clipping the resulting sequences were filtered
by sequence size: All sequences that were shorter than $min_seq_length
nucleotides were removed. The other remaining reads were aligned to
the genome in a second round. After this mapping results were
combined. These combined entries were then filtered by the
A-content. All entries with more than $max_a_content\\% of A's were
removed. These cleaned results are then split by the genome files (of
chromosome and/or plasmids) they are mapped to.\\\\

If you look at the mapping results you might noticed that some reads
are listed more than once. We allowed only one best \\texttt{segemehl}
hit per read but for certain reads there were some hits sharing the
same best value. In such case two or more placements are accepted.

\\subsection{Used Genome files}

The following genome files were used as reference
set:

$genome_file_listing

\\texttt{segemehl} generates one combined index files of them.

\\subsection{First read mapping}

The output files of the first mapping can be found in the subfolder
named \\textit{$read_mapping_folder}.

\\subsection{Unmapped reads, poly-A clipping and size filtering}

The files which contain read that were not mapped in the first run are
put in the folter \\textit{$unmapped_read_first_run_folder}. The
clipped files based on these files are also there. Likewise the files
that contain the reads which are equal or greate the given cut-off
(they have a "gtoe" in the file name) and the file that contains the
reads that are too small after the clipping ((they have a "lt" in the
file name)).

\\subsection{Second read mapping}

\\subsection{Combining of the mapping results}

\\subsection{Filtering the results by A-content}

\\subsection{Tracing files}

\\subsection{Splitting result by genome file}

\\section{Preparing for visualization}

To visualize and compare the different libraries a genome browser can
be used. Usually \\texttt{Integrated Genome Browser} (IGB) is a
suitable tool for this purpose.\\\\

\\subsection{GR files}

Based on the combined and filtered mapping results GR file that can be
viewed in the IGB were generated. These GR files describe the amount
of mapped reads per nucleotide position for a each genome file. Above
it was mentioned that some reads are equally well mapped to different
locations on a genome. In such a case we divide the raw counting
(i.e. 1) by the number of mappings. E.g. if a read is mapped to 4
different places each position only get 1/4 point counted for this
read.\\\\

For each genome file (chromosome/plasmid) two different kind of sets of
GR files were generated. One set represents the raw countings of
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
1000,000 and then multiply by 1000,000 (yes, for the libs with the
lowest number of mapped reads the whole normalization and
multiplication process results in the raw countings).\\\\

For each normalized/not-normalized library-genome-file-pair you have
two different GR files. One for the plus strand, one for the minus
strand. You can recognize the type by the file name. The plus strand
files end with "\\_plus.gr", the minus strand files with
"\\_minus.gr". While the plus strand file should contain only positive
intensities, the minus strand file should contain only negative
intensities.

\\subsection{IGB usage in a nutshell}

ftp://ftp.ncbi.nih.gov/genomes/Bacteria/ \\\\
ftp://ftp.ncbi.nih.gov/genbank/genomes/Bacteria/

\\section{Counting of annotation overlaps}

\\end{document}
""")
