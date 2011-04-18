import csv
from datetime import datetime
from string import Template

from rapl.inputstats import InputStats
from rapl.parameters import Parameters
from rapl.paths import Paths

class Reporter(object):

    def __init__(self, rapl_instance):
        """Create instance

        Arguments:
        - `rapl_instance`: The rapl instance that is reported about
                           
        """
        self.paths = Paths()
        self.parameters = Parameters()
        self.input_stats = InputStats()


    def report(self):
        """
        """
        template = Template(self._template_text())
        return(template.substitute(
                date = datetime.today().strftime("%d.%m.%Y"),
                input_folder = self._latex_safe(
                    self.paths.input_folder),
                RNA_lib_folder = self._latex_safe(
                    self.paths.rna_seq_folder),
                genome_file_folder = self._latex_safe(
                    self.paths.genome_folder),
                annotation_file_folder = self._latex_safe(
                    self.paths.annotation_folder),
                output_folder = self._latex_safe(
                    self.paths.output_folder),
                read_lib_listing = self._read_lib_listing(),
                tracing_summary = self._latex_safe(
                    self.paths.tracing_summary_file),
                rna_seq_file_stats = self._latex_safe(
                    self.paths.read_file_stats),
                trace_summary_table = self._trace_summary_table(),
                lib_genome_summary_table = (
                    self._lib_genome_read_summary_table()),
                lib_genome_summary_file = self._latex_safe(
                    self.paths.lib_genome_read_mapping_summary),
                genome_file_listing = self._genome_file_listing(),
                index_file_name = self._latex_safe(
                    self.paths.segemehl_index_name()),
                index_folder = self._latex_safe(
                    self.paths.read_mapping_index_folder),
                min_seq_length = self.parameters.min_seq_length,
                max_a_content = self.parameters.max_a_content,
                read_mapping_folder = self._latex_safe(
                    self.paths.read_mappings_first_run_folder),
                read_mapping_folder_second_run = self._latex_safe(
                    self.paths.read_mappings_second_run_folder),
                unmapped_read_first_run_folder = self._latex_safe(
                    self.paths.umapped_reads_of_first_mapping_folder),
                unmapped_read_second_run_folder = self._latex_safe(
                    self.paths.umapped_reads_of_second_mapping_folder),
                combined_mappings_folder = self._latex_safe(
                    self.paths.combined_mappings_folder),
                tracs_file_folder = self._latex_safe(
                    self.paths.read_tracing_folder),
                combined_mapping_split_folder = self._latex_safe(
                    self.paths.combined_mapping_split_folder),
                genome_file_stats = self._latex_safe(
                    self.paths.genome_file_stats),
                gr_folder = self._latex_safe(
                    self.paths.gr_folder)
               ))
    
    def _trace_summary_table(self):
        raw_table = open(self.paths.tracing_summary_file).read()
        row_counter = 0 
        table_header = []
        table_data = []
        for row in raw_table.split("\n"):
            if row_counter == 0:
                table_header = row[1:].split("\t")
            else:
                table_data.append(row.split("\t"))
            row_counter += 1
        table_string = "\\begin{tabular}{l%s}\n" % ("r" * (len(table_header)-1))
        table_string += "%s\\\\\n" % (" & ".join(
                ["\\rotatebox{90}{%s}" % head_item.replace("_", " ")
                 for head_item in table_header]))
        for row in table_data:
            table_string += "%s\\\\\n" % (
                " & ".join([self._latex_safe(field) for field in row]))
        table_string += "\\end{tabular}\n"
        return(table_string)

    def _lib_genome_read_summary_table(self):
        rows = []
        for line in open(self.paths.lib_genome_read_mapping_summary):
            rows.append(line[:-1].split("\t"))

        table_string = "\\begin{tabular}{l%s}\n" % ("r" * (len(rows[0])-1))

        table_string += " & ".join(
            ["\\rotatebox{90}{%s}" % (self._latex_safe((cell))) 
             for cell in rows[0]]) + "\\\\\n"
        for row in rows[1:]:
                    table_string += " & ".join(
                        [self._latex_safe(cell) for cell in row]) + "\\\\\n"
        table_string += "\\end{tabular}\n"
        return(table_string)

    def _read_lib_listing(self):
        """ """
        listing = "\\begin{itemize}\n"
        for read_file in self.paths.read_files:
            listing += "\\item %s with %s reads " % (
                self._latex_safe(read_file),
                self.input_stats._number_of_fasta_entries(
                    self.paths.read_file(read_file)))
        listing += "\\end{itemize}\n"
        return(listing)

    def _genome_file_listing(self):
        listing = "\\begin{itemize}\n"
        listing += "\n".join(
            ["\\item %s\n" % self._latex_safe(genome_file)
             for genome_file in self.paths.genome_files])
        listing += "\\end{itemize}\n"
        return(listing)

    def _latex_safe(self, input_string):
        return(input_string.replace("_", "\_"))

    def _template_text(self):
        return("""\\documentclass[12pt,a4paper,oneside]{article}
\\usepackage[latin1]{inputenc}
\\usepackage{graphicx}
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

This report intends to help the reader to understand the results of a
Solexa or 454 cDNA sequencing (also known as RNA-seq) and analyzing
process.\\\\

At first: Relax - advanced skills in bioinformatics/computation
biology are not required. Most of the files are plain text files which
can be opened with simple text viewer/ editor programs like
\\texttt{less}, \\texttt{gedit}, \\texttt{emacs}, \\texttt{vi} in
GNU/Linux, Mac OS and other Unix-like operating systems. In Windows
\\texttt{notepad} is a sufficient tool but the files can also be
loaded in more sophisticated office applications like
\\texttt{OpenOffice Doc} or \\texttt{Microsoft Word}. Some files are
so called CSV files (character separeted files) in which columns are
separated by tabs. These files can also be loaded in spreadsheet
application like \\texttt{OpenOffice Calc} or \\texttt{Microsoft
Excel} and then saved in the native format of the application. The
graphics are stored as PDF files and are generated using the
programming language \\texttt{R}. The names of the R scripts end with
\\texttt{.R} and can be modified with the above mentioned text
editors.

\\subsection{Organization of files and folders}

To make the navigation to required files and the generated output easy
and intuitive there are two folder in the root of the project
folder. The folder \\texttt{$input_folder} contains the files required
for the run. In the subfolder \\texttt{$RNA_lib_folder} the read files
are located, the subfolder \\texttt{$genome_file_folder} contains the
genome files (i.e. the FASTA files of all chromosom and plasmids)
while the subfolder \\texttt{$annotation_file_folder} harbours
annotation files.\\\\

The created files are stored in subfolder of the folder
\\texttt{$output_folder}. Detailed description of theses subfolder and
their content can be found in the following sections.

\\section{Mapping the reads to the reference genome}

\\subsection{Results}

For the impatient among us who do not want to be bothered with details
- here is the mapping result overview:\\\\
{\\tiny
$trace_summary_table
}
The data of this table can also be found in the file\\\\
\\texttt{$tracing_summary}.\\\\

Number of successfully mapped reads by read library and genome file:\\\\
{\\tiny
$lib_genome_summary_table
}
\\ \\\\
The data of this table can also be found in the file\\\\
\\texttt{$lib_genome_summary_file}.\\\\


\\subsection{Description of the mapping process}

The mapping against the reference genome was done using the program
\\texttt{segemehl} (Hoffmann \\textit{et. al}, \\textit{PLoS
Computational Biology}, 2009) and took place in two steps. In the
first step it was tried to align all raw reads to the reference
genome. Usually this is not successful for all of them e.g due to
sequencing errors. Some read cannot be mapped as they contain a poly-A
tail. To tackle this problem all reads for which the mapping failed in
the first run were undergoing a procedure in which a potential poly-A
tail was searched and clipped. After the clipping the resulting
sequences were filtered by sequence size: All sequences that were
shorter than $min_seq_length nucleotides were removed. The other
remaining reads were aligned to the genome in a second round. After
this the mapping results of the two runs were combined. These reads of
these combined files were then filtered by their A-content. All
entries with more than $max_a_content\\% of A's were removed. These
cleaned results are then split by the genome files (of chromosome
and/or plasmids) they are mapped to.\\\\

If you look at the mapping results you might noticed that some reads
are listed more than once. We allowed only one best \\texttt{segemehl}
hit per read but for certain reads there were some hits sharing the
same best value. In such case two or more placements are accepted.

\\section{Used libraries}

The first step of the process is the actual sequencing of the cDNA
libraries. The result of this sequencing is a set of FASTA files which
contain the raw reads. The barcode sequenced that are used for
multiplexed runs (more than one lib per sequencing lane) are already
removed from these. For your project we used the following read files
stored in the folder \\texttt{$RNA_lib_folder}:

$read_lib_listing

Basic statistics (e.g. number of reads) about these file are stored in
the file \\texttt{$rna_seq_file_stats}.

\\subsection{Used genome files}

The following genome files were used as reference
set:

$genome_file_listing

Basic statistics (e.g. number of lines) about these file are stored in
the file \\texttt{$genome_file_stats}.\\\\


\\texttt{segemehl} generates one combined index files of them save as\\\\
\\texttt{$index_file_name} \\\\
in the folder \\texttt{$index_folder}. 

\\subsection{First read mapping}

The output files of the first mapping can be found in the subfolder
named \\texttt{$read_mapping_folder}.

\\subsection{Unmapped reads, poly-A clipping and size filtering}

The files which contain reads that were not mapped in the first run
are put in the folter \\texttt{$unmapped_read_first_run_folder}. The
clipped files based on these files are also there. Likewise the files
that contain the reads which are equal or greate the given cut-off
(they have a "gtoe" in the file name) and the files that contains the
reads that are too small after the clipping ((they have a "lt" in the
file name)).

\\subsection{Second read mapping}

The result files of the second mapping are located in the folder\\\\
\\texttt{$read_mapping_folder_second_run}.\\\\

Unmappable reads can be found in the files in folder\\\\
\\texttt{$unmapped_read_second_run_folder}.

\\subsection{Combining of the mapping results and filtering the
results by A-content}

The combined results of both read mapping steps are stored in the
files in the folder \\texttt{$combined_mappings_folder}. Also the files
which are generated in the A-content filtering step are stored here.

\\subsection{Tracing files}

For each library a CVS file is generated that traces the handling of each
read. It shows at which stage a read is filter out or mapped. These
files are stored in the folder \\texttt{$tracs_file_folder}.

\\subsection{Splitting result by genome files}

The file that contain the mapping of libraries split by target genome
file are stored in folder \\texttt{$combined_mapping_split_folder}.

\\section{Interactive visualization}

To visualize and compare the different libraries a genome browser can
be used. Usually \\texttt{Integrated Genome Browser} (IGB) is a
suitable tool for this purpose.\\\\

\\subsection{GR files}

Based on the combined and A-content filtered mapping results GR file
that can be viewed in the IGB were generated. These GR files describe
the amount of mapped reads per nucleotide position for a each genome
file. Above it was mentioned that some reads are equally well mapped
to different locations on a genome. In such a case we divide the raw
counting (i.e. 1) by the number of mappings. E.g. if a read is mapped
to 4 different places each position only get 1/4 point counted for
this read.\\\\

For each genome file (chromosome/plasmid) two different kind of sets
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
multiplication process results in the raw countings).\\\\

For each normalized/not-normalized library-genome-file-pair you have
two different GR files. One for the plus strand, one for the minus
strand. You can recognize the type by the file name. The plus strand
files end with "\\_plus.gr", the minus strand files with
"\\_minus.gr". While the plus strand file should contain only positive
intensities, the minus strand file should contain only negative
intensities.\\\\

The GR files are located in the folder \\texttt{$gr_folder}.

\\subsection{IGB usage in a nutshell}

The Integrated Genome Browser is a open-source application that can be
downloaded from http://www.bioviz.org/igb/. It requires Java.\\\\

To make sense of the coordinates of the GR files data of the reference
genome has to be supplied. This data is taken from a FASTA file (.fna)
and a gff file of the analysed chromosome/plasmid. These file can be
downloaded from ftp://ftp.ncbi.nih.gov/genomes/Bacteria/ or\\\\
ftp://ftp.ncbi.nih.gov/genbank/genomes/Bacteria/

For example: To analyse mapped date of \\textit{E. coli} K12 MG1655
the files \\texttt{NC\_000913.fna} and \\texttt{NC\_000913.gff} are
downloaded.  Unfortunately some manual work is needed to modify the
FASTA so that IGB can handle it. At first the ending has to be changed
to ".fa". Additionally the head of the Fasta file has to be
shorted. In the example we change it to "$$>$$NC\_000913.2". Once the
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
your data.

%\\section{Counting of annotation overlaps}

\\end{document}
""")
