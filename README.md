# Vaxpack: Fast and Efficient Population Genetic Analysis Gadget
Elijah Martin, Myo Naung
January, 2019

Introduction
------------

Extensive use of sequencing technology to solve population related
biological problems demands efficient, and comphrehensive tools for data
analyses. Different software packages for population genetics/genomics
such as PopGenome, R package (Pfeifer, B. et al 2014) necessitate the
basic knowledge of command line tools. Programs with the graphical user
interface such as DnaSP require extensive pre-processing and
post-processing of the data. Here, **vaxpack** is a powerful
user-friendly R package for fast, easy and comphrehensive analyses of
population genomic/genetic data for a gene of interest of haploid
organisms. It is designed for biologists with minimal command line tool
knowledge. **vaxpack** includes a wide range of polymorphism and
population genetic diversity metrics, automatic translation of input DNA
sequences into the amino acid sequences, automatic visualization of
output on ggplot2, diversity and neutrality statistics with users'
defined minor haplotype/ allele frequency thresholds, as well as users'
defined sliding scale options. In addition, the output can be saved in
users’ defined formats for further manipulations, and analyses.
Graphical outputs are automatically labelled, and are ready to use.

The following sections explain how to use the vaxpack R package.

Installing vaxpack
------------------

The package can be installed directly from the Github using
**githubinstall()** function from the **githubinstall package** using
the packagename. vaxpack is built for R (&gt;= 3.1.0). See more about
helpful ways to install R packages from Github
[here](https://cran.r-project.org/web/packages/githubinstall/vignettes/githubinstall.html).

    install.packages("githubinstall")
    library(githubinstall)
    githubinstall("vaxpack")
    library(vaxpack)

Loading Data
------------

In the meantime, vaxpack only accepts aligned genomic (DNA) sequences of
the same length in “.fasta”, ".fas", ".fa", and ".seq" extensions as
input. Therefore, preprocessing of input data is required. Flexible
input formats such as from VCF will be included in the future release.
Coding sequences without introns are expected for accurate translation
of input into amino acids if the interest includes population genetic
analyses of amino acid sequences.

The input data has to be inside a folder similar to **readData()** from
PopGenome package. If the analyses include more than one populations of
the same gene/coding region, build different fasta files and put inside
a folder. 

**vaxpack\_input()** is an all in one function for conducting population
genetic analyses. It uses a prompt/response method for input, and needs
the following

1.  A file path to a folder containing aligned DNA fasta files to be
    analysed.

2.  A file path to a reference in fasta format (not inside folder).

3.  The name of the gene you are analysing for graphical purposes.

<!-- -->

    > vaxpack_input()
    The following analyses are intended for haploid organisms such as Plasmodium
    For this function to work you will need:
    1 - A folder containing all files to be analysed
         These need to be aligned, the same length, and with no gaps! (e.g. with MEGA)
         Only accepts bases 'A/a', 'T/t', 'C/c', 'G/g' 
         Gaps, or bases that are not AaTtCcGg, are called as reference, so accuracy is lost
         These replacements will be counted in the results table as "invalid sites"

         Make a different file for different populations to compare them
         e.g. 'Asymptomatic.fasta, Symptomatic.fasta', or 'Brazil.fasta, Peru.fasta'

    2 - A reference file for the gene you are analysing
         (This cannot be in the same folder as the other files!)

    Accepted extensions are ".fasta", ".fas", ".fa", and ".seq"

    Enter the file path to your folder with a "/" at the end 
    e.g. "/users/me/documents/r files/my fasta folder/"

     My folder path <-"/Users/Desktop/test_vaxpack/Fasta/"
     
      Enter the file path to your reference file
     e.g. "/users/me/documents/r files/reference.fasta"
     My ref .fasta file path <- "/Users/Desktop/test_vaxpack/reference CDS_.fasta"
     
     What is the name of the gene youre analysing ? <- genename

    Completed! Step: 23 of 23  
    vaxpack_input() took 42 secs 
    Now use vaxpack_output() to get your results! 

Normally, input files of 2kb gene from 500 samples took less than one
minute. Now, core calculation is done for vaxpack, and is ready for
outputs. Bigger sample size and larger gene might take longer.

Obtaining Output
----------------

Outputs can be instantly accessed via **vaxpack\_output()** functions.

    vaxpack_output()

    Choose your output:
    TABLES
    1  - Results Table 
    2  - Haplotype Table
    3  - Minor Allele Frequency Table (Nucleotides)
    4  - Minor Allele Frequency Table (Amino Acids)

    GRAPHS
    5  - Haplotype Population Pie Chart
    6  - AA Variant Percentage Column Graph
    7  - Phylogenetic Tree
    8  - Haplotype Accumulation Plot

    SLIDING SCALE
    9  - Polymorphism
    10 - Tajimas D
    11 - Nucleotide Diversity
    12 - All Sliding Scale Graphs Overlapped 

    13 - Sliding Scale data table 

    Enter a number from the selection above - 1

Specific output from above can be chosen easily. If option - 1: Results
Table is chosen, users can further manipulate threshold/cut-off of the
specific parameter.For example,

    Minimal haplotypes will be calculated using only the segregation sites where polymorphism
    is found in at least 'x' percent of the population at that site
    'x' <-  1 

    Saved as "vp.RESULTS.TABLE", use write.csv() to save to excel 
    e.g. write.csv(vp.RESULTS.TABLE, file = "my.results.table.in.excel")

1 - **Results Table** gives us sample size, length of sequence, invalid
sites, number of segregation size, number of single nucleotide
polymorphisms (SNPs), average tajima's D (Tajima, F, 1989), average
nucleotide diversity (Nei & Li, 1979), number of total haplotype based
on nucleotide, minimal haplotypes based on amino acid, and user defined
cut-off points.

***If more than one population is included in the analyses, the
following outputs use total of all input populations. ***

2 - **Haplotype Table** gives the information for total haplotypes
compositions according to their frequencies, and differences from the
reference amino acid sequence. Haplotypes are listed in descending order
of frequency, with Rank 1 representing the most common haplotype.

3 - **Minor Allele Frequency Table (Nucleotides)** gives us polymorphic
sites, reference nucleotid at the site, and composition of minor and
major allele frequencies.

4 - **Minor Allele Frequency Table (Amino Acids)** gives us polymorphic
sites, reference amino acid at the site, and composition of minor and
major allele frequencies.

5 - **Haplotype Population Pie Chart** displays distribution of
haplotypes above a specific cut-off value. It is built using **plotly
package**. The size of the fragment reflects the relative frequency of
haplotype found in the population.

6 - **AA Variant Percentage Column Graph** displays publication-quality
amino acid changes and positionsabove a specific user defined cut-off
values.

7 - **Phylogenetic Tree** is calculated using unrooted, neighor joining
(NJ) method using distance matrix and **ape package**. Different input
populations will be labelled in different colors. See more information
about ape package [here](http://ape-package.ird.fr/).

*Note: Bootstrapped values are not displayed in the tree in the
meantime, but will be included in future releases.*

8 - **Haplotype Accumulation Plot** is calculated based on "rarefaction
method" with 100 permutations using imported **specaccum()** function
from the vegan package. For detail information about vegan package, [see
here](https://cran.r-project.org/web/packages/vegan/index.html).

9 - **Polymorphism**, the publication-quality plot displays polymorphic
sites across the gene of interest under users' defined sliding window
scale.

10 - **Tajima's D**, the publication-quality plot displays Tajima's D
values across the gene of interest under users' defined sliding window
scale. p-values to test neutral theory of molecular evolution are
determined by original computer simulation from Fumio Tajima on a
variety of sample size following beta distrubution (Tajima, F, 1989).

11 - **Nucleotide Diversity**, the publication-quality plot displays
nucleotide diversity values across the gene of interest under users'
defined sliding window scale.

12 - **All Sliding Scale Graphs Overlapped**, the overlapped plot for
polymorphisms, nucleotide diversity, and Tajima's D values. It is scaled
to Tajima's D values.

13 - **Sliding Scale data table**, table output for segregation sites,
nucleotide diversity, tajima's D under users' defined sliding window
scales.

Additional Aspects
------------------

vaxpack can accept multiallelic positions. vaxpack accept SNP data
formatted in fasta, but will not produce accurate statistical outputs.
All genetic diversity metrics from vaxpack will be particularly useful
for antigenic diversity study. Future release will expand to diploid
organisms.

Handling missing data (Gap and ambigious bases)
-----------------------------------------------

vaxpack currently only accepts bases 'A/a', 'T/t', 'C/c', 'G/g'. Gaps or
ambigious bases are set as reference, and they are documented as invalid
sites.

Session Info
------------

    > sessionInfo()
    R version 3.4.0 (2017-04-21)
    Platform: x86_64-apple-darwin15.6.0 (64-bit)
    Running under: OS X El Capitan 10.11.6

    Matrix products: default
    BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
    LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

    locale:
    [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods  
    [7] base     

    loaded via a namespace (and not attached):
    [1] compiler_3.4.0 tools_3.4.0    yaml_2.2.0  

References
----------

1.  Gotellli, N.J. & Colwell, R.K. (2001). "Quantifying biodiversity:
    procedures and pitfalls in measurement and comparison of species
    richness". Ecology Letters *4*, 379–391.

2.  Nei, M.; Li, W. (1979). "Mathematical Model for Studying Genetic
    Variation in Terms of Restriction Endonucleases". PNAS. 76 (10):
    5269–73. <doi:10.1073/pnas.76.10.5269>. PMC 413122. PMID 291943.

3.  Nei, M.; Tajima, F. (1981), "DNA polymorphism detectable by
    restriction endonucleases", Genetics 97:145

4.  Pfeifer, B. et al. (2014) "PopGenome: An Efficient Swiss Army Knife
    for Population Genomic Analyses in R". Molecular Biology and
    Evolultion 31(7): 1929-1936.<doi:10.1093/molbev/msu136>

5.  Tajima, F. (1989). "Statistical method for testing the neutral
    mutation hypothesis by DNA polymorphism". Genetics. 123 (3): 585–95.
    PMC 1203831. PMID 2513255
