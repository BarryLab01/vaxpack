\name{vaxpack_output}
\alias{vaxpack_output}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Accessing Results From Populuation Genetic Analyses
%%  ~~function to do ... ~~
}
\description{vaxpack_output() is to access the results from a range of selection
from table, graphs, and sliding window scale. Users can define minimial haplotype
frequency, and minimal allele frequency to calcuate diversity and neutrality statistics.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vaxpack_output()
}
%- maybe also 'usage' for other objects documented here.

}
\details{Haplotype diversity, nucleotide diversity (Nei & Li, 1979),
and Tajima's D (Tajima, F, 1989) are used to calculate neutrality and
diversity statistic. The significant test for Tajima's D are determined
by original computer simulation from Fumio Tajima on a variety of sample
size (see more on reference). Graphical output are plotted on ggplot2. Input DNA sequences are automatically translated into amino acid using implemented function. So,the coding sequence is expected for the input for accurate reading frame.
Phylogenetic tree is drawn based on neighbour joining (NJ), unrooted method using ape package, and is only intended for initial exploration of samples, and therefore bootstrapped values are not included at the moment. Uses can define
their own sliding window, and step size. Option-5 from vaxpack output,
"Haplotype Population Pie chart" is plotted using plotly package.
Option-12 from vaxpack output "All Sliding Scale Graphs Overlapped" is scaled to Tajima's D. Haplotype accumulation plot is calculated based on "rarefaction method" with 100 permutations using imported specaccum() function from the vegan package.
%%  ~~ If necessary, more details than the description above ~~
}
\references{1. Gotellli, N.J. & Colwell, R.K. (2001). "Quantifying biodiversity: procedures and pitfalls in measurement and comparison of species richness". Ecology Letters 4, 379–391.

2. Nei, M.; Li, W. (1979). "Mathematical Model for Studying Genetic Variation in Terms of Restriction Endonucleases". PNAS. 76 (10): 5269–73. doi:10.1073/pnas.76.10.5269. PMC 413122. PMID 291943.

3. Nei, M.; Tajima, F. (1981), "DNA polymorphism detectable by restriction endonucleases", Genetics 97:145

4. Tajima, F. (1989). "Statistical method for testing the neutral mutation hypothesis by DNA polymorphism". Genetics. 123 (3): 585–95. PMC 1203831. PMID 2513255



%% ~put references to the literature/web site here ~
}
\author{Elijah Martin, Myo T. Naung (Barry Laboratory)
%%  ~~who you are~~
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{vaxpack_input()}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\note{The results can be stored in .csv or excel formats.
%%  ~~further notes~~
}
\examples{
vaxpack_output()
Enter a number from the selection above - 1
Minimal haplotypes will be calculated using only the segregation sites where
polymorphism is found in at least 'x' percent of the population at that site
'x' <- 1
-------------------------------------
vaxpack_output() #for more options
}

