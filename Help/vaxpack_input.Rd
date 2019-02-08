\name{vaxpack_input}
\alias{vaxpack_input}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Input For Population Genetic Analysis
%%  ~~function to do ... ~~
}
\description{vaxpack_input() is an all in one function for conducting population genetic analyse. It uses a prompt/response method for input, and needs the following:

(1) A file path to a folder containing aligned DNA fasta files to be analysed.

(2) A file path to a reference .fasta not in folder file.

(3) The name of the gene you are analysing for graphical purposes.

(4) Only accepts bases 'A/a', 'T/t', 'C/c', 'G/g'. Gaps or ambigious bases are set as reference, and they are documented as invalid sites in results table.

(5) Coding regions are preferred so the package can acutally translate into AA sequences.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\details{
The package accepts extensions such as ".fasta"", ".fas", ".fa", and ".seq".
If you want to compare different datasets of the same gene/coding region, make different fasta files and put inside a folder.
%%  ~~ If necessary, more details than the description above ~~
}
\value{vaxpack_output() to obtain results
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{Elijah Martin, Myo T. Naung (Barry Laboratory)
%%  ~~who you are~~
}
\note{vaxpack_input may be slow to process large numbes of sequences in a single file.
%%  ~~further notes~~
}

\seealso{
\code{\link{vaxpack_output()}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\examples{
library(vaxpack)
vaxpack_input()
}



