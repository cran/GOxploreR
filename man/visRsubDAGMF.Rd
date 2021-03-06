% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visdag.R
\name{visRsubDAGMF}
\alias{visRsubDAGMF}
\title{Visualise the GO-DAG for a specific organism based on certain GO-terms}
\usage{
visRsubDAGMF(goterm, organism)
}
\arguments{
\item{goterm}{A vector of molecular function GO-terms}

\item{organism}{The organism whose DAG we want to visualise based on the GO-terms supplied. If this argument is "MF"
the general reduced GO MF tree is used.}
}
\value{
Returns a plot highlighting the important categories i.e. the categories containing the GO-terms of interest
}
\description{
A category in the organism-specific DAG is regarded as important if it contains at least one GO-terms from
             the vector supplied. The function does not visualise the entire organism's DAG, it accepts a range of GO-terms
              and DAG categories that contains
             these GO-terms are visualise. We refer to these categories as the important and everything else
             is faded out. The supported organisms are "Homo sapiens / Human", "Rattus norvegicus / Rat",
             "Mus musculus / Mouse", "Danio rerio / Zebrafish",
             "Caenorhabditis elegans / Worm" ,"Arabidopsis thaliana / Cress", "Saccharomyces cerevisiae / Yeast",
              "Schizosaccharomyces pombe / Fission yeast",
             "Drosophila melanogaster / Fruit fly", "Escherichia coli / E.coli".
}
\examples{
Terms <- c("GO:0034040","GO:0008374" ,"GO:0003777","GO:0003674","GO:0015440","GO:0005464")
visRsubDAGMF(Terms, "Human")
}
