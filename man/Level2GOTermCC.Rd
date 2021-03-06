% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gotermlevelcc.R
\name{Level2GOTermCC}
\alias{Level2GOTermCC}
\title{GO cellular component (CC) terms on a level}
\usage{
Level2GOTermCC(level, organism = NULL)
}
\arguments{
\item{level}{A numeric value for the GO-level}

\item{organism}{The organism of interest. If NULL the results will be from the general GO CC tree}
}
\value{
A list of GO-terms from the level
}
\description{
Gets all the biological process (CC) GO-terms on a GO-level.
             The supported organisms are "Homo sapiens / Human", "Rattus norvegicus / Rat", "Mus musculus / Mouse",
             "Danio rerio / Zebrafish", "Caenorhabditis elegans / Worm" ,"Arabidopsis thaliana / Cress",
             "Saccharomyces cerevisiae / Yeast", "Schizosaccharomyces pombe / Fission yeast",
             "Drosophila melanogaster / Fruit fly", "Escherichia coli / E.coli".
}
\note{
The Gene Ontology (GO) cellular component (CC) tree was built using the root node GO:0005575
}
\examples{
# Gene association GO-terms for organism Rat on level 4
Level2GOTermCC(level = 4, organism = "Rat")

# Gene association GO-terms for the General GO CC tree
Level2GOTermCC(level = 7)

}
