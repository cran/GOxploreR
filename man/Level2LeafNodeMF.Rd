% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/degreeDistMF.R
\name{Level2LeafNodeMF}
\alias{Level2LeafNodeMF}
\title{GO molecular function (MF) level leaf nodes getter}
\usage{
Level2LeafNodeMF(level, organism = NULL)
}
\arguments{
\item{level}{A numeric value for the GO-level}

\item{organism}{An organism of interest from the list of supported organism. If NULL the results are from the general
GO-DAG (default).}
}
\value{
All leaf nodes on the level
}
\description{
Derive all the leaf nodes from a GO MF level. The supported organism are "Homo sapiens / Human",
             "Rattus norvegicus / Rat", "Mus musculus / Mouse", "Danio rerio / Zebrafish", "Caenorhabditis elegans / Worm",
             "Arabidopsis thaliana / Cress", "Saccharomyces cerevisiae / Yeast", "Schizosaccharomyces pombe / Fission yeast",
             "Drosophila melanogaster / Fruit fly", "Escherichia coli / E.coli".
}
\examples{
# level 9 leaf nodes for Human GO-DAG
Level2LeafNodeMF(9, "Human")

# level 5 leaf nodes for Mouse GO-DAG
Level2LeafNodeMF(5, "Mouse")

}
