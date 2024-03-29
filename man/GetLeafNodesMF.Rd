% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allgotree.R
\name{GetLeafNodesMF}
\alias{GetLeafNodesMF}
\title{All the MF GO-DAG leaf nodes}
\usage{
GetLeafNodesMF(organism = NULL)
}
\arguments{
\item{organism}{An object of class character that represents an organism.}
}
\value{
A two-coloum matrix of all the leaf nodes in the organism-specific DAG and their respective levels
}
\description{
Returns all the leaf nodes from a particular organism-specific GO-DAG. If empty (or if "MF"), all
             the general MF GO-DAG leaf nodes are returned. The supported organisms are "Homo sapiens / Human",
             "Rattus norvegicus / Rat", "Mus musculus / Mouse", "Danio rerio / Zebrafish", "Caenorhabditis elegans / Worm",
             "Arabidopsis thaliana / Cress", "Saccharomyces cerevisiae / Yeast", "Schizosaccharomyces pombe / Fission yeast",
             "Drosophila melanogaster / Fruit fly", "Escherichia coli / E.coli".
}
\examples{
\donttest{
# Mouse MF GO-DAG leaf nodes
GetLeafNodesMF(organism = "Human")

# Arabidopsis thaliana MF GO-DAG leaf nodes
GetLeafNodesMF(organism = "Arabidopsis thaliana")

# Drosophila melanogaster MF GO-DAG leaf nodes
GetLeafNodesMF(organism = "Drosophila melanogaster")
}
}
