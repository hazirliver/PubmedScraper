% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/similar_articles.R
\name{download_similar}
\alias{download_similar}
\title{Find similar articles via parsing pubmed}
\usage{
download_similar(PMID, year_left, year_right = substr(Sys.Date(), 1, 4))
}
\arguments{
\item{PMID}{Current PMID. Default - current year.}

\item{year_left}{Shows left time interval for output}

\item{year_right}{Shows right time interval for output}
}
\value{
A data.frame with similar articles
}
\description{
A new feature of Pabmed is used to find similar articles.
}
