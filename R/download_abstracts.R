#' Downloading abstracts by PMID list
#'
#' @param PMID_list A character vector PMIDs
#'
#' @importFrom rentrez entrez_fetch
#' @importFrom XML xpathApply xmlValue xmlChildren
#'
#' @return data.frame
#' @export
#'
download_abstracts <- function(PMID_list){
  PMID_list <- unique(PMID_list)
  fetch.pubmed <- entrez_fetch(db = "pubmed", id = PMID_list,
                               rettype = "xml", parsed = T)

  abstracts <-  xpathApply(fetch.pubmed, '//PubmedArticle//Article', function(x)
    xmlValue(xmlChildren(x)$Abstract))


  names(abstracts) <- your.ids


  abstracts.df <- stack(abstracts)
  colnames(abstracts.df) <- c("Abstract", "PMID")

  return(abstracts.df[c("PMID", "Abstract")])
}
