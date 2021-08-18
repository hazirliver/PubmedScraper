#' @title Main function
#'
#' @param PMID_list A character vector of PMIDs of interesting papers.
#' @param filename_base A character contains path to output and filename base
#' @param venv A character specifies your python's venv in which installed biopython and pandas.
#' @param n_cites_threshold An integer. The minimum number of papers that a given article should refer to.
#' @param n_refers_treshold An integer. The minimum number of papers that a given article should refer to.
#'
#' @return List of two data.frame's. First -- about papers that cites given \code{PMID_list}, second -- about papers on which refer given \code{PMID_list}.
#' @export
#'
scrap_pubmed <- function(PMID_list, filename_base, venv, n_cites_threshold = 2, n_refers_treshold = 2)
{
  return(check_pcg())
  #cit_tab <- graph_foo_citation(PMID_list, filename_base, n_cites_threshold)
  #cit_tab_medline <- add_medline(cit_tab, venv)



  #references_tab <- graph_foo_references(PMID_list, filename_base, n_refers_treshold)
  #references_tab_medline <- add_medline(references_tab, venv)

  return(list(cit_tab_medline, references_tab_medline))

}
