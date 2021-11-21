#' @title Find similar articles via parsing pubmed
#' @description A new feature of Pabmed is used to find similar articles.
#'
#' @importFrom rvest read_html html_nodes html_text
#' @importFrom stringr str_split
#'
#' @param PMID Current PMID. Default - current year.
#' @param year_left Shows left time interval for output
#' @param year_right Shows right time interval for output

#' @return A data.frame with similar articles
#' @export
#'
download_similar <- function(PMID,
                             year_left,
                             year_right = substr(Sys.Date(), 1, 4)){

  page_base <- paste0("https://pubmed.ncbi.nlm.nih.gov/?format=pmid&filter=years.",
                      year_left, "-" , year_right,
                      "&size=200&linkname=pubmed_pubmed&from_uid=")

  page <- read_html(paste0(page_base, PMID))

  pmids <- page %>%
    html_nodes(".search-results-chunk") %>%
    html_text()

  pmids <- unlist(str_split(pmids, "\\r\\n", simplify = F))
  pmids <- pmids[pmids != PMID]

}



#' @title Similar_articles function
#' @description Download data about papers that are similar to given \code{PMID_list} from PubMed
#' by PMID character vector with \code{\link{download_similar}} and returns data.frame.
#'
#' @importFrom dplyr bind_rows
#' @importFrom reshape2 melt
#'
#' @param PMID_list A character vector of PMIDs of interests papers.
#' @param year_left Shows left time interval for output
#' @param year_right Shows right time interval for output. Default - current year.
#'
#' @return If available data is not \code{NULL}.
#' then returns data.frame with data about papers that are similar to \code{PMID_list}.
#' @export
similar_articles  <- function(PMID_list,
                              year_left,
                              year_right = substr(Sys.Date(), 1, 4))
{
  # Добавить проверку на NULL выдачи download_data
  list_of_similar_matrices <- sapply(PMID_list, download_similar, year_left, year_right)
  list_of_similar_matrices <- list_of_similar_matrices[lengths(list_of_similar_matrices) != 0]

  similar_matrices_table <- melt(list_of_similar_matrices)
  colnames(similar_matrices_table) <- c("pmid_similar", "parent_work")

  return(similar_matrices_table)
}
