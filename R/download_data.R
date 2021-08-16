##### download_data #######################################################################3

#' @title Download "cited_in" or "references" data from PubMed
#' @description Download data about papers that cites input paper from PubMed by its PMID.
#'
#' @importFrom rentrez entrez_link
#' @importFrom reshape2 melt
#'
#' @param PMID A single character or integer of interests paper.
#' @param type One of "cited_by" or "references" to specify output info.
#'
#' @return If available data is not \code{NULL} then returns corresponding to
#' \code{type} character vector.
#' @export
download_data <- function(PMID, type = c("cited_by", "references"))
{
  cit_ref_data <- entrez_link(db='pubmed', dbfrom='pubmed',
                         retmode='xml', id=PMID,
                         cmd='neighbor')$links
  if (type == "cited_by"){
    citedby <- cit_ref_data$pubmed_pubmed_citedin

    if (!is.null(citedby)){
      return(citedby)
    }else{
      return(NULL)
    }
  }
  else{
    if (type == "references"){
      references_ <- cit_ref_data$pubmed_pubmed_refs
      if (!is.null(references_)){
        return(references_)
      }else{
        return(NULL)
      }
    }
    else
    {
      stop("Invalid type")
    }
  }



}


##### cited_by ########################################################################################33

#' @title Cited_by function
#' @description Download data about papers that cites given \code{PMID_list} from PubMed
#' by PMID character vector with \code{\link{download_data}} and returns data.frame.
#'
#' @importFrom dplyr bind_rows
#'
#' @param PMID_list A character vector of PMIDs of interests papers.
#'
#' @return If available data is not \code{NULL}.
#' then returns data.frame with data about papers that cites \code{PMID_list}.
#' @export
cited_by  <- function(PMID_list)
{
  # Добавить проверку на NULL выдачи download_data
  list_of_citation_by_matrices <- sapply(PMID_list, download_data, "cited_by")
  list_of_citation_by_matrices <- list_of_citation_by_matrices[lengths(list_of_citation_by_matrices) != 0]

  citation_by_table <- melt(list_of_citation_by_matrices)
  colnames(citation_by_table) <- c("pmid_cited_by", "parent_work")

  return(citation_by_table)
}






##### references ############################################################################################

#' @title Cited_by function
#' @description Download data about papers on which refers given \code{PMID_list} from PubMed
#' by PMID character vector with \code{\link{download_data}} and returns data.frame.
#'
#' @importFrom dplyr bind_rows
#'
#' @param PMID_list A character vector of PMIDs of interests papers.
#'
#' @return If available data is not \code{NULL}.
#' then returns data.frame with data about papers on which refers \code{PMID_list}.
#'
#' @export
references  <- function(PMID_list)
{
  # Добавить проверку на NULL выдачи download_data
  list_of_references_matrices <- sapply(PMID_list, download_data, "references")
  list_of_references_matrices <- list_of_references_matrices[lengths(list_of_references_matrices) != 0]

  references_table <- melt(list_of_references_matrices)
  colnames(references_table) <-c("pmid_reference", "parent_work")


  return(references_table)
}
