#' @title Get full text of PMC papers
#' @import tidypmc dplyr
#'
#' @param PMC_list A character vector of PMCs.
#' @param Collapse Need to collapse sections into a single text?
#'
#' @return List of full texts of \code{PMC_list}.
#' @export
full_text <- function(PMC_list, Collapse = F){
  get_full_text <- function(PMC){
    doc <- pmc_xml(PMC)
    txt <- pmc_text(doc)

    if (Collapse)
    {
      text_sec <- txt %>%
        summarise(text=paste(text,collapse='')) %>%
        pull(text)
    } else{
      text_sec <- txt %>%
        group_by(section) %>%
        summarise(text=paste(text,collapse=''))
    }

  }

  full_text_out <- lapply(PMC_list, get_full_text)
  names(full_text_out) <- PMC_list
  full_text_out
}



#' @title Get pictures captures PMC papers
#' @import tidypmc dplyr
#'
#' @param PMC_list A character vector of PMCs.
#' @param Collapse Need to collapse sections into a single text?
#'
#' @return List of pictures captures of \code{PMC_list}.
#' @export
pic_captures <- function(PMC_list, Collapse = F){
  get_pic_captures <- function(PMC){
    doc <- pmc_xml(PMC)
    cap <- pmc_caption(doc)

    if (Collapse)
    {
      cap_sec <- cap %>%
        summarise(text=paste(text,collapse='')) %>%
        pull(text)
    } else{
      cap_sec <- cap %>%
        group_by(label) %>%
        summarise(text=paste(text,collapse=''))
    }

  }

  pic_captures_out <- lapply(PMC_list, get_pic_captures)
  names(pic_captures_out) <- PMC_list
  pic_captures_out
}




#' @title Collapse uncollapsed dfs from \code{full_text} or \code{pic_captures}
#' @import dplyr
#'
#' @param dfs_list A character vector of PMCs.
#'
#' @return List of pictures captures of \code{PMC_list}.
#' @export
collapse_df <- function(dfs_list){

  collapse_df <- function(df){
    cap_sec <- df %>%
      summarise(text=paste(text,collapse='')) %>%
      pull(text)

  }

  collapsed_out <- lapply(dfs_list, collapse_df)
  names(collapsed_out) <- dfs_list
  collapsed_out
}
