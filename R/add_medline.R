#' @title Adding additional Medline info
#' @description Supplements the \code{data_table} with additional information from the Medline database.
#'
#' @importFrom utils read.csv write.table
#' @importFrom reticulate use_python source_python
#' @importFrom dplyr mutate
#'
#' @param data_table A data.frame that will supplement with additional Medline information.
#'
#' @return A data.frame with merged additional Medline info
#' @export
#'
add_medline <- function(data_table)
{
  local_PMID <- data_table %>%
    pull(PMID)

  user_name <- system('echo "$USER"', intern = T)
  venv1 <- paste0("/home/", user_name, "/user_venv/bin/python")
  use_python(venv1, required = T)


  source_python(system.file("additional_inf.py", package = "PubmedScraper"))
  data_table_add_inf <- main_foo(local_PMID)


  # Загружаем полученную информацию и мерджим как дополнительные поля к нашей таблице
  data_table_with_add_inf <- merge(data_table, data_table_add_inf)

  data_table_with_add_inf <- data_table_with_add_inf %>%
    mutate(DP = anydate(DP),
           link_pmid = paste0("https://pubmed.ncbi.nlm.nih.gov/", as.character(PMID)),
           link_pmc = ifelse(as.character(PMC) == 'NaN', "",
                             paste0("https://www.ncbi.nlm.nih.gov/pmc/articles/", as.character(PMC))))





  return(data_table_with_add_inf)
}
