#' @title Adding additional Medline info
#' @description Supplements the \code{data_table} with additional information from the Medline database.
#'
#' @importFrom utils read.csv write.table
#'
#' @param data_table A data.frame that will supplement with additional Medline information.
#' @param output_filename_base A character same as filename_base from \code{\link{scrap_pubmed}}
#' @param type Either cited_by either references to specify data.frame.
#'
#' @return A data.frame with merged additional Medline info
#'
#'
add_medline <- function(data_table, output_filename_base, type)
{
  data_table %>%
    pull(PMID) %>%
    write.table(paste0(output_filename_base, "_", type, ".lst"), quote = F, row.names = F, col.names = F)

  # Обращаемся к питону из виртуального окружения, куда установлен pandas и biopython
  user_name <- system('echo "$USER"', intern = T)
  add_inf_cmd  <- paste0("/home/", user_name, "/user_venv/bin/python ./inst/additional_inf.py -f ", output_filename_base, "_", type, ".lst -o ", output_filename_base, "_", type, ".tsv")

  # Запускаем скрипт, подгружающий дополнительную информацию о публикации
  system(add_inf_cmd)

  # Загружаем полученную информацию и мерджим как дополнительные поля к нашей таблице
  data_table_add_inf <- read.csv(paste0(output_filename_base, "_", type, ".tsv"), sep = '\t')
  data_table_with_add_inf <- merge(data_table, data_table_add_inf)

  # Чистим временные файлы
  rm_cmd1 <- paste0("rm ", output_filename_base, "_", type, ".tsv")
  rm_cmd2 <- paste0("rm ", output_filename_base, "_", type, ".lst")

  system(rm_cmd1)
  system(rm_cmd2)

  return(data_table_with_add_inf)
}
