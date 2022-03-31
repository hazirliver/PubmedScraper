#' @title References main function
#' @description Creates references graph and corresponding data.frame with
#' additional info about the number of citations of each paper.

#' @importFrom igraph graph.data.frame degree delete.vertices plot.igraph V `V<-` layout_nicely
#' @import dplyr
#' @importFrom stringr str_detect
#' @importFrom grDevices svg dev.off
#'
#' @param PMID_list A character vector of PMIDs.
#' @param output_filename_base A character with path to output folder plus filename base.
#' @param n_refers_treshold An integer. The minimum number of papers that a given article should refer to.
#'
#' @return A data.frame with info about papers that refer to a given \code{PMID_list}.
#' @export
graph_foo_references <- function(PMID_list, output_filename_base, n_refers_treshold = 2)
{
  # Качаем информацию о том, на какие статьи ссылается каждая статья из заданного пула
  references_table <- references(PMID_list)

  # Строим граф цитирований.
  # Считаем степень вершины, как количество статей, которые процитировали заданную статью
  # Количество стрелочек Входящих в вершину
  df.g_references <- graph.data.frame(d = references_table, directed = T)
  deg_references <- degree(df.g_references, mode="out")


  # ### Построение рисунка графа
  # svg(paste0(output_filename_base, "_references.svg"), width = 20, height = 25)
  #
  # # Прозрачная заливка
  # # Красная граница -- изначальный пул, черная -- новые статьи
  # V(df.g_references)$color <- NA
  # bc <- ifelse(str_detect(names(deg_references),paste0(PMID_list, collapse = "|")), "red", "black")
  #
  #
  #
  # plot.igraph(df.g_references, vertex.frame.color=bc, layout=layout_nicely,
  #             vertex.label = NA, vertex.size=deg_references, edge.arrow.size = 0.5)
  #
  # dev.off()
  # ###


  # Для дальнейшей работы используется вся информация (без удаления вершин)
  # Таблица: PMID, степень вершины
  # Показывает, сколько статей из указанного пула ссылаются на заданную статью
  deg_to_csv <- data.frame(PMID=names(deg_references), n_refers=deg_references, row.names=NULL) %>%
    arrange(desc(n_refers))


  # Фильруем по минимальному количеству статей, на которые должна ссылаться заданная статья
  # Получаем список отобранных PMID для дальнейшего прогона
  double_info_pmid_list <- deg_to_csv %>%
    dplyr::filter(n_refers > n_refers_treshold) %>%
    pull(PMID)

  # Получаем аналогичную data_table табличку для отобранных ранее PMID
  double_info_data_table <- cited_by(double_info_pmid_list)

  # Группируем и сортируем статьи по цитируемости
  top_papers_double <- double_info_data_table  %>%
    group_by(parent_work) %>%
    summarise(n_cit = n()) %>%
    arrange(desc(n_cit))

  colnames(top_papers_double) <- c("PMID", "n_cited")

  # Соединияем две полученные таблицы
  # Если Какую-то статью из отобранных на первом шаге никто не процитировал, у нее значение 0
  top_papers_double_to_csv <- merge(top_papers_double, deg_to_csv, by = "PMID")  %>%
    replace(is.na(.), 0)  %>%
    arrange(desc(n_refers))

  top_papers_double_to_csv$is_in_initial_pull <- ifelse(str_detect(top_papers_double_to_csv$PMID,
                                                                   paste0(PMID_list, collapse = "|")),
                                                        "initial", "new")

  return(top_papers_double_to_csv)
}
