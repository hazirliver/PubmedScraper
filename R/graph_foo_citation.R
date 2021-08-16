#' @title Citation_by main function
#' @description Creates citation_by graph and corresponding data.frame with
#' additional info about the number of citations of each paper.

#' @importFrom igraph graph.data.frame degree delete.vertices plot.igraph V `V<-` layout_nicely
#' @import dplyr
#' @importFrom stringr str_detect
#' @importFrom grDevices svg dev.off
#'
#' @param PMID_list A character vector of PMIDs.
#' @param output_filename_base A character with path to output folder plus filename base.
#' @param n_cites_threshold An integer. The minimum number of papers that a given article should refer to.
#'
#' @return A data.frame with info about papers that refer to a given \code{PMID_list}.
#' @export
graph_foo_citation <- function(PMID_list, output_filename_base, n_cites_threshold = 2)
{
  # Качаем информацию о том, какие статьи цитируют заданный пул
  data_table <- cited_by(PMID_list)

  # Строим граф цитирований.
  # Считаем степень вершины, как количество статей, которые цитирует заданная
  # Количество стрелочек ВЫходящих из вершины
  df.g0 <- graph.data.frame(d = data_table, directed = T)
  deg0 <- degree(df.g0, mode="out")

  # Для построения более красивого графика убираем вершины:
  # 1) Степень вершины = 1 и при этом это не вершина из изначального пула
  # 2) Исходящая степень вершины = 0 и при этом входящая степень = 0
  df.g <- delete.vertices(df.g0,
                          (!str_detect(names(deg0),paste0(PMID_list, collapse = "|")) & deg0 == 1) |
                            (deg0 == 0 & degree(df.g0, mode="in") == 0))
  deg <- degree(df.g, mode = "out")


  ### Построение рисунка графа
  svg(paste0(output_filename_base, "_cited_by.svg"), width = 20, height = 25)

  # Прозрачная заливка
  # Красная граница -- изначальный пул, черная -- новые статьи
  V(df.g)$color <- NA
  bc <- ifelse(str_detect(names(deg),paste0(PMID_list, collapse = "|")), "red", "black")



  plot.igraph(df.g, vertex.frame.color=bc, layout=layout_nicely,
              vertex.label = NA, vertex.size=deg, edge.arrow.size = 0.5)

  dev.off()
  ###


  # Для дальнейшей работы используется вся информация (без удаления вершин)
  # Таблица: PMID, степень вершины
  # Показывает, на сколько статей из указанного пула ссылается заданная статья
  deg_to_csv <- data.frame(PMID=names(deg0), n_cites=deg0, row.names=NULL) %>%
    arrange(desc(n_cites))


  # Фильруем по минимальному количеству статей, на которые должна ссылаться заданная статья
  # Получаем список отобранных PMID для дальнейшего прогона
  double_info_PMID_list <- deg_to_csv %>%
    filter(n_cites >= n_cites_threshold) %>%
    select(PMID) %>%
    pull()

  # Получаем аналогичную data_table табличку для отобранных ранее PMID
  double_info_data_table <- cited_by(double_info_PMID_list)


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
    arrange(desc(n_cites))

  top_papers_double_to_csv$is_in_initial_pull <- ifelse(str_detect(top_papers_double_to_csv$PMID,
                                                                   paste0(PMID_list, collapse = "|")),
                                                        "initial", "new")

  return(top_papers_double_to_csv)
}
