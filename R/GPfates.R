#' Run GPfates
#'
#' @param counts The non-normalised expression counts
#' @param nfates The number of fates
#' @param ndims The number of dimensions
#' @param log_expression_cutoff The log expression cutoff
#' @param min_cells_expression_cutoff The min expression cutoff
#'
#' @export
#' @importFrom glue glue
#' @importFrom readr read_csv
#' @importFrom utils write.table
GPfates <- function(
  counts,
  nfates = 1,
  ndims = 2,
  log_expression_cutoff = 2,
  min_cells_expression_cutoff = 2
) {
  temp_folder <- tempdir()

  utils::write.table(t(counts), paste0(temp_folder, "expression.csv"), sep="\t")

  utils::write.table(
    data.frame(cell_id = rownames(counts), row.names = rownames(counts)),
    paste0(temp_folder, "cellinfo.csv"),
    sep="\t"
  )

  system2(
    "/bin/bash",
    args = c(
      "-c",
      shQuote(glue::glue(
        "cd {find.package('GPfates')}/venv",
        "source bin/activate",
        "python3 {find.package('GPfates')}/wrapper.py {temp_folder} {log_expression_cutoff} {min_cells_expression_cutoff} {nfates} {ndims}",
        .sep = ";"))
    )
  )

  pseudotime <- readr::read_csv(glue::glue("{temp_folder}pseudotimes.csv"), col_names = c("cell_id", "time"))

  phi <- readr::read_csv(glue::glue("{temp_folder}phi.csv"), col_names = c("cell_id", glue::glue("M{seq_len(nfates)}")), skip = 1)

  dr <- readr::read_csv(glue::glue("{temp_folder}dr.csv"), col_names = c("cell_id", glue::glue("Comp{seq_len(5)}")), skip = 1)

  # remove temporary output
  unlink(temp_folder, recursive = TRUE)

  list(
    pseudotime = pseudotime,
    phi = phi,
    dr = dr
  )
}
