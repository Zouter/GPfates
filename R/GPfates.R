#' Run GPfates
#'
#' @param counts The non-normalised expression counts
#' @param nfates The number of fates
#' @param ndims The number of dimensions
#' @param log_expression_cutoff The log expression cutoff
#' @param min_cells_expression_cutoff The min expression cutoff
#' @param verbose Print output produced by method
#' @param num_cores The number of cores GPfates is allowed to use
#'
#' @export
#' @importFrom glue glue
#' @importFrom readr read_csv
#' @importFrom utils write.table
#' @importFrom methods formalArgs
#' @importFrom jsonlite toJSON
#' @importFrom dynutils run_until_exit
GPfates <- function(
  counts,
  nfates = 1,
  ndims = 2,
  log_expression_cutoff = 2,
  min_cells_expression_cutoff = 2,
  num_cores = 1,
  verbose = FALSE
) {
  # create temporary folder
  temp_folder <- tempfile()
  dir.create(temp_folder, recursive = TRUE)

  counts_file <- paste0(temp_folder, "/expression.tsv")
  cellinfo_file <- paste0(temp_folder, "/cellinfo.tsv")
  json_file <- paste0(temp_folder, "/params.json")

  tryCatch({
    # write counts to temporary file
    utils::write.table(t(counts), counts_file, sep="\t")

    # write cell info to temporary file
    cell_info <- data.frame(cell_id = rownames(counts), row.names = rownames(counts))
    utils::write.table(cell_info, cellinfo_file, sep="\t")

    # write parameters to temporary folder
    params <- as.list(environment())[methods::formalArgs(GPfates)]
    params <- params[names(params) != "counts"]
    write(jsonlite::toJSON(params, auto_unbox = TRUE), json_file)

    # execute python script
    if (!is.null(num_cores)) {
      num_cores_str <- glue::glue(
        "export MKL_NUM_THREADS={num_cores};",
        "export NUMEXPR_NUM_THREADS={num_cores};",
        "export OMP_NUM_THREADS={num_cores}"
      )
    } else {
      num_cores_str <- "echo 'no cores'"
    }

    commands <- glue::glue(
      "cd {find.package('GPfates')}/venv",
      "source bin/activate",
      "{num_cores_str}",
      "python3 {find.package('GPfates')}/wrapper.py {temp_folder}",
      .sep = ";"
    )

    output <- processx::run("/bin/bash", c("-c", commands), echo=TRUE)

    if (verbose) cat(output$output, "\n", sep="")

    # read output
    pseudotime <- readr::read_csv(
      glue::glue("{temp_folder}/pseudotimes.csv"),
      col_types = readr::cols(cell_id = "c", time = "d"),
      col_names = c("cell_id", "time")
    )
    phi <- readr::read_csv(
      glue::glue("{temp_folder}/phi.csv"),
      col_names = c("cell_id", paste0("M", seq_len(nfates))),
      col_types = readr::cols(.default = "d", cell_id = "c"),
      skip = 1
    )
    dr <- readr::read_csv(
      glue::glue("{temp_folder}/dr.csv"),
      col_names = c("cell_id", paste0("Comp", seq_len(5))),
      col_types = readr::cols(.default = "d", cell_id = "c"),
      skip = 1
    )

  }, finally = {
    # remove temporary output
    unlink(temp_folder, recursive = TRUE)
  })

  list(
    pseudotime = pseudotime,
    phi = phi,
    dr = dr
  )
}
