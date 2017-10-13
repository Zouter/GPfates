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

  # write counts to temporary file
  utils::write.table(t(counts), paste0(temp_folder, "/expression.tsv"), sep="\t")

  # write cell info to temporary file
  cell_info <- data.frame(cell_id = rownames(counts), row.names = rownames(counts))
  utils::write.table(cell_info, paste0(temp_folder, "/cellinfo.tsv"), sep="\t")

  # write parameters to temporary folder
  params <- as.list(environment())[methods::formalArgs(GPfates)]
  params <- params[names(params) != "counts"]
  write(jsonlite::toJSON(params, auto_unbox = TRUE), paste0(temp_folder, "/params.json"))

  # execute python script
  cores_string <-
    if (is.null(num_cores)) {
      ""
    } else {
      glue::glue(
        "export MKL_NUM_THREADS={num_cores};",
        "export NUMEXPR_NUM_THREADS={num_cores};",
        "export OMP_NUM_THREADS={num_cores};"
      )
    }
  output <- system2(
    "/bin/bash",
    args = c(
      "-c",
      shQuote(glue::glue(
        "cd {find.package('GPfates')}/venv;",
        "source bin/activate;",
        "{cores_string}",
        "python3 {find.package('GPfates')}/wrapper.py {temp_folder}"
      ))
    )
  )

  if (verbose) cat(output, "\n", sep="")

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

  # remove temporary output
  unlink(temp_folder, recursive = TRUE)

  list(
    pseudotime = pseudotime,
    phi = phi,
    dr = dr
  )
}
