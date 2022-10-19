#' @title Convert one or more data sets to taxmap
#'
#' @description This function is a modification of metacoder's
#'   \code{lookup_tax_data}, which looks up taxonomic data from NCBI sequence
#'   IDs, taxon IDs, or taxon names that are present in a table, list, or
#'   vector. Also can incorporate additional associated datasets. This edit
#'   automatically resubmits the query to NCBI if the request times out.
#'
#' @param tax_data A table, list, or vector that contain sequence IDs, taxon
#'   IDs, or taxon names. * tables: The `column` option must be used to specify
#'   which column contains the sequence IDs, taxon IDs, or taxon names. * lists:
#'   There must be only one item per list entry unless the `column` option is
#'   used to specify what item to use in each list entry. * vectors: simply a
#'   vector of sequence IDs, taxon IDs, or taxon names.
#' @param type What type of information can be used to look up the
#'   classifications. Takes one of the following values: * `"seq_id"`: A
#'   database sequence ID with an associated classification (e.g. NCBI accession
#'   numbers). * `"taxon_id"`: A reference database taxon ID (e.g. a NCBI taxon
#'   ID) * `"taxon_name"`: A single taxon name (e.g. "Homo sapiens" or
#'   "Primates") * `"fuzzy_name"`: A single taxon name, but check for
#'   misspellings first. Only use if you think there are misspellings. Using
#'   `"taxon_name"` is faster.
#' @param column (`character` or `integer`) The name or index of the column that
#'   contains information used to lookup classifications. This only applies when
#'   a table or list is supplied to `tax_data`.
#' @param datasets Additional lists/vectors/tables that should be included in
#'   the resulting `taxmap` object. The `mappings` option is use to specify how
#'   these data sets relate to the `tax_data` and, by inference, what taxa apply
#'   to each item.
#' @param mappings (named `character`) This defines how the taxonomic
#'   information in `tax_data` applies to data in `datasets`. This option should
#'   have the same number of inputs as `datasets`, with values corresponding to
#'   each dataset. The names of the character vector specify what information in
#'   `tax_data` is shared with info in each `dataset`, which is specified by the
#'   corresponding values of the character vector. If there are no shared
#'   variables, you can add `NA` as a placeholder, but you could just leave that
#'   data out since it is not benefiting from being in the taxmap object. The
#'   names/values can be one of the following: * For tables, the names of
#'   columns can be used. * `"{{index}}"` : This means to use the index of
#'   rows/items * `"{{name}}"`  : This means to use row/item names. *
#'   `"{{value}}"` : This means to use the values in vectors or lists. Lists
#'   will be converted to vectors using [unlist()].
#' @param database (`character`) The name of a database to use to look up
#'   classifications. Options include "ncbi", "itis", "eol", "col", "tropicos",
#'   and "nbn".
#' @param include_tax_data (`TRUE`/`FALSE`) Whether or not to include `tax_data`
#'   as a dataset, like those in `datasets`.
#' @param use_database_ids (`TRUE`/`FALSE`) Whether or not to use downloaded
#'   database taxon ids instead of arbitrary, automatically-generated taxon ids.
#' @param ask  (`TRUE`/`FALSE`) Whether or not to prompt the user for input.
#'   Currently, this would only happen when looking up the taxonomy of a taxon
#'   name with multiple matches. If `FALSE`, taxa with multiple hits are treated
#'   as if they do not exist in the database. This might change in the future if
#'   we can find an elegant way of handling this.
#'
#' @section Failed Downloads: If you have invalid inputs or a download fails for
#'   another reason, then there will be a "unknown" taxon ID as a placeholder
#'   and failed inputs will be assigned to this ID. You can remove these using
#'   [filter_taxa()] like so: `filter_taxa(result, taxon_ids != "unknown")`. Add
#'   `drop_obs = FALSE` if you want the input data, but want to remove the
#'   taxon.
#'
#' @family parsers
#'
#' @examples \dontrun{
#'
#'   # Look up taxon names in vector from NCBI
#'   lookup_tax_data(c("homo sapiens", "felis catus", "Solanaceae"),
#'                   type = "taxon_name")
#'
#'   # Look up taxon names in list from NCBI
#'   lookup_tax_data(list("homo sapiens", "felis catus", "Solanaceae"),
#'                   type = "taxon_name")
#'
#'   # Look up taxon names in table from NCBI
#'   my_table <- data.frame(name = c("homo sapiens", "felis catus"),
#'                          decency = c("meh", "good"))
#'   lookup_tax_data(my_table, type = "taxon_name", column = "name")
#'
#'   # Look up taxon names from NCBI with fuzzy matching
#'   lookup_tax_data(c("homo sapienss", "feles catus", "Solanacese"),
#'                   type = "fuzzy_name")
#'
#'   # Look up taxon names from a different database
#'   lookup_tax_data(c("homo sapiens", "felis catus", "Solanaceae"),
#'                   type = "taxon_name", database = "ITIS")
#'
#'   # Prevent asking questions for ambiguous taxon names
#'   lookup_tax_data(c("homo sapiens", "felis catus", "Solanaceae"),
#'                   type = "taxon_name", database = "ITIS", ask = FALSE)
#'
#'   # Look up taxon IDs from NCBI
#'   lookup_tax_data(c("9689", "9694", "9643"), type = "taxon_id")
#'
#'   # Look up sequence IDs from NCBI
#'   lookup_tax_data(c("AB548412", "FJ358423", "DQ334818"),
#'                   type = "seq_id")
#'
#'   # Make up new taxon IDs instead of using the downloaded ones
#'   lookup_tax_data(c("AB548412", "FJ358423", "DQ334818"),
#'                   type = "seq_id", use_database_ids = FALSE)
#'
#'
#'   # --- Parsing multiple datasets at once (advanced) ---
#'   # The rest is one example for how to classify multiple datasets at once.
#'
#'   # Make example data with taxonomic classifications
#'   species_data <- data.frame(tax = c("Mammalia;Carnivora;Felidae",
#'                                      "Mammalia;Carnivora;Felidae",
#'                                      "Mammalia;Carnivora;Ursidae"),
#'                              species = c("Panthera leo",
#'                                          "Panthera tigris",
#'                                          "Ursus americanus"),
#'                              species_id = c("A", "B", "C"))
#'
#'   # Make example data associated with the taxonomic data
#'   # Note how this does not contain classifications, but
#'   # does have a varaible in common with "species_data" ("id" = "species_id")
#'   abundance <- data.frame(id = c("A", "B", "C", "A", "B", "C"),
#'                           sample_id = c(1, 1, 1, 2, 2, 2),
#'                           counts = c(23, 4, 3, 34, 5, 13))
#'
#'   # Make another related data set named by species id
#'   common_names <- c(A = "Lion", B = "Tiger", C = "Bear", "Oh my!")
#'
#'   # Make another related data set with no names
#'   foods <- list(c("ungulates", "boar"),
#'                 c("ungulates", "boar"),
#'                 c("salmon", "fruit", "nuts"))
#'
#'   # Make a taxmap object with these three datasets
#'   x = lookup_tax_data(species_data,
#'                       type = "taxon_name",
#'                       datasets = list(counts = abundance,
#'                                       my_names = common_names,
#'                                       foods = foods),
#'                       mappings = c("species_id" = "id",
#'                                    "species_id" = "{{name}}",
#'                                    "{{index}}" = "{{index}}"),
#'                       column = "species")
#'
#'   # Note how all the datasets have taxon ids now
#'   x$data
#'
#'   # This allows for complex mappings between variables that other functions use
#'   map_data(x, my_names, foods)
#'   map_data(x, counts, my_names)
#' }
#' @export

lookup_tax_data <- function(tax_data, type, column = 1, datasets = list(),
                            mappings = c(), database = "ncbi",
                            include_tax_data = TRUE, use_database_ids = TRUE,
                            ask = TRUE) {

  # Check for zero-length inputs
  if (length_of_thing(tax_data) <= 0) {
    return(metacoder::taxmap(data = list(tax_data = tax_data)))
  }

  # Check that a supported database is being used
  supported_databases <- names(metacoder::database_list)
  database <- tolower(database)
  if (! database %in% supported_databases) {
    stop(paste0('The database "', database,
                '" is not a valid database for looking up that taxonomy of ',
                'sequnece ids. Valid choices include:\n',
                limited_print(supported_databases, type = "silent")))
  }

  # Check that column exists
  check_class_col(tax_data, column)

  # Hidden parameters
  batch_size <- 100
  max_print <- 10
  internal_class_sep <- "||||"
  internal_class_name <- "___class_string___"

  # Define lookup functions
  format_class_table <- function(class_table) {
    # Complain about failed queries
    failed_queries <- is.na(class_table)
    if (sum(failed_queries) > 0) {
      failed_names <- unique(names(failed_queries[failed_queries]))
      error_msg <- paste0('The following ', length(failed_names),
                          ' unique inputs could not be looked up:\n  ',
                          limited_print(failed_names,
                                                   type = "silent"))
      if (ask) {
        error_msg <- paste0(error_msg, "\n",
                            'This probably means they dont exist in the database "', database, '".')
      } else {
        error_msg <- paste0(error_msg, "\n",
                            'This probably means they dont exist in the database "', database,
                            '" or have multiple matches. ',
                            'Use "ask = TRUE" to specify which is the correct match when multiple matches occur.')
      }
      warning(error_msg, call. = FALSE)
    }

    # Rename columns of result
    col_names <- c(paste0(database, "_name"),
                   paste0(database, "_rank"),
                   paste0(database, "_id"))
    class_table[! failed_queries] <- lapply(class_table[! failed_queries],
                                            function(x) stats::setNames(x, col_names))

    # Add placeholders for failed requests
    class_table[failed_queries] <- lapply(seq_len(sum(failed_queries)),
                                          function(i) {
                                            out <- dplyr::tibble(a = "unknown taxon",
                                                                 b = "unknown rank",
                                                                 c = "unknown")
                                            colnames(out) <- col_names
                                            return(out)
                                          })
    return(class_table)
  }


  use_taxon_id <- function(ids) {
    message("Looking up classifications for ", length(unique(ids)),
            ' unique taxon IDs from database "', database, '"...')

    # Look up classifications
    lookup_all <- function(ids) {
      lookup_one <- function(id) {
        taxize::classification(id, ask = FALSE, rows = 1, db = database,
                               message = FALSE)
      }
      output <- progress_lapply(ids, lookup_one)
      return(output)
    }
    tryCatch(msgs <- utils::capture.output(raw_result <- map_unique(ids, lookup_all),
                                           type = "message"),
             error = function(e) stop(e))


    # Remove repeated messages (e.g. no NCBI API key)
    on.exit(message(paste0(unique(msgs), collapse = "\n")))

    # Reformat result
    result <- stats::setNames(unlist(raw_result, recursive = FALSE), ids)
    format_class_table(result)
  }

  use_seq_id <- function(ids) {
    # Check that a supported database is being used
    supported_databases <- c("ncbi")
    if (! database %in% supported_databases) {
      stop(call. = FALSE,
           paste0('The database "', database,
                  '" is not a valid database for looking up that taxonomy of ',
                  'sequnece ids. Valid choices include:\n',
                  limited_print(supported_databases, type = "silent")))
    }

    # Look up classifications
    message("Looking up classifications for ", length(unique(ids)), " unique sequence IDs from NCBI...")
    lookup_all <- function(ids) {
      lookup_one <- function(id) {
        # Re-attempt query if HTTP request timeout error
        attempt <- 0
        result <- NULL
        while(is.null(result) && attempt <= 5){
          attempt <- attempt + 1
          try(result <- taxize::classification(taxize::genbank2uid(id)[1], db = database))

        }
        result
      }



      output <- progress_lapply(ids, lookup_one)
      return(output)
    }
    tryCatch(msgs <- utils::capture.output(raw_result <- map_unique(ids, lookup_all),
                                           type = "message"),
             error = function(e) stop(e))

    # Remove repeated messages (e.g. no NCBI API key)
    on.exit(message(paste0(unique(msgs), collapse = "\n")))

    # Reformat result
    result <- stats::setNames(unlist(raw_result, recursive = FALSE), ids)
    format_class_table(result)
  }

  use_taxon_name <- function(my_names) {
    message("Looking up classifications for ", length(unique(my_names)),
            ' unique taxon names from database "', database, '"...')

    # Look up classifications
    lookup_all <- function(my_names) {
      lookup_one <- function(my_name) {
        taxize::classification(my_name, ask = ask, db = database,
                               messages = FALSE)
      }
      output <- progress_lapply(my_names, lookup_one)
      return(output)
    }
    tryCatch(msgs <- utils::capture.output(raw_result <- map_unique(my_names, lookup_all),
                                           type = "message"),
             error = function(e) stop(e))

    # Remove repeated messages (e.g. no NCBI API key)
    on.exit(message(paste0(unique(msgs), collapse = "\n")))

    # Reformat result
    result <- stats::setNames(unlist(raw_result, recursive = FALSE), my_names)
    format_class_table(result)
  }

  use_taxon_name_fuzzy <- function(my_names) {
    message("Looking up classifications for ", length(unique(my_names)),
            ' unique taxon names from database "', database, '" using fuzzy name matching...')

    # Look up similar taxon names
    corrected <- map_unique(my_names, correct_taxon_names, database = database)

    # Check for not found names
    not_found <- unique(names(corrected[is.na(corrected)]))
    if (length(not_found) > 0) {
      warning(call. = FALSE,
              "No taxon name was found that was similar to the following ",
              length(not_found), " inputs:\n  ",
              limited_print(type = "silent", not_found))

    }

    # Replace not found values with original input
    corrected[is.na(corrected)] <- names(corrected[is.na(corrected)])

    # Check for changed names
    changed <- tolower(names(corrected)) != tolower(corrected) & ! is.na(corrected) & ! is.na(names(corrected))
    if (any(changed)) {
      before <- names(corrected[changed])
      after <- corrected[changed]
      change_char <- unique(paste0('"', before, '" -> "', after, '"'))
      message("The following ", length(change_char), " names were corrected before looking up classifications:\n  ",
              limited_print(type = "silent", change_char))
    }

    # Run standard taxon name lookup
    use_taxon_name(corrected)
  }

  lookup_funcs <- list("seq_id" = use_seq_id,
                       "taxon_id" = use_taxon_id,
                       "taxon_name" = use_taxon_name,
                       "fuzzy_name" = use_taxon_name_fuzzy)

  # Get query information
  if (is.data.frame(tax_data)) { # is table
    query <- as.character(tax_data[[column]])
  } else if (is.list(tax_data)) { # is list
    query <- vapply(tax_data,
                    function(x) as.character(x[[column]]),
                    character(1))
  } else if (is.vector(tax_data)) { # is vector
    query <- as.character(tax_data)
  }

  # Look up taxonomic classifications
  if (! type %in% names(lookup_funcs)) {
    stop(paste0('Invalid "type" option. It must be one of the following:\n  ',
                paste0(names(lookup_funcs), collapse = ", ")))
  }
  classifications <- lookup_funcs[[type]](query)
  class_strings <- unlist(lapply(classifications, function(x) {
    lapply(seq_len(nrow(x)), function(i) {
      paste0(as.character(unlist(x[1:i, 1])), "___", as.character(unlist(x[1:i, 2])), collapse = internal_class_sep)
    })
  }))
  combined_class <- do.call(rbind, unname(classifications))
  internal_class_frame <- stats::setNames(data.frame(class_strings,
                                                     stringsAsFactors = FALSE),
                                          internal_class_name)
  combined_class <- cbind(internal_class_frame,
                          combined_class,
                          stringsAsFactors = FALSE)

  # Add mapping columns to classfication data
  tax_data_indexes <- cumsum(vapply(classifications, nrow, numeric(1)))
  mappping_cols <- unique(c(names(mappings), "{{index}}", "{{name}}"))
  if (!is.data.frame(tax_data)) {
    mappping_cols <- c(mappping_cols, "{{value}}")
  }
  for (col in mappping_cols) {
    combined_class[tax_data_indexes, col] <- get_sort_var(tax_data, col)
  }

  # Add input data to datasets included in the resulting object
  if (include_tax_data) {
    datasets <- c(list(query_data = tax_data), datasets)
    mappings <- c("{{index}}" = "{{index}}", mappings)
  }

  # Make taxmap object
  output <- metacoder::parse_tax_data(tax_data = combined_class,
                                      datasets = datasets,
                                      class_cols = 1,
                                      class_sep = internal_class_sep,
                                      class_key = c('taxon_name', 'taxon_rank'),
                                      class_regex = '^(.+)___(.+)$',
                                      mappings = mappings,
                                      include_tax_data = include_tax_data)
  output$data$class_data <- NULL

  if (include_tax_data) {
    # Remove mapping columns from output
    output$data$tax_data[mappping_cols] <- NULL

    # Remove class column from output
    output$data$tax_data[internal_class_name] <- NULL

    # Fix incorrect taxon ids for tax_data
    #   This is due to the "{{index}}" being interpreted as a column name,
    #   which is needed for the user-defined data sets to be parsed right.
    output$data$tax_data$taxon_id <- output$input_ids

    # Remove duplicates from `tax_data`
    output$data$tax_data <- output$data$tax_data[!duplicated(output$data$tax_data), ]
  }

  # Replace standard taxon ids with database taxon ids
  if (use_database_ids) {
    taxon_id_col <- paste0(database, "_id")
    # I am not sure why the following line works...
    new_ids <- unique(combined_class[[taxon_id_col]])[match(output$taxon_ids(),
                                                            unique(output$input_ids))]
    output$replace_taxon_ids(new_ids)
  }

  return(output)
}
