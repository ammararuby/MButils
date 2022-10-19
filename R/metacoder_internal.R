#' Print a subset of a character vector
#'
#' Prints the start and end values for a character vector. The number of values
#' printed depend on the width of the screen by default.
#'
#' @param chars (`character`) What to print.
#' @param prefix (`character` of length 1) What to print before `chars`, on the
#'   same line.
#' @param sep What to put between consecutive values
#' @param mid What is used to indicate omitted values
#' @param trunc What is appended onto truncated values
#' @param max_chars (`numeric` of length 1) The maximum number of characters to
#'   print.
#' @param type (`"error"`, `"warning"`, `"message"`, `"cat"`, `"print"`,
#'   `"silent"`, `"plain"`)
#' @importFrom crayon strip_style
#'
#' @return `NULL`
#' @references https://github.com/grunwaldlab/metacoder
#' @keywords internal
#' @noRd
limited_print <- function(chars, prefix = "", sep = ", ", mid = " ... ",
                          trunc_char = "[truncated]",
                          max_chars = getOption("width") - nchar(prefix) - 5,
                          type = "message") {

  # https://stat.ethz.ch/pipermail/r-help/2006-March/101023.html
  interleave <- function(v1,v2) {
    ord1 <- 2*(1:length(v1))-1
    ord2 <- 2*(1:length(v2))
    c(v1,v2)[order(c(ord1,ord2))]
  }

  truncate <- function(x, max_chars = 30) {
    if (nchar(x) > max_chars) {
      x <- paste0(substr(x, 0, max_chars - nchar(crayon::strip_style(trunc_char))), trunc_char)
    }
    return(x)
  }

  # Remove colsole fonts
  raw_chars <- chars
  chars <- crayon::strip_style(chars)

  # Convert NA to "NA"
  chars[is.na(chars)] <- "NA"

  #
  if (length(chars) == 0) {
    output <- prefix
    return(invisible(NULL))
  }

  #
  q = "'"
  interleaved <- interleave(chars[1:(length(chars) / 2)],
                            rev(chars[(length(chars) / 2 + 1):length(chars)]))
  is_greater_than_max <- cumsum(nchar(interleaved) + nchar(crayon::strip_style(sep))) + 10 > max_chars
  if (all(! is_greater_than_max)) {
    max_printed <- length(chars)
  } else {
    max_printed <- which.max(is_greater_than_max) - 1
  }
  if (max_printed < length(chars)) {
    if (max_printed < 2) {
      first_part <- truncate(chars[1])
      second_part <- truncate(chars[length(chars)])
    } else {
      first_part <-  raw_chars[1:ceiling(max_printed / 2)]
      second_part <- raw_chars[(length(chars) - floor(max_printed / 2) + 1):length(chars)]
    }
    if (length(chars) > 1) {
      output <- paste0(paste0(collapse = sep, first_part),
                       mid,
                       paste0(collapse = sep, second_part),
                       "\n")
    } else {
      output <- paste0(paste0(collapse = sep, first_part),
                       "\n")

    }
  } else {
    output <- paste0(paste0(collapse = sep, raw_chars), "\n")
  }
  output <- paste(prefix, output, collapse = "")

  if (type == "error") {
    stop(output)
  } else if (type == "warning") {
    warning(output)
  } else if (type == "message") {
    message(output)
  } else if (type == "cat") {
    cat(output)
  } else if (type == "print") {
    print(output)
  } else if (type == "plain") {
    output <- crayon::strip_style(output)
  } else if (type != "silent") {
    stop("invalid type option")
  }
  return(invisible(output))
}

#' Check for name/index in input data
#'
#' Used by `lookup_tax_data` to check that columm/class_col is valid for the
#' input data
#'
#' @usage check_class_col(tax_data, column)
#'
#' @param tax_data A table, list, or vector that contain sequence IDs, taxon
#'   IDs, or taxon names. * tables: The 'column' option must be used to specify
#'   which column contains the sequence IDs, taxon IDs, or taxon names. * lists:
#'   There must be only one item per list entry unless the 'column' option is
#'   used to specify what item to use in each list entry. * vectors: simply a
#'   vector of sequence IDs, taxon IDs, or taxon names.
#' @param column ('character' or 'integer') The name or index of the column that
#'   contains information used to lookup classifications. This only applies when
#'   a table or list is supplied to 'tax_data'.
#' @references https://github.com/grunwaldlab/metacoder
#' @keywords internal
#' @noRd
check_class_col <- function (tax_data, column){
  if (is.data.frame(tax_data)) {
    if (is.numeric(column)) {
      if (column == 0 || abs(column) > ncol(tax_data)) {
        stop(call. = FALSE, "Column index \"", column,
             "\" out of bounds. Must be between 1 and ",
             ncol(tax_data), ".")
      }
    }
    else if (!column %in% colnames(tax_data)) {
      stop(call. = FALSE, "No column \"", column, "\" in input table. Valid columns include:\n  ",
           limited_print(colnames(tax_data), type = "silent"))
    }
  }
  else if (is.list(tax_data) || is.vector(tax_data)) {
    my_lengths <- vapply(tax_data, length, numeric(1))
    had_col_name <- vapply(tax_data, function(x) column %in%
                             names(x), logical(1))
    if (is.numeric(column)) {
      if (column < 1 || any(column > my_lengths)) {
        stop(call. = FALSE, "Column index \"", column,
             "\" out of bounds for inputs:\n", limited_print(which(column >
                                                                     my_lengths), type = "silent"))
      }
    }
    else if (!all(had_col_name)) {
      stop(call. = FALSE, "No item named \"", column, "\" in the following inputs:\n",
           limited_print(which(!had_col_name), type = "silent"))
    }
  }
  else {
    stop(call. = FALSE, "Cannot read input of class \"",
         class(tax_data)[1], "\". Input must be a table, list or vector.")
  }
}

#' Get indexes of a unique set of the input
#'
#' @description Get indexes of a unique set of the input
#'
#' @usage unique_mapping(input)
#'
#' @references https://github.com/grunwaldlab/metacoder
#' @keywords internal
#' @noRd
unique_mapping <- function (input) {
  unique_input <- unique(input)
  vapply(input, function(x) {
    if (is.na(x))
      which(is.na(unique_input))
    else which(x == unique_input)
  }, numeric(1))
}

#' Run a function on unique values of a iterable
#'
#' Runs a function on unique values of a list/vector and then reformats the
#' output so there is a one-to-one relationship with the input.
#'
#' @param input What to pass to \code{func}
#' @param func (\code{function})
#' @param ... passed to \code{func}
#' @references https://github.com/grunwaldlab/metacoder
#' @keywords internal
#' @noRd
map_unique <- function(input, func, ...) {
  input_class <- class(input)
  unique_input <- unique(input)
  class(unique_input) <- input_class
  func(unique_input, ...)[unique_mapping(input)]
}

#' Get a vector from a vector/list/table to be used in mapping
#'
#' Get a vector from a vector/list/table to be used in mapping
#'
#' @usage get_sort_var(data, var)
#'
#' @param data A vector/list/table
#' @param var What to get. * For tables, the names of columns can be used. *
#'   '"index"' : This means to use the index of rows/items * '"name"' : This
#'   means to use row/item names. * '"value"' : This means to use the values in
#'   vectors or lists.
#' @keywords internal
#' @references https://github.com/grunwaldlab/metacoder
#' @noRd
get_sort_var <- function (data, var) {
  if (is.data.frame(data) && var %in% colnames(data)) {
    return(data[[var]])
  }
  else if (var == "{{index}}") {
    if (is.data.frame(data)) {
      return(seq_len(nrow(data)))
    }
    else {
      return(seq_len(length(data)))
    }
  }
  else if (var == "{{name}}") {
    if (is.data.frame(data)) {
      return(rownames(data))
    }
    else {
      if (is.null(names(data))) {
        return(rep(NA_character_, length(data)))
      }
      else {
        return(names(data))
      }
    }
  }
  else if (var == "{{value}}") {
    if (is.data.frame(data)) {
      stop("The `{{value}}` setting of the `mappings` option cannot be used with data.frames.")
    }
    else {
      return(unlist(data))
    }
  }
  else {
    stop(paste0("No column named \"", var, "\".\""))
  }
}

#' Check length of thing
#'
#' Check the length of an object, be it list, vector, or table.
#'
#' @param obj
#'
#' @return \code{numeric} of length 1.
#' @keywords internal
#' @references https://github.com/grunwaldlab/metacoder
#' @noRd
length_of_thing <- function(obj) {
  if (is.data.frame(obj)) {
    return(nrow(obj))
  } else {
    return(length(obj))
  }
}

#' lapply with progress bars
#'
#' Imitates lapply with optional progress bars
#'
#' @usage progress_lapply(X, FUN, progress = interactive(), ...)
#' @param X The thing to iterate over
#' @param FUN The function to apply to each element
#' @param progress (logical of length 1) Whether or not to print a progress bar.
#'   Default is to only print a progress bar during interactive use.
#' @param ... Passed to function
#'
#' @return \code{numeric} of length 1.
#' @keywords internal
#' @references https://github.com/grunwaldlab/metacoder
#' @noRd
progress_lapply <- function (X, FUN, progress = interactive(), ...)
{
  if (progress) {
    progress_bar <- utils::txtProgressBar(min = 0, max = length(X),
                                          style = 3)
    one_iteration <- function(index) {
      output <- FUN(X[[index]], ...)
      utils::setTxtProgressBar(progress_bar, index)
      return(output)
    }
    output <- lapply(seq_len(length(X)), one_iteration)
    close(progress_bar)
  }
  else {
    output <- lapply(X, FUN, ...)
  }
  return(output)
}
