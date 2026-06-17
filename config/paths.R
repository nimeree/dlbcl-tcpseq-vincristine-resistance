# Portable path helpers for thesis analysis scripts.
#
# Set these environment variables before running scripts, or edit the defaults
# to match your local data layout. The GitHub repository does not include raw
# sequencing files or generated analysis outputs.

analysis_root <- function() {
  Sys.getenv(THESIS_ANALYSIS_DIR, unset = file.path(data, analysis))
}

input_root <- function() {
  Sys.getenv(THESIS_INPUT_DIR, unset = file.path(data, input))
}

project_resource_root <- function() {
  Sys.getenv(THESIS_RESOURCE_DIR, unset = file.path(data, external))
}

external_root <- function() {
  Sys.getenv(THESIS_EXTERNAL_ROOT, unset = file.path(data, external))
}

analysis_path <- function(...) file.path(analysis_root(), ...)
input_path <- function(...) file.path(input_root(), ...)
project_resource_path <- function(...) file.path(project_resource_root(), ...)
external_path <- function(...) file.path(external_root(), ...)
