% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_metadata_from_file.R
\name{load_metadata_from_file}
\alias{load_metadata_from_file}
\title{load_metadata_from_file(xlsxFile = NULL, phenotable_tsv = NULL, phenolabels_tsv = NULL)}
\usage{
load_metadata_from_file(
  opt = NULL,
  xlsxFile = NULL,
  phenotable_tsv = NULL,
  phenolabels_tsv = NULL,
  class_to_ignore = "N_A"
)
}
\description{
Loads an Excel(TM) file or tsv files and returns a list containing 1) the metadata itself and 2) description of label types for the metadata, provided the Excel file contains these data on sheets 1 and 2, respsctively.
}
