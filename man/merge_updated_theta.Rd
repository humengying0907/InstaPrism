% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InstaPrism.R
\name{merge_updated_theta}
\alias{merge_updated_theta}
\title{merge cell.state level information to cell.type level}
\usage{
merge_updated_theta(InstaPrism_updated_obj)
}
\arguments{
\item{InstaPrism_updated_obj}{an InstaPrism_updated_obj from InstaPrism_update() function}
}
\value{
fraction estimation at cell.type level
}
\description{
InstaPrism_updated_obj contains cell.state level fraction estimation when keep.phi = 'phi.cs',
    this function merges cell.state level theta to cell.type level.
}
