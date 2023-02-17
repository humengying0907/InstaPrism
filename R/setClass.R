
#' An S4 class to represent input for InstaPrism when input.type='raw'
#' @slot phi.cs reference matrix for cell.states
#' @slot phi.ct reference matrix for cell.types
#' @slot bulk_mixture bulk RNA-seq to deconvolute
#' @slot map a list to store the correspondence between cell states and cell types
#' @keywords internal
#' @noRd
setClass("bpPrepare",
         slots = c(phi.cs='matrix',
                   phi.ct='matrix',
                   bulk_mixture='matrix',
                   map='list'))

#' An S4 class to represent full posterior information
#' @slot theta a cell state/type abundance matrix
#' @slot Z an array of cell state/type specific gene expression in each individual
#' @keywords internal
#' @noRd

setClass("posterior",
         slots = c(theta='matrix',
                   Z='array'))

#' An S4 class to represent posterior information of theta only
#' @slot theta a cell state/type abundance matrix
#' @keywords internal
#' @noRd
setClass('theta',
         slots=c(theta='matrix'))

setClassUnion("posterior.obj",
              c("theta","posterior"))

#' An S4 class to represent InstaPrism results under default setting (update=F)
#' @slot Post.ini.cs cell state abundance matrix (theta) and/or an S4 posterior object containing both theta and Z
#' @slot Post.ini.ct cell type abundance matrix (theta) and/or an S4 posterior object containing both theta and Z
#' @keywords internal
#' @noRd
#'
setClass("InstaPrism",slots = c(Post.ini.cs='posterior.obj',
                                Post.ini.ct='posterior.obj'))

#' An S4 class to represent extended InstaPrism results
#' @slot Post.ini.cs cell state abundance matrix (theta) and/or an S4 posterior object containing both theta and Z
#' @slot Post.ini.ct cell type abundance matrix (theta) and/or an S4 posterior object containing both theta and Z
#' @slot Post.updated.ct cell.type fraction estimates using updated reference
#' @keywords internal
#' @noRd
#'
setClass('InstaPrismExtra',slots = c(Post.ini.cs='posterior.obj',
                                     Post.ini.ct='posterior.obj',
                                     Post.updated.ct='posterior.obj'))




