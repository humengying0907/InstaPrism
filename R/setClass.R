
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



#' S4 class to store (non-malignant) reference matrix phi or psi (if no malignant cells)
#'
#' @slot phi a matrix of dimension K*G,
#'		rownames are cell state/type names; colnames are gene IDs/names
#' @slot pseudo.min the desired minimum value used to normalize phi.
#' @keywords internal
#' @noRd
setClass("refPhi",
         slots = c(
           phi = "matrix",
           pseudo.min ="numeric"
         ),
         prototype = list(
           phi = matrix(),
           pseudo.min = NA_real_
         ),
         validity = function(object){
           errors <- character()

           if(length(object@pseudo.min) != 1 | !is.numeric(object@pseudo.min)){
             msg <- paste("invalid pseudo.min")
             errors <- c(errors, msg)
           }

           phi.min <- min(object@phi)
           if(phi.min < 0){
             msg <- paste("reference contain negative values")
             errors <- c(errors, msg)
           }
           if(!is.na(object@pseudo.min) & phi.min != object@pseudo.min)
             warning("Warning: pseudo.min does not match min(phi)")
           if (length(errors) == 0) TRUE else errors
         }
)

#' S4 class to store reference matrix psi for malignant cells
#'
#' @slot psi_mal a matrix of dimension N*G, to denote updated sample-specific profiles (for malignant cells)
#'		rownames are bulk sample IDs; colnames are gene IDs/names
#' @slot psi_env a matrix of dimension (K-1)*G, to denote updated shared profiles (for non-malignant cells)
#'		rownames are non-malignant cell types; colnames are gene IDs/names
#' @slot key a character variable to denote the names for malignant cells
#' @slot pseudo.min the desired minimum value used to normalize phi.
#' @keywords internal
#' @noRd
setClass("refTumor",
         slots = c(
           psi_mal = "matrix",
           psi_env = "matrix",
           key = "character",
           pseudo.min ="numeric"
         ),
         prototype = list(
           psi_mal = matrix(),
           psi_env = matrix(),
           key = NA_character_,
           pseudo.min = NA_real_
         ),
         validity = function(object){
           errors <- character()

           if(length(object@key) != 1 | !is.character(object@key)){
             msg <- paste("invalid key")
             errors <- c(errors, msg)
           }

           psi_mal_genes <- colnames(object@psi_mal)
           psi_env_genes <- colnames(object@psi_env)
           if(!identical(psi_mal_genes, psi_env_genes)){
             msg <- paste("Gene names of psi_mal and psi_env do not match")
             errors <- c(errors, msg)
           }

           if(length(object@pseudo.min) != 1 | !is.numeric(object@pseudo.min)){
             msg <- paste("invalid pseudo.min")
             errors <- c(errors, msg)
           }

           phi.min <- min(object@psi_mal, object@psi_env)
           if(phi.min < 0){
             msg <- paste("reference contain negative values")
             errors <- c(errors, msg)
           }
           if(!is.na(object@pseudo.min) & phi.min != object@pseudo.min)
             warning("Warning: pseudo.min does not match min(phi)")
           if (length(errors) == 0) TRUE else errors
         }

)
