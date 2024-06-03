
#' Input Phi_cs object (cell states only) for InstaPrism
#' @slot phi.cs reference matrix for cell.states
#' @slot map a list with mapping information from cell.states to cell.types
#' @keywords internal
#' @noRd
setClass("refPhi_cs",
         slots = c(phi.cs='matrix',
                   map='list'))

#' Input Phi object for InstaPrism
#' @slot phi.cs reference matrix for cell.states
#' @slot phi.ct reference matrix for cell.types
#' @slot map a list with mapping information from cell.states to cell.types
#' @keywords internal
#' @noRd
setClass("refPhi",
         slots = c(phi.cs='matrix',
                   phi.ct='matrix',
                   map='list'))

#' An S4 class to represent output of bpPrepare function
#' @slot phi.cs reference matrix for cell.states
#' @slot phi.ct reference matrix for cell.types
#' @slot bulk_mixture bulk RNA-seq to deconvolute
#' @slot map a list with mapping information from cell.states to cell.types
#' @keywords internal
#' @noRd
setClass("bpPrepareObj",
         slots = c(phi.cs='matrix',
                   phi.ct='matrix',
                   bulk_mixture='matrix',
                   map='list'))

#' posterior information class
#' @slot theta a cell state/type abundance matrix
#' @slot Z an array of cell state/type specific gene expression in each individual
#' @keywords internal
#' @noRd
setClass("posterior",
         slots = c(theta='matrix',
                   Z='array'))

setClass('theta',
         slots=c(theta='matrix'))

setClassUnion("posterior.obj",
              c("theta","posterior"))


#' an S4 class object with scRNA based reference matrix
#'
#' @slot phi.cs matrix.
#' @keywords internal
#' @noRd
#'
setClass('initial_reference',slots = c(phi.cs = 'matrix'))

#' an S4 class object with scRNA based reference matrix
#'
#' @slot phi.cs matrix.
#' @slot phi.ct matrix.
#' @keywords internal
#' @noRd
#'
setClass('initial_reference_full',slots = c(phi.cs = 'matrix',
                                       phi.ct = 'matrix'))

setClassUnion("initial.reference.obj",
              c("initial_reference","initial_reference_full"))


#' InstaPrism output class
#' @slot Post.ini.cs cell state abundance matrix (theta)
#' @slot Post.ini.ct cell type abundance matrix (theta)
#' @slot map a list with mapping information from cell.states to cell.types
#' @slot initial.reference scRNA based reference matrix used for deconvolution
#' @slot initial.scaler a matrix with scaler information
#' @keywords internal
#' @noRd
#'
setClass("InstaPrism",slots = c(Post.ini.cs='posterior.obj',
                                Post.ini.ct='posterior.obj',
                                map='list',
                                initial.reference = 'initial.reference.obj',
                                initial.scaler = 'matrix'))

#' InstaPrism_update output class
#' @slot theta cell state/type abundance matrix
#' @slot psi_mal a gene * sample reference matrix for malignant cells,
#'    where each column represents the malignant reference for a corresponding sample;
#'    the psi_mal will be a NA matrix if no malignant cells present in the problem
#' @slot psi_env reference matrix for the non-malignant cells, which is shared by all the samples
#' @slot scaler a gene * sample matrix to reconstruct any cell.type specific expression of interest
#' @slot map a list with mapping information from cell.states to cell.types
#' @slot updated.cell.types a vector indicating cell.types with updated reference;
#'    equals to NULL if no environmental (non-malignant) cell.types being updated
#' @keywords internal
#' @noRd
#'
setClass('InstaPrism_update',slots = c(theta = 'matrix',
                                       psi_mal ='matrix',
                                       psi_env ='matrix',
                                       scaler = 'matrix',
                                       map = 'list',
                                       updated.cell.types = 'updated.cell.types',
                                       key = 'key'))

setClassUnion("updated.cell.types",
              c("character","NULL"))

setClassUnion("key",
              c("character","logical"))


############### S4 object related with InstaPrism_legacy function  #################
#' InstaPrism_legacy output class
#' @slot Post.ini.cs cell state abundance matrix (theta) and/or an S4 posterior object containing both theta and Z
#' @slot Post.ini.ct cell type abundance matrix (theta) and/or an S4 posterior object containing both theta and Z
#' @slot map a list with mapping information from cell.states to cell.types
#' @slot initial.reference scRNA based reference matrix used for deconvolution
#' @keywords internal
#' @noRd
#'
setClass("InstaPrism_legacy",slots = c(Post.ini.cs='posterior.obj',
                                       Post.ini.ct='posterior.obj',
                                       map='list',
                                       initial.reference = 'initial.reference.obj'))

#' S4 class to store (non-malignant) reference matrix phi or psi (if no malignant cells)
#' @description this corresponds to the refPhi class in BayesPrism package
#' @slot phi a matrix of dimension K*G,
#'		rownames are cell state/type names; colnames are gene IDs/names
#' @slot pseudo.min the desired minimum value used to normalize phi.
#' @keywords internal
#' @noRd
setClass("bpRefPhi",
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
#' @description this corresponds to the refTumor class in BayesPrism package
#' @slot psi_mal a matrix of dimension N*G, to denote updated sample-specific profiles (for malignant cells)
#'		rownames are bulk sample IDs; colnames are gene IDs/names
#' @slot psi_env a matrix of dimension (K-1)*G, to denote updated shared profiles (for non-malignant cells)
#'		rownames are non-malignant cell types; colnames are gene IDs/names
#' @slot key a character variable to denote the names for malignant cells
#' @slot pseudo.min the desired minimum value used to normalize phi.
#' @keywords internal
#' @noRd
setClass("bpRefTumor",
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


setClassUnion("updated_reference",
              c("bpRefPhi","bpRefTumor"))


#' An S4 class to represent extended InstaPrism results
#' @slot Post.ini.cs cell state abundance matrix (theta) and/or an S4 posterior object containing both theta and Z
#' @slot Post.ini.ct cell type abundance matrix (theta) and/or an S4 posterior object containing both theta and Z
#' @slot Post.updated.ct cell.type fraction estimates using updated reference
#' @slot initial.reference an initial_reference object with initial reference (scRNA reference phi)
#' @slot updated.reference updated reference object
#' @keywords internal
#' @noRd
#'
setClass('InstaPrismExtra',slots = c(Post.ini.cs='posterior.obj',
                                     Post.ini.ct='posterior.obj',
                                     Post.updated.ct='posterior.obj',
                                     map='list',
                                     initial.reference = 'initial.reference.obj',
                                     updated.reference = 'updated_reference'))

