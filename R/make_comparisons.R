#' Make vecotr of comparisons for EdgeR
#'
#' Generate a vector of all comparisons in a list for edgeR from a vector of all
#' qualities to compare by, seperated by an underscore
#'
#' @param grouper vector to generate comaprisons from
#' @keywords EdgeR
#' @export
#' @return vector of comparisons for use with edgeR
#' @details Proud of this!
#' MakeComparisons()
make_comparisons <- function(grouper){

  unlist(sapply(unique(unlist(str_split(grouper, ###For every cell type or date
                                       "_",
                                       n = 2))),
               FUN = function(λ){ ##Filter groups for that input
                 filtered <- unique(grouper[str_detect(grouper,λ)])

                 sapply(filtered, ## For these filtered inputs
                        FUN = function(þ) { ##Make a comparison with all other options
                          sapply(filtered, FUN = function(ω){
                            paste(þ,
                                   ω,
                                   sep = " - "
                            )
                          })
                        })
               }
  ))
}
