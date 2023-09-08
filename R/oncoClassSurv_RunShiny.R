#' oncoClassSurv_RunShiny
#'
#' Perform convenient interactive operation through shiny application.
#'
#' @return Classifications summary, heat map plot and/or tables from prediction results.
#' @export
#' @import shiny
#' @importFrom data.table fread
#' @importFrom data.table data.table
#' @importFrom data.table is.data.table
#' @importFrom tibble column_to_rownames
#' @importFrom stats prcomp
#' @importFrom stats as.formula
#' @importFrom stats predict
#' @importFrom randomForest randomForest
#' @import e1071
#' @import ggplot2
#' @import ggfortify
#' @import patchwork
#' @import survival
#' @import survminer
#' @import BiocManager
#' @import limma
#' @import sva
#' @examples
#' \dontrun{
#' #Perform an interactive operation:
#' oncoClassSurv_RunShiny()
#' }
oncoClassSurv_RunShiny<-function(){
  shiny::runApp("./oncoClassSurvShinyAPP",display.mode = "auto")
}
