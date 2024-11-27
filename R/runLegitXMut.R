#' Launch Shiny App for LegitXMut
#'
#' A function that launches the Shiny app for the LegitXMut package.
#' The app provides a user-friendly interactive interface for alignment and visulization
#' using functions from the LegitXMut package.
#'
#' @return No return value; launches a Shiny app in the browser.
#'
#' @examples
#' \dontrun{
#' LegitXMut::runLegitXMut()
#' }
#'
#' @references
#' Chang, W., J. Cheng, J. Allaire, C. Sievert, B. Schloerke, Y. Xie, J. Allen, J. McPherson, A. Dipert, B. Borges (2024).
#' shiny: Web Application Framework for R. R package version 1.9.1, https://CRAN.R-project.org/package=shiny
#'
#' Silva, Anjali. TestingPackage. GitHub, https://github.com/anjalisilva/TestingPackage. Accessed 5 Nov. 2024. \href{https://shiny.rstudio.com/}{Link}
#'
#' @export
#' @importFrom shiny runApp

runLegitXMut <- function() {
  # Locate the Shiny app directory in the package
  appDir <- system.file("shiny-scripts", package = "LegitXMut")

  # Check if the directory exists
  if (appDir == "") {
    stop("Could not find Shiny app directory. Please ensure the 'shiny-scripts' folder exists in the LegitXMut package.")
  }

  # Run the Shiny app
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")
  return(actionShiny)
}
# [END] written by Zhenghao Xiao
