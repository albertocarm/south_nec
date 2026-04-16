#' Example dataset: anonymised Ki67 / peri-operative therapy cohort
#'
#' A cleaned, anonymised `data.frame` used by the Shiny app and the README
#' examples. Columns are:
#' \describe{
#'   \item{time_efs}{Event-free survival time (months).}
#'   \item{event_status}{Event indicator: 1 = event, 0 = censored.}
#'   \item{primary_surgery}{Factor (No/Yes): whether primary surgery was
#'     performed.}
#'   \item{disease_stage}{Factor (I/II/III).}
#'   \item{primary_site}{Factor (Colorectal/Other/Pancreas).}
#'   \item{ki67_percent}{Numeric 0-100, Ki67 proliferation index.}
#'   \item{periop_therapy}{0/1, receipt of peri-operative therapy.}
#' }
#'
#' Use [prep_data()] to derive the numeric dummies consumed by
#' [cure_model()] and the Shiny app.
#'
#' @format A data.frame with 174 rows and 7 columns.
#' @source Internal research dataset (anonymised).
"df_datos_app"
