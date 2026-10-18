\name{GetMetacellSEACells}
\alias{GetMetacellSEACells}
\title{Access metacell results stored by MetacellsByGroupsSEACells}
\usage{
GetMetacellSEACells(seurat_obj, wgcna_name = "SEACells")
}
\arguments{
\item{seurat_obj}{A Seurat object processed by
\code{\link{MetacellsByGroupsSEACells}}.}
\item{wgcna_name}{Name under which results were stored. Default
\code{"SEACells"}.}
}
\value{
The stored list: \code{metacell_object}, \code{soft_assignments},
\code{qc}, \code{fit_stats}, \code{zerorate_overall},
\code{dubrate_overall}, \code{score_overall}, \code{params},
\code{skipped_groups}, \code{call}.
}
\description{
Access metacell results stored by MetacellsByGroupsSEACells.
}
