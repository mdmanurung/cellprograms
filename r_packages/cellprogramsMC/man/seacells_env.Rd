\name{seacells_env}
\alias{seacells_env}
\docType{data}
\title{Basilisk environment for SEACells}
\format{A \code{basilisk::BasiliskEnvironment} object.}
\source{
Created at package load; see \code{cellprogramsMC-package}.
}
\usage{
seacells_env
}
\description{
The Python environment used by \code{\link{MetacellsByGroupsSEACells}} to
run the SEACells aggregation backend. The environment is created
automatically on first use (requires network access at that point).

SEACells 0.3.3 (PyPI) requires Python >= 3.11; its core dependencies are
CPU-only (numpy, pandas, scipy, scanpy, anndata, palantir, numba,
scikit-learn, pyranges, matplotlib, seaborn, tqdm) and are resolved by pip.
}
\keyword{datasets}
