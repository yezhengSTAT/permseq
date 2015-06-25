\docType{methods}
\name{Prior-method}
\alias{plot,Prior-method}
\alias{print,Prior-method}
\alias{show,Prior-method}
\alias{summary,Prior-method}
\alias{names,Prior-method}
\alias{[,Prior-method}
\alias{[[,Prior-method}
\alias{[<-,Prior-method}
\alias{[[<-,Prior-method}
\title{"Prior" class methods}
\usage{
\S4method{plot}{Prior}(x, y, ...)
\S4method{print}{Prior}(x)
\S4method{show}{Prior}(object)
\S4method{summary}{Prior}(object)
\S4method{names}{Prior}(x)
}
\arguments{
  \item{x}{Prior object}
  \item{y}{Other parameter}
  \item{...}{Other parameters}
  \item{object}{Prior object}
}

\description{

The package defines the following methods for the "Prior" class:
\describe{
  \item{\code{plot}}{Plot the DNase (or Histone selected to be used as DNase) read counts versus average over uni-read ChIP-seq read counts.}
  \item{\code{print}}{Print out all the information within slots of a "Prior" class.}
  \item{\code{show}}{Show the selected important information within a "Prior" class.}
  \item{\code{summary}}{Alignment details report of DNase, Histone and/or ChIP-seq files from bowtie.}
  \item{\code{names}}{Names of the individual slots in "Prior" object.}
  \item{\code{[}}{ Extract part of the Prior object.}
  \item{\code{[[}}{ Extract part of the Prior object.}
\item{\code{[<-}}{Set value to the extract part of the Prior object.}
\item{\code{[[<-}}{Set value to the extract part of the Prior object.}
}

}

\author{
Ye Zheng, M. Constanza Rojo-Alfaro.
}
