\name{geo.legend}
\alias{geo.legend}
\alias{geo.palette}
\title{Adds a geological (or other temporal) legend to a plotted tree}
\usage{
geo.legend(leg=NULL, colors=NULL, alpha=0.2, ...)
geo.palette()
}
\arguments{
	\item{leg}{a matrix with the starting & ending point of each plotted era in rows, & names of the time periods as rownames.}
	\item{colors}{a vector of colors for the time periods of the rows in \code{leg}.}
	\item{alpha}{transparency level to apply to \code{colors}.}
	\item{...}{optional arguments.}
}
\description{
	Adds a geological legend to a plotted tree.
}
\details{
	The function \code{geo.legend} adds a geological (or other temporal) legend to a plotted tree.
	
	The function \code{geo.palette} returns a geological time color palette to the user.
}
\value{
	\code{geo.legend} adds a visual element to a plotted tree and invisible returns an object of class \code{geo.legend} containing the time periods and colors of the painted legend.
	
	\code{geo.palette} simply returns a geological timescale color palette as an object of class \code{"geo.palette"}.
}
\references{
	Revell, L. J. (2024) phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). \emph{PeerJ}, \bold{12}, e16505.
}
\examples{
## load tree
data(anoletree)

## rescale tree to 50 ma total depth
anoletree<-rescale(anoletree,model="depth",depth=50)

## plot phylogeny
plotTree(anoletree,ylim=c(-0.16,1)*Ntip(anoletree),
    ftype="i",lwd=1,fsize=0.5)

## add geological color legend
obj<-geo.legend()

## reset par to default
par(mar=c(5.1,4.1,4.1,2.1))
}
\author{Liam Revell \email{liam.revell@umb.edu}}
\keyword{phylogenetics}
\keyword{plotting}
