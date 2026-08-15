\name{Random values simulation from the zero adjusted Dirichlet distribution}
\alias{rzad}
\title{
Random values simulation from the zero adjusted Dirichlet distribution
}

\description{
Random values simulation from the zero adjusted Dirichlet distribution.
}

\usage{
rzad(n, phi, mu, patterns)
}

\arguments{
\item{n}{
The number of compositional vectors to simulate.
}
\item{phi}{
The concentration parameter.
}
\item{mu}{
The mean vector.
}
\item{patterns}{
A matrix with the patterns of zeros, where the value of 0 indicates the presence of a zero,
and the last column contains the percentage of occurrence each pattern.
}
}

\details{
Random values from the zero adjusted Dirichlet distribution are generated.
}

\value{
A matrix with values generated from the zero adjusted Dirichlet distribution.
}

\references{
Tsagris M. and Stewart C. (2018). A Dirichlet regression model for compositional data with zeros.
Lobachevskii Journal of Mathematics, 39(3): 398--412.

Preprint available from https://arxiv.org/pdf/1410.5011.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{zad.mle}, \link{zadr}
}
}

\examples{
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
y[sample(1:450, 15) ] <- 0
y <- y / rowSums(y)
mod <- zad.mle(y)
phi <- mod$phi
mu <- mod$mu
patterns <- mod$patterns
x <- rzad(100, phi, mu, patterns)
}

