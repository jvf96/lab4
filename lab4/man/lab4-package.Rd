\name{lab4-package}
\docType{package}
\title{
\packageTitle{lab4}
}
\description{
\packageDescription{lab4}
}
\details{

The DESCRIPTION file
\packageDESCRIPTION{lab4}
\packageIndices{lab4}
create a package to handle linear regression models
}
\author{ packageAuthor {Jorge Villar Fuentes}

Maintainer packageMaintainer{Jorge Villar Fuentes}
}

\keyword{ package }
\seealso{https://en.wikipedia.org/wiki/Greatest_common_divisor.
}
\examples{
linreg_ex <- linreg$new(Petal.Length~Sepal.Width+Sepal.Length, data=iris)
linreg_ex$summary()

