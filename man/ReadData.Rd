\name{ReadData}
\alias{ReadData}
\title{
Function for preparing data object
}

\description{
ReadData takes data such as matrix, labels, and platform information, and produce data object to be used in the down stream analysis, such as filtering genes.}
\usage{
ReadData(Data, Labels, Platform = NULL, verbose = TRUE)
}

\arguments{
  \item{Data}{
a dataframe, matrix, or ExpressionSet with values to be used in the down stream analysis. Samples as columns  and rows genes/features as rows. Matrix should has column names and row names.}
  \item{Labels}{
a vector indicating the classes of the samples.
Should be with the same length of the columns number in data.
This can be a variable name stored in the ExpressionSet if ExpressionSet is used.
}
  \item{Platform}{
Optional, vector with the same length of labels indicating.
This can be a variable name stored in the ExpressionSet if ExpressionSet is used.
}
 \item{verbose}{
a logical value indicating whether processing messages will be printed or not. Default is TRUE.
}
}

\value{
data object multiclassPairs_object
 \item{Data}{
dataframe (gene as rows and samples as columns)
}
 \item{Labels}{
 a vector containing classes information
}
 \item{Platform}{
 	a vector containing Platform information, or NULL if no input is used
}

}

\author{
Nour-al-dain Marzouka <nour-al-dain.marzouka at med.lu.se>
}

\examples{
Data <- matrix(runif(100000), nrow=100, ncol=100,
               dimnames = list(paste0("G",1:100), paste0("S",1:100)))

L <- sample(x = c("A","B","C"), size = 100, replace = TRUE)

P <- sample(x = c("P1","P2"), size = 100, replace = TRUE)


table(P,L)

object <- ReadData(Data = Data,
                   Labels = L,
                   Platform = P,
                   verbose = FALSE)
object
}
