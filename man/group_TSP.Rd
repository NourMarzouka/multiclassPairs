\name{group_TSP}
\alias{group_TSP}
\title{
Internal function: for grouping labels for one-vs-rest usage}
\description{
Used to convert labels to factor to be used by switchBox package.
}
\usage{
group_TSP(label, my_group)
}
\arguments{
  \item{label}{
a vector indicating multi classes}
  \item{my_group}{
character indicate the wanted class
}
}

\value{
a factor contains two levels one is the wanted class and the other is "rest" that represent any other class other than the wanted class
}

\author{
Nour-al-dain Marzouka <nour-al-dain.marzouka at med.lu.se>
}

\examples{
L <- sample(x = c("A","B","C","D"), size = 1000, replace = TRUE)

}
