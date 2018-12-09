\name{titanic_comp}
\docType{data}
\alias{titanic_comp}
\title{Titanic passenger data}
\description{Data about the Titanic passenger, complete cases only.}

\usage{data(titanic_comp)}
\format{
  A data frame with the following columns:

   \describe{
    \item{Survived}{Whether a passenger survived.}
    \item{Pclass}{The class of the passenger's ticket. A factor with 3 levels.}
    \item{Sex}{Male or female. A factor with 2 levels.}
    \item{Age}{The age of the passenger.}
    \item{SibSp}{The number of siblings/spouse aboard.}
    \item{Parch}{The number of parents/children aboard.}
    \item{Fare}{The passenger's fare.}
    \item{Embarked}{The passenger's port of embarkation. A factor with 3 levels.}
  }
}


\references{https://www.kaggle.com/c/titanic/data}
\examples{
data(titanic_comp)
}

