\name{BikeSharing}
\docType{data}
\alias{BikeSharing}

\title{Bike Sharing Dataset (Processed)}

\description{
A processed version of the UCI Bike Sharing Dataset (hourly data). The data include
derived variables for part of day, quarter, and Fourier terms for cyclic effects
of hour and month.
}

\usage{
data("BikeSharing")
}

\format{
A data frame with 17,379 observations and 25 variables (original plus derived):

\describe{
  \item{instant}{Record index.}
  \item{dteday}{Date.}
  \item{season}{Season (1: spring, 2: summer, 3: fall, 4: winter).}
  \item{yr}{Year (0: 2011, 1: 2012).}
  \item{mnth}{Month (1--12).}
  \item{hr}{Hour (0--23).}
  \item{holiday}{Whether the day is a holiday (0/1).}
  \item{weekday}{Day of week (0--6).}
  \item{workingday}{Working day (0/1).}
  \item{weathersit}{Weather situation (1--4).}
  \item{temp}{Normalized temperature.}
  \item{atemp}{Normalized feeling temperature.}
  \item{hum}{Normalized humidity.}
  \item{windspeed}{Normalized wind speed.}
  \item{casual}{Count of casual users.}
  \item{registered}{Count of registered users.}
  \item{cnt}{Total count (casual + registered).}
  \item{hr_num}{Hour as numeric (same as \code{hr}).}
  \item{month_num}{Month as integer.}
  \item{part_of_day}{Factor: Night (0--5h), Morning (6--11h), Afternoon (12--17h), Evening (18--23h).}
  \item{quarter}{Factor: Q1--Q4.}
  \item{hr_sin}{Sine term for 24-hour cycle.}
  \item{hr_cos}{Cosine term for 24-hour cycle.}
  \item{mon_sin}{Sine term for 12-month cycle.}
  \item{mon_cos}{Cosine term for 12-month cycle.}
}
}

\source{
UCI Machine Learning Repository: Bike Sharing Dataset.
\url{https://archive.ics.uci.edu/ml/machine-learning-databases/00275/}
}

\examples{
data("BikeSharing")
head(BikeSharing)
dim(BikeSharing)
}
\keyword{datasets}
\concept{Bike sharing}
\concept{Count regression}
