% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/auto_units.R
\name{auto_units}
\alias{auto_units}
\title{Automatically select new time units}
\usage{
auto_units(y, from, target = 10, period_units = time_units)
}
\arguments{
\item{y}{A numeric vector of time values}

\item{from}{The original units of `y`}

\item{target}{The target value (order of magnitude) for the midpoint of rescaled time values. Default 10.}

\item{period_units}{A list of acceptable/understood time units. See Details. Default `time_units`.}
}
\value{
Character: Automatically-selected new time units, which will be one of `period_units`.
}
\description{
Given a vector of time values in original units, this function selects new
time units such that, when time is rescaled to the new units, the midpoint of
the time vector is as close to 10 as possible.
}
\details{
# Acceptable/understood time units in `period_units`

```
c("picoseconds",
"nanoseconds",
"microseconds",
"milliseconds",
"seconds",
"minutes",
"hours",
"days",
"weeks",
"months",
"years")
```
}
\author{
Caroline Ring
}
