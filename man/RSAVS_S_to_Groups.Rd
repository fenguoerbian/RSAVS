% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pre_define_functions.r
\name{RSAVS_S_to_Groups}
\alias{RSAVS_S_to_Groups}
\title{Convert pair-wise difference vector to subgroup structure}
\usage{
RSAVS_S_to_Groups(s_vec, n)
}
\arguments{
\item{s_vec}{the s vector(pair-wise difference vector), length \code{n * (n - 1) / 2} and 
\eqn{s_{ij} = \mu_i - \mu_j}.}

\item{n}{number of observations.}
}
\value{
a list containing the grouping result.
}
\description{
This function converts the augmented pair-wise difference vector(the s vector in the algorithm) 
to the subgrouping result.
}
\details{
Detailed definition of the s vector(length \code{n * (n - 1) / 2}) can be found in the paper.
}
\section{Note}{

  It's possible that there's logical contradiction in the s vector, 
  e.g. \code{s_{12} = s_{13} = 0}, but \code{s_{23} != 0}, especially when the algorithm is
  configured with a loose tolerance. In this function, \code{i} and \code{j} would be put in the
  same subgroup as long as there is a path that connects them. In previous example, 1st, 2nd and
  3rd observations will be classified as in the same subgroup even though \code{s_{23} != 0}. 
  
  This strategy would presumably provide a more concentrate results, hence less number of subgroups. 
  But of course it tends to merge subgroups into big ones and leaving some alone observations.
  
  For large scale data, especially when the number observation is big. It will be difficult to save
  all \code{s_vec}s during the algorithm considering the algorithm have to search over the 
  \code{lam1_length * lam2_length} grid. During the algorithm, the \code{s_vec} is utilized to
  improve the subgrouping results then discarded.
}

\examples{
n <- 10    # number of observations
group_center <- c(-1, 0, 1)    # three group centers
# subgroup effect vector    
alpha_true <- sample(group_center, size = n, replace = TRUE)
d_mat <- RSAVS:::RSAVS_Generate_D_Matrix(n)    # pair-wise difference matrix
s_vec <- d_mat \%*\% alpha_true
RSAVS:::RSAVS_S_to_Groups(s_vec, n)

# if there's contradiction in s_vec
s_vec[1] <- 0    # s_{12} = 0
s_vec[2] <- 0    # s_{13} = 0
s_vec[n] <- 1    # s_{23} != 0
# 1, 2, 3 will still be in the same subgroup
RSAVS:::RSAVS_S_to_Groups(s_vec, n)
}
