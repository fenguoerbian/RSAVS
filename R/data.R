#' Portuguese student performance dataset
#'
#' A collection of datasets containing the grading point of evaluation of 1044
#' Portuguese students from two core classes(Mathematics and Portuguese).
#'
#' @format 
#' \itemize{
#'   \item \code{mat_df} is a data frame with 395 observations.
#'   \item \code{por_df} is a data frame with 649 observatoins.
#'   \item \code{full_df} is a combination of \code{mat_df} and \code{por_df}, resulting 1044 observations.
#' }
#' 
#' These data frames all have 33 variables:
#' \describe{
#'   \item{school}{character vector. Student's school 
#'   
#'   binary: \code{"GP"} - Gabriel Pereira or \code{"MS"} - Mousinho da Silveira}
#'   \item{sex}{character vector. Student's sex 
#'   
#'   binary: \code{"F"} - female or \code{"M"} - male}
#'   \item{age}{Interger, student's age, from 15 to 22}
#'   \item{address}{Character, student's home address type. 
#'   
#'   binary: \code{"U"} - urban or \code{"R"} - rural}
#'   \item{famsize}{Character vector, family size 
#'   
#'   binary: \code{"LE3"} - less or equal to 3 or \code{"GT3"} - greater than 3}
#'   \item{Pstatus}{Character vector, parent's cohabitation status 
#'   
#'   binary: "T" - living together or "A" - apart}
#'   \item{Medu}{Integer, mother's education 
#'   
#'     multinomial: \itemize{
#'     \item 0 - none,  
#'     \item 1 - primary education (4th grade), 
#'     \item 2 – 5th to 9th grade, 
#'     \item 3 – secondary education, 
#'     \item 4 – higher education
#'     }
#'   }
#'   \item{Fedu}{Integer, father's education 
#'   
#'     nominal: \itemize{
#'     \item 0 - none,  
#'     \item 1 - primary education (4th grade), 
#'     \item 2 – 5th to 9th grade, 
#'     \item 3 – secondary education, 
#'     \item 4 – higher education
#'     }
#'   }
#'   \item{Mjob}{Character vector, mother's job, 
#'   
#'   nominal: \itemize{
#'     \item "teacher" - teacher, 
#'     \item "health" - health care related, 
#'     \item "services" - civil services (e.g. administrative or police), 
#'     \item "at_home" - at home, 
#'     \item "other" - other
#'   }
#'   }
#'   \item{Fjob}{Character vector, father's job, 
#'   
#'   nominal: \itemize{
#'     \item "teacher" - teacher, 
#'     \item "health" - health care related, 
#'     \item "services" - civil services (e.g. administrative or police), 
#'     \item "at_home" - at home, 
#'     \item "other" - other
#'   }
#'   }
#'   \item{reason}{Character vector, reason to choose this school 
#'   
#'   nominal: \itemize{
#'     \item "home" - close to home, 
#'     \item "reputation" - school reputation, 
#'     \item "course"  - course preference, 
#'     \item "other" - other
#'   }}
#'   \item{guardian}{Character vector, student's guardian 
#'   
#'   nominal: "mother", "father" or "other"}
#'   \item{traveltime}{Integer, home to school travel time 
#'   
#'   nominal/incremental: \itemize{
#'     \item 1 - less than 15 min., 
#'     \item 2 - 15 to 30 min., 
#'     \item 3 - 30 min. to 1 hour, 
#'     \item 4 - greater than 1 hour
#'   }}
#'   \item{studytime}{Integer, weekly study time 
#'   
#'   nominal/incremental: \itemize{
#'   \item 1 --- less than 2 hours, 
#'   \item 2 --- 2 to 5 hours, 
#'   \item 3 --- 5 to 10 hours, 
#'   \item 4 --- greater than 10 hours}}
#'   \item{failures}{Integer, number of past class failures 
#'   
#'   nominal/incremental: n if 1 <= n < 3, else 4}
#'   \item{schoolsup}{Character vector, extra educational support 
#'   
#'   binary: "yes" or "no"}
#'   \item{famsup}{Character vector, family educational support. 
#'   
#'   binary: "yes" or "no"}
#'   \item{paid}{Character vector, extra paid classes within the course subject (Math or Portuguese)
#'   
#'   binary: "yes" or "no"}
#'   \item{activites}{Character vector,  extra-curricular activities 
#'     
#'   binary: "yes" or "no"}
#'   \item{nursery}{Character vector, attended nursery school 
#'     
#'   binary: "yes" or "no"}
#'   \item{higher}{Character vector, wants to take higher education 
#'     
#'   binary: "yes" or "no"}
#'   \item{internet}{Character vector, internet access at home 
#'     
#'   binary: "yes" or "no"}
#'   \item{romantic}{Character vector, with a romantic relationship 
#'     
#'   binary: "yes" or "no"}
#'   \item{famrel}{Integer, quality of family relationships 
#'   
#'   numeric: from 1 - very bad to 5 - excellent}
#'   \item{freetime}{Integer, free time after school 
#'   
#'   numeric: from 1 - very low to 5 - very high}
#'   \item{goout}{Integer, going out with friends 
#'   
#'   numeric: from 1 - very low to 5 - very high}
#'   \item{Dalc}{Integer, workday alcohol consumption 
#'   
#'   numeric: from 1 - very low to 5 - very high}
#'   \item{Walc}{Integer, weekend alcohol consumption 
#'   
#'   numeric: from 1 - very low to 5 - very high}
#'   \item{health}{Integer, current health status 
#'   
#'   numeric: from 1 - very bad to 5 - very good}
#'   \item{absences}{Integer, number of school absences 
#'   
#'   numeric: from 0 to 93}
#'   \item{G1}{Integer, first period grade
#'   
#'   numeric: from 0 to 20}
#'   \item{G2}{Integer, second period grade
#'   
#'   numeric: from 0 to 20}
#'   \item{G3}{Integer, final grade
#'   
#'   numeric: from 0 to 20}
#' }
#' 
#' @section Reference paper:
#'   P. Cortez and A. Silva. Using Data Mining to Predict Secondary School Student Performance. 
#'   In A. Brito and J. Teixeira Eds., Proceedings of 5th FUture BUsiness TEChnology Conference (FUBUTEC 2008) 
#'   pp. 5-12, Porto, Portugal, April, 2008, EUROSIS, ISBN 978-9077381-39-7.
#'   
#'   This paper is available at \url{http://www3.dsi.uminho.pt/pcortez/student.pdf}
#' @source \url{https://archive.ics.uci.edu/ml/datasets/Student+Performance}
#' @name Student_performance_dataset
#' @aliases full_df mat_df por_df
NULL
