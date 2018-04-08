#################################################################################
##  Comparision of Numerical Integration by Trapzoidal method and Simpsonʻs 1/3 rule
##
##  two main functions: trap_area and simpson_area
##  inputs: x = a vector
##   		fx = a vector representing a given function of x
##  return values: area, a numerical value 
##
##  utility functions: 
##			funkx: returns vector fx for any user-supplied expression (e.g. 1/(x+2))
##			subsample_xii: 	input: x= data vector, n= number of intervals desired
##							ouput: index positions for the n intervals
##			get_int: called by subsample_xii
##					input: x, n
##					output: returns number of index positions (width) to produce n 
##						intervals given x
##			find_deltas: returns a vector of differences between adjacent values
##					input: vector x
##					output: x[n]-x[n-1], ... , x[3]-x[2], x[2]-x[1]
##							ouput vector is length(x)-1
#################################################################################

############################ Function Definitions ###############################

############################ funcx
## Function to obtain f(x) values for y=1/(2+x)
## need to write one for each function tested
############################
funcx <- function( x ) { # x = 0 to 2
	fx = 1/(2+x)
	return (fx)
}

############################ funkx
## Alternatively can write a generic function to return fx for any user-supplied
## function - use "eval()" within a function definition
#
# Must assign x within function or it doesn't get passed to local scope
############################
funkx <- function( expression, x ) {
	  x <- x             ## assignment inside function introduces x to local scope            
	  return(eval(expression))  # evaluates fx for xʻs, returns fx
}  	
x <- seq(0,2, by=.001)  
fx <- funkx(1/(2+x), x) # an example of using funkx


############################ get_int and subsample_xii
## Finds the index values for taking n intervals along 
## vectors x and fx
############################
get_int <- function(x, n=2) { # for a given number of intervals n, return # indexʻs
	return( (length(x)-1)/n )  
}

subsample_xii <- function(x, n=2) { # returns indices of x, fx for n intervals
	int <- get_int( x, n )
	ii <- round(seq(1, length(x), by=int))
	return(ii)
}


############################ find_deltas
## Calculates difference between each x and its 
## previous value (x[4]-x[3], x[3]-x[2], x[2]-x[1], etc.)
## returns a vector one shorter than original vector
## used in several functions
############################
find_deltas <- function( x ) {
  iio <- 1: (length(x)-1)  # index from 1 to n-1
  oii <- 2: length(x)      # index from 2 to n
  deltas <- x[oii] - x[iio]		
  return(deltas)
}

############################ trap_area
## Calculates area using trapezoidal rule
## Strategy is to make two vectors: 
## right side of trapezoid (fx[oii]) & the left side (fx[iio])
## h is the width (x) of the intervals
############################
trap_area <- function( x, fx ) {
    iio <- 1: (length(x)-1)  # index from 1 to n-1
    oii <- 2: length(x)      # index from 2 to n
    h <- x[oii] - x[iio]     # widths of each interval
    h <- find_deltas(x)		# widths of each interval
    area <- h * (fx[iio] + fx[oii])/2 # widths of intervals * average heights
    return( sum(area) )
}

############################ simpson13_area
## simpson13_area: Calculates area using Simpsons 1/3 rule
## which calculates the parabola fit between each set of 3 adjacent points
## Strategy is to make coefficient vector (c) depending on number
## of mid-intervals (the number of intervals after removing 1st & last) 
## and multiplying c to fx. These are predicted by parabola geometry 
## h is the width (x) of the intervals
############################
simpson13_area <- function( x, fx ) {
    oiio <- 2: (length(x)-1)  # index from 2 to n-1
    nmid_intervals <-  (length(oiio)-1)/2 # number of middle intervals
    c <- c( 1, 4, rep(c(2,4), times=nmid_intervals), 1) # coefficients for fxʻs
    									# 1,4,1 
    									# 1,4,2,4,1 
    									# 1,4,2,4,2,4,1 etc.    
    h <- find_deltas(x)		# widths of each interval
    area <- mean(h/3) * ( fx*c ) # widths of intervals/3 * fx * parabola coefficients
    return( sum(area) )
}

ii <- subsample_xii(x, n=2) # run once for n=2, 4, 8, 16, 32, etc.
trap_area(x[ii], fx[ii])
simpson13_area(x[ii], fx[ii])

############################ compare_methods
## Can replace above with the sapply() below over the intervals
## indicated by a vector c(2,4,...128)
## and put in a function to be able to compare different input functions
## This function runs the analyses and creates the output table
############################
compare_methods <- function( x, fx ) { # runs trapezoidal and Simpsons methods 
  out <- sapply( c(2,4,8,16,32,64,128), function(n) { # over 2,4,...128 intervals
  	ii <- subsample_xii(x, n)
	xx <- x[ii]
	fxx <- fx[ii]
	t_area = trap_area(xx, fxx)
	s_area = simpson13_area(x[ii], fx[ii])
	c( 	n=n, 
		trap_area = t_area,
		simpson_area = s_area,
		diff = t_area-s_area
	  ) 							# collect results in output vector
	 }
    )

	out <- data.frame( t(out))	# converts output to dataframe
	out$delta_trap_area <- c(NA, find_deltas( out$trap_area)) # incremental improvement
	out$delta_simpson_area <- c(NA, find_deltas( out$simpson_area)) 
	return(out)
}

############################ Begin of Execution ###############################

x <- seq(0,2, by=.001) 			# x = 0 to 2
fx <- funkx(1/(2+x), x)
out <- compare_methods(x, fx) 	# area 0.69315
cat("y=1/(2+x)","\n", file="HW6_integration.csv")
write.table(out, file="HW6_integration.csv", sep=",", row.names=FALSE, append=TRUE)

x <- seq(0,1, by=.001) 			# x = 0 to 1
fx <- funkx(1/(1+x^2), x)
out <- compare_methods(x, fx) 	# area pi/4=0.7853982
cat("y=1/(1+x^2)","\n", file="HW6_integration.csv", append=TRUE)
write.table(out, file="HW6_integration.csv", sep=",", row.names=FALSE, append=TRUE)

x <- seq(0,.4, by=.0005) 		# x=0 to 0.4
fx <- funkx(exp(-x^2/2), x)
out <- compare_methods(x, fx) 	# area 0.389
cat("y=exp(-x^2/2)","\n", file="HW6_integration.csv", append=TRUE)
write.table(out, file="HW6_integration.csv", sep=",", row.names=FALSE, append=TRUE)

x <- seq(-1,1, by=.001) 					# x=-1 to 1
fx <- funkx((x^2 - cos(x))*exp(-x), x)
out <- compare_methods(x, fx) 			# area -1.054 
cat("y=(x^2 - cos(x))*exp(-x)","\n", file="HW6_integration.csv", append=TRUE)
write.table(out, file="HW6_integration.csv", sep=",", row.names=FALSE, append=TRUE)
