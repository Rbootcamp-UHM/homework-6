## Mechanicʻs rule to calculate square roots
## 
mysqrt <- function ( x, guess=1, tol=0.001, debug=F) {
	error <- tol+1
	while( tol < error) {
	  newguess <- x/guess
	  if (debug) print(newguess)	  
	  guess <- (guess+newguess)/2
	  error <- abs( newguess^2-x )
	}
	return(newguess)
}

## Calculating the derivative based on a vector of xʻs and y=fx ʻs
## using the central difference formla. Returns a vector of
## the derivative function fʻx of the same length as x and y
#
calc_deriv <- function ( x, fx ) {
## letʻs get the derivative for internal points 1 from start to 1 from end
   fxminus1 <- fx[1:(length(fx)-2)]
   fxplus1 <- fx[3:length(fx)]  
   xminus1 <- x[1:(length(x)-2)]
   xplus1 <- x[3:length(x)]  
   deriv <- (fxplus1 - fxminus1) / (xplus1 - xminus1) 
   return( c( deriv[1], deriv, deriv[length(deriv)] ) )
}

## An internal function that finds for a given value of x, the index 
## position corresponding to that particular value of x 
# 
find_index <- function( x, point ) {
  ii <- which(x > point)[1]
  if ( is.na(ii) )  ii <- length(x)
  deltaii <- abs(point - x[ii])
  deltaii_1 <- abs(point - x[ii-1])
  if ( deltaii_1 < deltaii )  ii <- ii-1
  return(ii)
}

## Newtonʻs method for finding roots using the first derivative
## 
find_root <- function( guessx, x=x, fx=fx ) {
		
  fpx <- calc_deriv( x, fx )
  i <- find_index(x, guessx)
  oldi = -1

  while( oldi != i ) {  ## test for convergence on same index value
    oldi <- i
    x2 <- x[i] - fx[i]/fpx[i]  ## Newtonʻs root-finding formula 
    i <- find_index(x, x2)
  }	
#  return(data.frame(xest=x2, x=x[i], fx=fx[i], i=i))
   return(x2)
}


## Finds all roots calling find_root on multiple guesses
## ideally spread over the range of x, fx
##
return_roots <- function( guesses, x, fx) { 
	roots <- sapply( guesses, function(i) find_root(guessx=i, x, fx))
	return(unique(roots))
}


## Example of use to find roots 
x <- seq(-5, 5, by=.001)
fx <- x^3 - 8*x - 4
roots <- return_roots(-5:5, x, fx)
residue = roots^3 - 8*roots - 4
print(cbind( roots,residue) )
#          roots       residue
# [1,] -2.534070 -3.739315e-08
# [2,] -0.517304 -1.430850e-07
# [3,]  3.051374  1.282059e-06


## Part 2 - Comparing Mechanicʻs Rule and Newtonʻs Rule
###############################################################
## We will find roots of 
##		numbers <- c(64, 10, .3, 3.14159, 9, 4, 2, 1, 28, 50)
##  	along the range 0-10:  x <- seq(0, 10, by=.001)
##

## myfunc:
## Creating a square or cube root function so that we can use 
## Newtonʻs method to find the roots
## e.g., to find square root of 64: set 64 = x^2, fx = 64 - x^2
## 64 = number is the number of interest
## 2 = is the power to solve for (e.g., square)
## x is the root which is found when fx=0 

myfunc <- function( number, x, exp) {
	return (number - x^exp)
}

## Creates one line for a table of solutions for each root, 
## comparing Newtonʻs to Mechanicʻs rule for square roots
## and Newtonʻs to the built-in R function ^(1/exp) for other roots

compare_roots <- function ( number, x, guesses, exp=2, file="rtout.txt" ) {
	fx <- myfunc( number, x, exp)
	roots <- return_roots(guesses, x, fx)
    if (exp==2) comproot=mysqrt(number) else comproot=number^(1/exp) 
	cat(number,"\t", 
		format(roots, scientific=T),"\t",
		format(comproot, scientific=T),"\t",
		format(roots - comproot, scientific=T), "\n",
		file=file, sep="", append=TRUE
	)
}	

## Print the header for the table
##
print_header <- function ( title="Square roots", comp="Mechanic", fileout="rtout.txt") {
	cat(title, "\n", file=fileout, append=TRUE)
	cat("Number","\t", 
		"Newton","\t",
		comp,"\t",
		"difference", "\n",
		file=fileout, sep="", append=TRUE
	)	
}

## Comparing Newtonʻs rule and Mechanicʻs rule for square roots
## By repetitively calling the function

print_header( "Square roots", "Mechanic", file="HW6_rtout.txt")
x <- seq(0, 10, by=.001)
compare_roots( 64, x, 1:10, 2)	## Example: sqrt of 64, range x=1:10 
compare_roots( 10, x, 1:10, 2)	
compare_roots( .3, x, 1:10, 2)	
compare_roots( 3.14159, x, 1:10, 2)	
compare_roots( 9, x, 1:10, 2)	
compare_roots( 4, x, 1:10, 2)	
compare_roots( 2, x, 1:10, 2)	
compare_roots( 1, x, 1:10, 2)	
compare_roots( 28, x, 1:10, 2)	
compare_roots( 50, x, 1:10, 2)	

## Comparing Newtonʻs rule and R-functions for cube roots
## By repetitively calling the function

print_header( "\nCube Roots", "x^(1/3)", file="HW6_rtout.txt")
compare_roots( 64, x, 1:10, 3)	
compare_roots( 10, x, 1:10, 3) ## Example: cube root of 10, range x=1:10 
compare_roots( .3, x, 1:10, 3)	
compare_roots( 3.14159, x, 1:10, 3)	
compare_roots( 9, x, 1:10, 3)	
compare_roots( 4, x, 1:10, 3)	
compare_roots( 2, x, 1:10, 3)	
compare_roots( 1, x, 1:10, 3)	
compare_roots( 28, x, 1:10, 3)	
compare_roots( 50, x, 1:10, 3)	

## Same as above but using sapply
## Comparing Newtonʻs rule vs Mechanicʻs rule or x^(1/3) 

numbers <- c(64, 10, .3, 3.14159, 9, 4, 2, 1, 28, 50)
x <- seq(0, 10, by=.001)

print_header( "Square roots", "Mechanic", file="HW6_rtout2.txt")
sapply(numbers, function(number) compare_roots(number, x, 1:10, 2, file="HW6_rtout2.txt" ))

print_header( "\nCube Roots", "x^(1/3)", file="HW6_rtout2.txt")
sapply(numbers, function(number) compare_roots(number, x, 1:10, 3, file="HW6_rtout2.txt" ))
