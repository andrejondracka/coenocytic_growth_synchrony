library(dplyr)

initcon <- c(rep(0,328),rep(1,1116),rep(2,270),rep(3,12),rep(4,3),rep(5,1),rep(6,1),rep(7,2)) ###sets a pool of initial conditions from the distribution identical to t=0 in the real data
vari <- 0.50
mean <- 11

###function that generates truncated normal distribution (to avoid having negative durations of intervals)
meanfunction <- function(mean=11, var) {
  randommean = rnorm(1,mean,mean*var)
  while (randommean<0) {
    randommean = rnorm(1,mean,mean*var)}
  return(randommean)
}

###function that generates a time history of one cell. sets up starting number of nuclei from the specified initial condition status, and sets the duration of the nucler cycle based on specified mean and variance (time units = hours). in the end, records 50 hours of cell history. 
generatecellhistory <- function(initcon=initcon,mean=mean,var=var) {

a <- sample(initcon,1)
d1 <- rep(a,round(meanfunction(11,var)))
d2 <- rep(a+1,round(meanfunction(11,var)))
d3 <- rep(a+2,round(meanfunction(11,var)))
d4 <- rep(a+3,round(meanfunction(11,var)))
d5 <- rep(a+4,round(meanfunction(11,var)))
d6 <- rep(a+5,round(meanfunction(11,var)))
d7 <- rep(a+6,round(meanfunction(11,var)))
d8 <- rep(a+7,round(meanfunction(11,var)))
d9 <- rep(a+8,round(meanfunction(11,var)))
d10 <- rep(a+9,round(meanfunction(11,var)))
d11 <- rep(a+10*a,round(meanfunction(11,var)))

cell <- c(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11)
cell <- cell[1:50]

return(cell)

}

###function that generates a population of 5000 cells (by generating independent cell histories), calculates mean nuclear content and standard deviation of nuclear content at each time point (sd used as a measure of synchrony)
generatepopulation <- function(initcon=initcon,mean=mean,var=var) {

test2 <- replicate(5000,generatecellhistory(initcon,mean,var))

test3 <- as.data.frame(test2)
test3$mean <-  apply(test3, 1, function(x) {mean(x)})
test3$SD <- apply(test3, 1, function(x) {sd(x)})

return(test3$SD)

}



#####simulate 100 times for various variances of nuclear division times, subtract the SD from the t = 0, and calculate mean and standard deviation of the 100 independent simulations
###var = 0.3
time <- c(1:50)
pop <- replicate(100,generatepopulation(initcon,mean,var=0.3))

pop03 <- as.data.frame(pop)
pop03$mean <- apply(pop03,1,function(x) {mean(x)})
pop03$std <- apply(pop03,1,function(x) {sd(x)})
pop03 <- mutate(pop03, meannorm <- mean-mean[1])

interval03 <- pop03$meannorm + outer(pop03$std, c(1,-1))

###var = 0.1
pop <- replicate(100,generatepopulation(initcon,mean,var=0.1))

pop01 <- as.data.frame(pop)
pop01$mean <- apply(pop01,1,function(x) {mean(x)})
pop01$std <- apply(pop01,1,function(x) {sd(x)})
pop01 <- mutate(pop01, meannorm <- mean-mean[1])

interval01 <- pop01$meannorm + outer(pop01$std, c(1,-1))

###var = 0.5
pop <- replicate(100,generatepopulation(initcon,mean,var=0.5))

pop05 <- as.data.frame(pop)
pop05$mean <- apply(pop05,1,function(x) {mean(x)})
pop05$std <- apply(pop05,1,function(x) {sd(x)})
pop05 <- mutate(pop05, meannorm <- mean-mean[1])

interval05 <- pop05$meannorm + outer(pop05$std, c(1,-1))

###var = 0.2
pop <- replicate(100,generatepopulation(initcon,mean,var=0.2))

pop02 <- as.data.frame(pop)
pop02$mean <- apply(pop02,1,function(x) {mean(x)})
pop02$std <- apply(pop02,1,function(x) {sd(x)})
pop02 <- mutate(pop02, meannorm <- mean-mean[1])

interval02 <- pop02$meannorm + outer(pop02$std, c(1,-1))

###var = 0.05
pop <- replicate(100,generatepopulation(initcon,mean,var=0.05))

pop005 <- as.data.frame(pop)
pop005$mean <- apply(pop005,1,function(x) {mean(x)})
pop005$std <- apply(pop005,1,function(x) {sd(x)})
pop005 <- mutate(pop005, meannorm <- mean-mean[1])

interval005 <- pop005$meannorm + outer(pop005$std, c(1,-1))


### function that makes pretty transparent colors
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

#generate the palette of colours to plot the shaded regions of plots for various variances
library(wesanderson)

pal2 <- wes_palette("Zissou", 10, type="continuous")
pal <- makeTransparent(pal2, 120)



time3 <- c(1,13,25,37) #times correspnding to the experimental data; the first point in the simulation is actually t=0)

###plots the curves for different CVs
plot(time3-1,pop05$meannorm[time3],type="l",ylim=c(-0.1,0.5), xlab = "time (h)", ylab = "increase of geometric SD")
polygon(c(time3-1, rev(time3-1)), c(interval05[time3,1],rev(interval05[time3,2])), border = NA, col=pal[1])
polygon(c(time3-1, rev(time3-1)), c(interval03[time3,1],rev(interval03[time3,2])), border = NA, col=pal[2])
polygon(c(time3-1, rev(time3-1)), c(interval02[time3,1],rev(interval02[time3,2])), border = NA, col=pal[3])
polygon(c(time3-1, rev(time3-1)), c(interval01[time3,1],rev(interval01[time3,2])), border = NA, col=pal[4])
lines(time3-1, pop03$meannorm[time3], type = "l", col="black")
lines(time3-1, pop02$meannorm[time3], type = "l", col="black")
lines(time3-1, pop01$meannorm[time3], type = "l", col="black")
lines(time3-1, pop05$meannorm[time3], type = "l", col="black")

###adding real experimental data to the plot; values of geometric SD from the two biological replicates at corresponding times
data <- c(0.6626, 0.6807, 0.6268, 0.6903)
data2 <- c(0.797, 0.845, 0.820, 0.883)

datam <-data-data[1]
datam2 <- data2-data2[1]
points(time3-1,datam,type="b", lwd = 3)
points(time3-1,datam2,type="b", lwd = 3)

legend (0,0.5, c("CV = 0.5","CV = 0.3","CV = 0.2","CV = 0.1"), cex=1, col=c(pal[1],pal[2],pal[3],pal[4]),pch=15,bty="n",y.intersp = 1.5)







