##6.1
data <- read.table("~/Desktop/lab6.1_data", header = T)
attach(data)
str(data)
avg <- mean(snow)
v <- var(snow)
NLL <- function(mu, k, x){
  -sum(dnbinom(x, size = k, mu = mu, log = TRUE))
}

kf <- function(vr, a){
  result = a/(vr/(a - 1))
  return (result)
}

k <- kf(v,avg)
snowmle <- mle2(NLL, start = list(mu = avg, k = k), data = list(x = data$snow))
kvec = seq(0.5, 1.5, by = 0.05)
muvec = seq(1, 3, by=0.05)
NLL.matrix = matrix(nrow=length(muvec), ncol=length(kvec))
for (i in 1:length(muvec))
{
for (j in 1:length(kvec))
{
NLL.matrix[i,j] = NLL(muvec[i], kvec[j], snow)
}
}
contour(muvec, kvec, NLL.matrix, xlab = "mu", ylab = "k")
points (2, 1, pch=2)

mu.mle = coef(snowmle)[1]
k.mle = coef(snowmle)[2]
snow.table = table (factor(snow, levels = 0:max(snow)))/length(snow)
b1 = barplot (snow.table, ylab="Frequency",
xlab="Number of snow days")
points (b1, dnbinom(x=0:max(snow), mu=mu.mle, size=k.mle), pch=16)

##
data <- read.csv("~/Desktop/ReefFish.txt", sep="")
attach(data)
##Function 1
NLL1 <- function(a){
  -sum(dbinom(recruits, settlers, prob = a, log = TRUE ))
}

##Function 2
NLL2 <- function(a, b){
  p = a/(1 + (a/b)*settlers)
  -sum(dbinom(x = recruits, settlers, prob = p, log = TRUE))
}


##Function 3
NLL3 <- function(a, b, d){
  p = a/(1 + ((a/b)*(settlers)^d))
  -sum(dbinom(x = recruits, settlers, prob = p, log = TRUE))
}

#Plot
plot(settlers, recruits, xlab = "settlers", ylab = "recruits")

recruits <- recruits[settlers > 0]
settlers <- settlers[settlers > 0]


#Optimized
f1 <-  mle2 (NLL1, start=list(a=0.5), method="L-BFGS-B", lower = 0.003)
f2 <- mle2(NLL2, start = list(a = 0.6, b = 10), method = "L-BFGS-B", lower = 0.003 )
f3 <- mle2 (NLL3, start=list(a=0.5, b=10, d=1), method="L-BFGS-B", 
  lower = c(a=0.003, b=1, d=0.003), upper=c(a=0.99, b=10, d=2))

#Plot
plot(settlers, recruits, xlab = "settlers", ylab = "recruits")
a1 = coef(f1)["a"]
x = seq(min(settlers), max(settlers), 1)
lines(x,(a1)*x, lty=3, lwd=2)

a.h = coef(ml.bh)["a"]
b.h = coef(ml.bh)["b"]
bev <- function(x){
  a.h*x/(1 + (a.h/b.h)*x)
}
lines(x, bev(x),lty=2, lwd=2)

a = coef(f3)["a"]
b = coef(f3)["b"]
d = coef(f3)["d"]
lines (x, a*x/(1 + (a/b)*x^d), lty=3, col=2, lwd=2)

nll = c(constant = -logLik(f1), BH = -logLik(f2), shep = -logLik(f3))
numparam = c(1, 2, 3)
aic = -2*nll + 2*numparam
AICtab(ml.const, ml.bh, ml.shep, weights = TRUE)

anova (ml.const, ml.bh, ml.shep)

confint(profile(f3, which="d"))
confint(profile(f3, which=c("a", "b")))
confint(profile(f3))

##
data <- read.table("~/Desktop/pistachio_data_allyears-1.txt", header=TRUE, quote="\"")
data <- na.omit(data)
attach(data)
library(bbmle)

sex <- as.numeric(as.factor(data$sex))

#Function 1
NLL.sepat <- function(a.m, a.f, t.m, t.f, data){
  a <- c(a.m, a.f)[sex]
  t <- c(t.m, t.f)[sex]
  n = number
  meanfound = a*n/(1 + a*t*n)
  -sum(dpois(captured, lambda = meanfound, log=T))
}

#Function 2
NLL.sepa <- function(a.m, a.f, t, data){
  a <- c(a.m, a.f)[sex]
  n = number
  meanfound = a*n/(1 + a*t*n)
  -sum(dpois(captured, lambda = meanfound, log=T))
}

#Funtion 3
NLL.sept <- function(a, t.m, t.f, data){
  t <- c(t.m, t.f)[sex]
  n = number
  meanfound = a*n/(1 + a*t*n)
  -sum(dpois(captured, lambda = meanfound, log=T))
}

#Function 4
NLL.none <- function(a, t, data){
  n = number
  meanfound = a*n/(1 + a*t*n)
  -sum(dpois(captured, lambda = meanfound, log=T))
}

#Optimize
sepat <- mle2(NLL.sepat, start = list(a.m = 0.875, a.f = 0.875, t.m = 0.125, t.f = 0.125), data = list (data = data))
sepa <- mle2(NLL.sepa, start = list(a.m = 0.875, a.f = 0.875, t = 0.125), data = list(data = data))
sept <- mle2(NLL.sept, start = list(a = 0.875, t.m = 0.125, t.f = 0.125), data = list(data = data))
nn <- mle2(NLL.none, start = list(a = 0.875, t = 0.125), method ="L-BFGS-B", lower=0.01)

#Plot
plot(data$number[data$sex == "F"], data$captured[data$sex == "F"], col = "pink", xlab = "Number available", ylab = "Number captured")
points(data$number[data$sex == "M"], data$captured[data$sex == "M"], col = "blue")

#Curving part
vec <- c(0:130)
a.m = coef(sepat)["a.m"]
a.f = coef(sepat)["a.f"]
t.m = coef(sepat)["t.m"]
t.f = coef(sepat)["t.f"]

func <- function(a,t,n){
  a*n/(1 + a*t*n)
}

lines(vec, func(a.m, t.m, vec), lty = 1, col = "black", lwd = 2)
lines(vec, func(a.f, t.f, vec), lty = 1, col = "red", lwd = 2)

a.m = coef(sepa)["a.m"]
a.f = coef(sepa)["a.f"]
t = coef(sepa)["t"]
lines(vec, func(a.m, t, vec), lty = 2, col = "black", lwd = 2)
lines(vec, func(a.f, t, vec), lty = 2, col = "red", lwd = 2)

a = coef(sept)["a"]
t.m = coef(sept)["t.m"]
t.f = coef(sept)["t.a"]
lines(vec, func(a, t.m, vec), lty = 3, col = "black", lwd = 3)
lines(vec, func(a, t.f, vec), lty = 3, col = "red", lwd = 3)

a = coef(nn)["a"]
t = coef(nn)["t"]
lines(vec, func(a,t,vec), lty = 4, col = "blue", lwd = 2)

legend (x="bottomright",
        legend=c("male", "female", "separate a and f", "separate a",
                 "separate tau", "no sex diff."),
        col=c("black", "red", "black", "black", "black", "blue"),
        lty=c(1,1,1,2,3,1), bty="n")


aic = AICtab(sepat, sepa, sept, nn, weights = TRUE)
aic
anova(sepat, sepa, nn)
anova(sepat, sept, nn)

sepat.ll <- logLik(sepat)
sepat.prof <- profile(sepat)
plot(sepat.prof)
confint(sepat.prof)

none.ll <- logLik(nn)
nn.ll <- profile(nn)
plot(nn.ll)
confint(nn)

##
data <- read.csv("~/Desktop/quiz6.3Data-1.txt", sep="")
attach(data)
library(bbmle)

##2
t.low <- time[light == "low"]
t.high <- time[light == "high"]
xmin <- min(time)
xmax <- max(time)
plot(t.low, 1/t.low, xlab = "Time to explosion", ylab ="Explosion Rate", xlim = c(xmin, xmax), col = "black")
points(t.high, 1/t.high, col = "red" )
legend("topright", legend = c("low light", "high light"), col = c("black", "red"), pch = 1)

summary(lm(1/t.low ~ t.low))
summary(lm(1/t.high ~ t.high))

a.high = 0.31
a.low = 0.19
b = - 0.004

NLL.none <- function(a, b, data){
  shape = var(data$time)/(mean(data$time)^2)
  ep = a + b*data$dbh
  -sum(dgamma(x = data$time, rate = ep, shape = shape, log = T))
}

nn <- mle2(NLL.none, start = list(a = 0.22, b = - 0.04), data = list(data = data))
a.n <- coef(nn)[1]
b.n <- coef(nn)[2]

plot(dbh, 1/time, xlab = "breast height", ylab = "explosion rate")
lines(dbh, a.n + b.n*dbh, lty = 1, col = "blue", lwd = 1)

NLL.sepa <- function(a.h, a.l, b, data){
  shape = var(data$time)/(mean(data$time)^2)
  lv = as.numeric(as.factor(data$light))
  a = c(a.h, a.l)[lv]
  ep = a + b*data$dbh
  -sum(dgamma(x = data$time, rate = ep, shape = shape, log = T))
}

sepa <- mle2(NLL.sepa, start = list(a.h = 0.31, a.l = 0.19, b = - 0.04), data = list(data = data))
a.l <- coef(sepa)["a.l"]
a.h <- coef(sepa)["a.h"]
b <- coef(sepa)["b"]
lines(dbh, a.h + b*dbh, lty = 1, col = "red", lwd = 1)
lines(dbh, a.l + b*dbh, lty = 1, col = "black")
legend("topright", legend = c("no diff a", "high light", "low light"), col = c("blue", "red", "black"), lty = 1)

aic = AICtab(nn, sepa)
aic
anova(nn, sepa)

##
library (bbmle)
set.seed(1001)
mu.true = 1
k.true = 0.4
d = rnbinom (50, mu = mu.true, size = k.true)
NLL.nbinom <- function (mu, k, dat) {
-sum(dnbinom(dat, mu = mu, size = k, log=T))
}

ml.mle = mle2(minuslogl = NLL.nbinom, start = list(mu = 1.2 * mu.true,
k = 1.2*k.true), data = list(dat = d))

coef(ml.mle)
mu.est = coef(ml.mle)["mu"]
k.est <- coef(ml.mle)["k"]
ll <- logLik(ml.mle)
ml.mle.prof = profile(ml.mle)
plot(ml.mle.prof)
confint(ml.mle.prof)

##
##Using negative binomial 
d <- read.csv("~/Desktop/quiz6.1Data-1.txt", sep="")
foo = tapply (d$numYes, d$weeks, mean)
plot (sort(unique(d$weeks)), foo, type="b")

## Remember that the mean number of sucesses will be p*10.  Dividing
## those numbers on the graph by 10, it looks to me like a = 0.1, b = (3.6-1)/(10*10) = 0.026

## NLL function
NLL = function (a, b, data) {
  p = a + b*data$weeks
  -sum(dbinom(data$numYes, prob=p, size=10))
}

## Do the estimate
NLL.mle = mle2 (minuslogl=NLL, start=list(a=0.1, b=0.026), data=list(data=d))

a.est = coef(NLL.mle)[1]
b.est = coef(NLL.mle)[2]

weekSeq = seq(1, 10, 0.01)
plot (d$weeks, jitter(d$numYes), xlab="Weeks of lessons",
      ylab="Number who said yes", cex.lab=1.4)
lines (weekSeq, 10*(a.est + b.est*weekSeq), col="red")

##
dat <- read.csv("~/Desktop/HvZdata_fall2011_nostarve.csv")
attach(dat)

stepfunc <- function(data){
  n = length(data$zombie)
  I = data$zombie[-n]
  S = data$human[-n]
  y.pred = I*S*exp(b*S*dt)/(S-I+I*exp(b*S*dt))
  -sum(dnorm(data$zombie[-1], y.pred, sd = sigma, log = T))
}

####Using negative binomial
data <- read.csv("~/Desktop/quiz6.2Data-1.txt", sep="")
attach(data)
plot(duration, after/before, xlab = "duration", ylab = "fraction of people who left")

summary(lm(after/before ~duration))


NLL <- function(a,b, data){
  p = a + b*data$duration
  -sum(dnbinom(x = data$duration , size = data$before  , prob = p, log = T))
}

mle <- mle2(NLL, start = list(a = 0.819, b = -0.017), data = list(data = data))

ll <- logLik(mle)
mle.prof = profile(mle)
plot(mle.prof)
confint(mle.prof)

plot(duration, after/before, xlab = "duration", ylab = "fraction of people who left")
lines(duration, - 0.017*duration + 0.819)
lines(duration, coef(mle)[1] + coef(mle)[2]*duration, col = "blue")

##
data <- read.csv("~/Desktop/quiz11Data-2.txt", sep="")
attach(data)
NLL <- function(a, b, e, dt, data) {
  l = length(data$pop)
  y2 = data$pop[-l]
  y.pred = y2 + (a*data$pop[-l]- b*data$pop[-l]^(1+e))*dt
  -sum(dpois(x = data$pop[-1], y.pred, log = T))
}

mle1 <- mle2(NLL, start = list(a = 1.8, b = 1, e = 0.1), data = list(data = data, dt = 1))

plot(time, pop, type = "b", col = "red")

pred <- function(a,b,e,dt){
  l = length(data$pop)
  y2 = data$pop[-l]
  y.pred = y2 + (a*data$pop[-l]- b*data$pop[-l]^(1+e))*dt
}

a.p <- coef(mle1)["a"]
b.p <- coef(mle1)["b"]
e.p <- coef(mle1)["e"]
points(pred(1.8, 1, 0.1, 1), col = "blue", type = "b")
points(pred(a.p, b.p, e.p, 1), col = "pink", type = "b")
legend("bottomright", legend = c ("population data", "predicited value", "optimized value"), col = c("red", "blue", "pink"), pch = 1)

