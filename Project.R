library(bbmle)
library(emdbook)
data <- read.csv("~/Desktop/RsylDDForRobinFixed.txt", sep="")
attach(data)

#Plot 
plot(EggDensity[CanopyBinary == 0],LarvalSurvival[CanopyBinary == 0], xlab = "Egg Density", ylab = "Larval Survival")
points(EggDensity[CanopyBinary ==  1], LarvalSurvival[CanopyBinary == 1], col = "red")
legend(x = "topright", col = c("black","red"), pch = c(1,1), legend = c("Close", "Open"), border = F, bty = "n")
egg= 0:350

#Rough Estimate 
lines(egg, 0.8/(1+(egg*0.6)), col="orange")
lines(egg, 0.4/(1+(egg*0.0015)), col= "yellow")
legend (x  = 190, y = 1.05,  col = c("orange", "yellow"), lty = c(1,1), legend = c("Rough Estimate 1", "Rough Estimate 2"), border = F, bty = "n")


#Models
NLL.one <- function(S, theta, data ){
   survival = S
  -sum(dbetabinom(x = data$metamorphs, prob = survival, size = data$InitialEggs, theta = theta, log = TRUE))
}

mle.one <- mle2(NLL.one, start = list(S = 0.4, theta = 1.5), data = list(data = data), method = "BFGS")
mle.one

plot(EggDensity[CanopyBinary == 0],LarvalSurvival[CanopyBinary == 0], xlab = "Egg Density", ylab = "Larval Survival")
points(EggDensity[CanopyBinary ==  1], LarvalSurvival[CanopyBinary == 1], col = "red")
legend(x = "topright", col = c("black","red"), pch = c(1,1), legend = c("Close", "Open"), border = F, bty = "n")

s.one = coef(mle.one)["S"]
abline(s.one, 0, col = "green")
prof <- profile(mle.one, method = "BFGS")
confint(prof)

NLL.two <- function(S, a, theta, data){
  survival= S/(1+(a*(data$EggDensity)))
  -sum(dbetabinom(x = data$metamorphs, prob = survival, size = data$InitialEggs, theta = theta, log = TRUE))
}

mle.two <- mle2(NLL.two, start=list(S = 0.4, a = 0.0015, theta = 1.5), data=list(data = data), method = "Nelder-Mead")

s.two = coef(mle.two)["S"]
a.two = coef(mle.two)["a"]

plot(EggDensity[CanopyBinary == 0],LarvalSurvival[CanopyBinary == 0], xlab = "Egg Density", ylab = "Larval Survival")
points(EggDensity[CanopyBinary ==  1], LarvalSurvival[CanopyBinary == 1], col = "red")
legend(x = 260, y = 1.05, col = c("black","red"), pch = c(1,1), legend = c("Close", "Open"), border = F, bty = "n")
lines(egg, s.two/(1+(egg*a.two)), col= "light blue")

legend(x = 250 ,y =  0.9 , col = "light blue", lty = 1, legend = "With So and a", border = F, bty = "n")

#Models Testing if Canopy Matters
plot(EggDensity[CanopyBinary == 0], LarvalSurvival[CanopyBinary == 0], col = "black", xlab = "Egg Density", ylab ="Larval Survival")
points(EggDensity[CanopyBinary == 1], LarvalSurvival[CanopyBinary == 1], col= "red")
legend(x = 245, y = 1.05, col = c("black","red"), pch = c(1,1), legend = c("Close", "Open"), bty = "n", border = F)

#s depending on canopy without a
NLL.three <- function(s.close, s.open, theta, data){
  code = as.factor(data$CanopyBinary)
  S = c(s.close, s.open)[code]
  survival = S
  -sum(dbetabinom(x = data$metamorphs, prob = survival, size = data$InitialEggs, theta = theta, log = TRUE))
}

mle.three <- mle2(NLL.three, start = list(s.open = 0.4, s.close = 0.4,theta = 1.5), data = list(data = data), method = "Nelder-Mead")
mle.three

s.close3 =  coef(mle.three)["s.close"]
s.open3 = coef(mle.three)["s.open"]
abline(s.close3, 0, col = "purple")
abline(s.open3, 0, col = "purple", lty= 2)
legend(x = 235, y = 0.9, col = "purple", lty = c(1,2), legend = c("So of close canopy", "So of open canopy"), border = F, bty ="n")

#a depends on canopy not S
plot(EggDensity[CanopyBinary == 0], LarvalSurvival[CanopyBinary == 0], col = "black", xlab = "Egg Density", ylab ="Larval Survival")
points(EggDensity[CanopyBinary == 1], LarvalSurvival[CanopyBinary == 1], col = "red")
legend(x = 250, y = 1.05, col = c("black","red"), pch = 1, legend= c("Close","Open"), bty = "n", border = F)

NLL.four <-  function(a.close, a.open, S, theta, data){
  code = as.factor(data$CanopyBinary)
  a = c(a.close, a.open)[code]
  survival = S/(1+(a*(data$EggDensity)))
-sum(dbetabinom(x = data$metamorphs, prob = survival, size = data$InitialEggs, theta = theta, log = TRUE))
}

mle.four <- mle2(NLL.four, start = list(S = 0.4, a.close = 0.0015, a.open = 0.0015, theta = 1.5), data = list(data = data), method = "Nelder-Mead")
mle.four

s.four = coef(mle.four)["S"]
a.close4 = coef(mle.four)["a.close"]
a.open4 = coef(mle.four)["a.open"]

lines(egg, s.four/(1+(egg*a.close4)), col = "orange")
lines(egg, s.four/(1+(egg*a.open4)), col = "orange", lty = 2)
legend(x = 240, y = 0.9, col = "black", lty = c(1,2), legend = c("a of close canopy", "a of open canopy"), border = F, bty ="n")


#S depends on canopy not a
plot(EggDensity[CanopyBinary == 0], LarvalSurvival[CanopyBinary == 0], col = "black", xlab = "Egg Density", ylab = "Larval Survival")
points(data$EggDensity[data$CanopyBinary == 1], data$LarvalSurvival[data$CanopyBinary == 1], col = "red")
legend(x = 245, y = 1.05, col = c("black","red"), pch = 1, legend = c("Close","Open"), bty = "n", border = F)

NLL.five <- function(s.open, s.close, a, theta, data){
code = as.factor(data$CanopyBinary)
  S = c(s.close, s.open)[code]
  survival= S/(1+(a*(data$EggDensity)))
  -sum(dbetabinom(x= data$metamorphs, prob = survival, size = data$InitialEggs, theta = theta, log = TRUE))
}

mle.five <- mle2(NLL.five, start = list(s.open = 0.4, s.close = 0.4, a = 0.0015, theta = 1.5), data = list(data = data), method = "Nelder-Mead")
mle.five

s.close5 = coef(mle.five)["s.close"]
s.open5 = coef(mle.five)["s.open"]
a.5= coef(mle.five)["a"]
lines(egg, s.close5/(1+(egg*a.5)), col = "pink")
lines(egg, s.open5/(1+(egg*a.5)), col = "pink", lty = 2)
legend(x = 235, y = 0.9, col = "pink", lty = c(1,2), legend = c("So of close canopy", "So of open canopy"), border = F, bty ="n")

#S and a depend on canopy

plot(EggDensity[CanopyBinary == 0], LarvalSurvival[CanopyBinary == 0], col = "black", xlab = "Egg Density", ylab = "Larval Survival")
points(EggDensity[CanopyBinary == 1], LarvalSurvival[CanopyBinary == 1], col = "red")
legend(x = 230, y = 1.05, col = c("black","red"), pch = 1, legend = c("Close","Open"), bty = "n", border = F)

NLL.six <- function(s.open, s.close, a.open, a.close, theta, data){
 code = as.factor(data$CanopyBinary)
 S = c(s.close, s.open)[code]
 a = c(a.close, a.open)[code]
 survival= S/(1+(a*(data$EggDensity)))
-sum(dbetabinom(x = data$metamorphs, prob = survival, size = data$InitialEggs, theta = theta, log = TRUE))
}

mle.six= mle2(NLL.six, start = list(s.open = 0.4, s.close = 0.4, a.open = 0.0015,a.close = 0.0015, theta = 1.5), data = list(data = data), method = "Nelder-Mead")
mle.six

s.open6 = coef(mle.six)["s.open"]
a.open6 = coef(mle.six)["a.open"]
s.close6 = coef(mle.six)["s.close"]
a.close6 = coef(mle.six)["a.close"]
lines(egg, s.close6/(1+(egg*a.close6)), col = "blue")
lines(egg, s.open6/(1+(egg*a.open6)), col = "blue", lty = 2)

legend(x = 220 , y = 0.9, col = "blue", lty = c(1,2), legend = c("So, a of close canopy","So, a of open canopy"), bty = "n", border = F)

plot(EggDensity[CanopyBinary == 0], LarvalSurvival[CanopyBinary == 0], col = "black", xlab = "Egg Density", ylab = "Larval Survival")
points(EggDensity[CanopyBinary == 1], LarvalSurvival[CanopyBinary == 1], col = "red")
##1
abline(s.one, 0, col = "green")

##2
lines(egg, s.two/(1+(egg*a.two)), col= "light blue")

##3
abline(s.close3, 0, col = "purple")
abline(s.open3, 0, col = "purple", lty= 2)

##4
lines(egg, s.four/(1+(egg*a.close4)), col = "black")
lines(egg, s.four/(1+(egg*a.open4)), col = "black", lty = 2)

##5
lines(egg, s.close5/(1+(egg*a.5)), col = "pink")
lines(egg, s.open5/(1+(egg*a.5)), col = "pink", lty = 2)

##6
lines(egg, s.close6/(1+(egg*a.close6)), col = "blue")
lines(egg, s.open6/(1+(egg*a.open6)), col = "blue", lty = 2)

legend(x = "topright", col = c("light blue", "green", "purple", "orange", "pink", "blue"), lty = 1, legend = c("So only", "a and So", "seperate So w/o a", "seperate a", "seperate So with a", "seperate a and So"), bty = "n", border = F)
legend (x = 193, y = 1.03, col = c("black", "red"), legend = c("Close", "Red"), bty = "n", border = F, pch = 1)

anova(mle.one, mle.three, mle.five, mle.six)
anova(mle.one, mle.two, mle.four, mle.six)

AICtab(mle.one, mle.two, mle.three, mle.four, mle.five, mle.six, weights = TRUE)

attach(data)
par(oma = c(4, 1, 1, 1))
par(mfrow=c(2,3))

plot(EggDensity[CanopyBinary == 0], LarvalSurvival[CanopyBinary == 0], col = "black", xlab = "Egg Density", ylab = "Larval Survival", main = "a. Model 1: with So only", cex.main = 0.8)
points(EggDensity[CanopyBinary == 1], LarvalSurvival[CanopyBinary == 1], col = "red")
abline(s.one, 0, col = "green", lty = 3)

plot(EggDensity[CanopyBinary == 0], LarvalSurvival[CanopyBinary == 0], col = "black", xlab = "Egg Density", ylab = "Larval Survival", main = "b. Model 2: With So and a", cex.main = 0.8)
points(EggDensity[CanopyBinary == 1], LarvalSurvival[CanopyBinary == 1], col = "red")
lines(egg, s.two/(1+(egg*a.two)), col= "light blue", lty =)

plot(EggDensity[CanopyBinary == 0], LarvalSurvival[CanopyBinary == 0], col = "black", xlab = "Egg Density", ylab = "Larval Survival", main = "c. Model 3: a-independent with seperate So", cex.main = 0.8)
points(EggDensity[CanopyBinary == 1], LarvalSurvival[CanopyBinary == 1], col = "red")
abline(s.close3, 0, col = "purple")
abline(s.open3, 0, col = "purple", lty= 2)

plot(EggDensity[CanopyBinary == 0], LarvalSurvival[CanopyBinary == 0], col = "black", xlab = "Egg Density", ylab = "Larval Survival", main ="d. Model 4: with seperate a", cex.main = 0.8)
points(EggDensity[CanopyBinary == 1], LarvalSurvival[CanopyBinary == 1], col = "red")
lines(egg, s.four/(1+(egg*a.close4)), col = "orange")
lines(egg, s.four/(1+(egg*a.open4)), col = "orange", lty = 2)

plot(EggDensity[CanopyBinary == 0], LarvalSurvival[CanopyBinary == 0], col = "black", xlab = "Egg Density", ylab = "Larval Survival", main = "e. Model 5: a-dependent with seperate So", cex.main = 0.8)
points(EggDensity[CanopyBinary == 1], LarvalSurvival[CanopyBinary == 1], col = "red")
lines(egg, s.close5/(1+(egg*a.5)), col = "pink")
lines(egg, s.open5/(1+(egg*a.5)), col = "pink", lty = 2)

plot(EggDensity[CanopyBinary == 0], LarvalSurvival[CanopyBinary == 0], col = "black", xlab = "Egg Density", ylab = "Larval Survival", main = "f. Model 6: with seperate So and a", cex.main = 0.8)
points(EggDensity[CanopyBinary == 1], LarvalSurvival[CanopyBinary == 1], col = "red")
lines(egg, s.close6/(1+(egg*a.close6)), col = "blue")
lines(egg, s.open6/(1+(egg*a.open6)), col = "blue", lty = 2)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottomright", legend = c("Close", "Open"), xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n", pch = c(1,1), col = c("black", "red"), cex = 1)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottomleft", legend = c("Close", "Open", "Canopy-independent"), xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty = "n", lty = c(1,2,3), col = "black", cex = 1)
