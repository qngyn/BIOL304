//6.1
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
