empdata <- read.table(file="emp_data.csv", header=TRUE, sep=";")
library(gamlss)
library(boot)
library(multicore)

dists <- c("NO", "LOGNO", "GA", "WEI", "IG", "exGAUS")

## Функция расчитывает AIC для разных типов распредления
## Возвращает AIC для каждого типа распределения
define.dist <- function(x, distname = FALSE) {
  rw <- NULL
  x <- as.numeric(na.omit(x))
  for (i in dists) {
    aic <- gamlss(formula = x ~ 1, family = i)$aic
    aic <- round(aic, digits = 2)
    rw <- c(rw, aic)
  }
  if (distname == TRUE) {
    distname <- dists[which.min(rw)]
    return(distname)
  }
  else {
    rw <- matrix(rw)
    rownames(rw) <- dists
    return(rw)
  }
}

## Функция для бутстрапа, расчитывает параметры заданного распределения (mu, sigma)
fit.dist <- function(data, indices, dist) {
  x <- as.numeric(na.omit(data[indices]))
  fit.out <- gamlss(x ~ 1, family = dist)
  mu <- as.numeric(fit.out$mu.coefficients)
  sigma <- abs(as.numeric(fit.out$sigma.coefficients))
  result <- c(mu, sigma)
  return(result)
}

## Функция бустрапит параметры распределения
## Взвращает вектор: mu, ниж. 95%, верх. 95%, сигма, ниж. 95%, верх. 95%
boot.dist <- function(x, dist) {
# Задаём параметры для геренатора (для воспроизводимости результатов)
  set.seed(1234)
# Расчитываем параметры распределения
  boot.out <- boot(x, statistic=fit.dist, dist = dist, R=1000, parallel = "multicore")
  mu <- boot.out$t0[1]
  sigma <- boot.out$t0[2]
  conf.mu <- boot.ci(boot.out, index = 1, type = "bca")$bca[4:5]
  conf.sigma <- boot.ci(boot.out, index = 2, type = "bca")$bca[4:5]
  result <- c(mu, conf.mu, sigma, conf.sigma)
  return(result)
}

# Таблица значений AIC для каждого типа распределения
dist.aic <- data.frame(mclapply(empdata[,1:3], FUN = define.dist, mc.silent = TRUE))
rownames(dist.aic) <- dists

# Таблица параметров распределния + доверительные интервалы
dist.params  <- data.frame(mclapply(empdata[,1:3], FUN = function(x) boot.dist(x, dist = define.dist(x, distname = TRUE)),  mc.silent = TRUE))
rownames(dist.params) <- c("mu", "mu - 95% low", "mu - 95 up", "sigma",  "sigma - 95% low", "sigma - 95 up")

estMean <- rowMeans(dist.params)[1]
estSD <- rowMeans(dist.params)[3]
# Печать результатов
print(dist.aic)
print(dist.params)