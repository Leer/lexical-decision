empdata <- read.table(file="emp_data.csv", header=TRUE, sep=";")
library(gamlss)
library(boot)
library(multicore)

dists <- c("NO", "LOGNO", "GA", "WEI", "IG", "exGAUS")

## Функция расчитывает AIC для разных типов распредления
## Возвращает таблицу: Распределение - AIC
define.dist <- function(x) {
  rw <- NULL
  x <- as.numeric(na.omit(x))
  for (i in dists) {
    aic <- gamlss(formula = x ~ 1, family = i)$aic
    aic <- round(aic, digits = 2)
    rw <- c(rw, aic)
  }
  return(rw)
}

rw <- data.frame(mclapply(empdata, FUN = define.dist, mc.silent = TRUE))
rownames(rw) <- dists
# Печать результатов


# Наиболее подходящий тип распределения
dist <- dists[(sapply(rw, FUN = which.min))]
#print(dist)

## Функция для бутстрапа, расчитывает параметры заданного распределения (mu, sigma)
fit.dist <- function(data, indices, dist) {rw
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
  boot.out <- boot(x, statistic=fit.dist, R=1000, parallel = "multicore", dist = dist)
  mu <- boot.out$t0[1]
  sigma <- boot.out$t0[2]
  conf.mu <- boot.ci(boot.out, index = 1, type = "bca")$bca[4:5]
  conf.sigma <- boot.ci(boot.out, index = 2, type = "bca")$bca[4:5]
  result <- c(mu, conf.mu, sigma, conf.sigma)
  return((result))
}

distparams  <- data.frame(mclapply(empdata[,1:3], FUN = boot.dist, dist = LOGNO, mc.silent = TRUE))
rownames(distparams) <- c("mu", "mu - 95% low", "mu - 95 up", "sigma",  "sigma - 95% low", "sigma - 95 up")

estMean <- rowMeans(distparams)[1]
estSD <- rowMeans(distparams)[3]
# Печать результатов
print(rw)
print(distparams)