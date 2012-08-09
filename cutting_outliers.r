#################### data
empdata <- read.table(file="emp_data.csv", header=TRUE, sep=";")

#################### pakaches
library('fitdistrplus')

#################### enter parametres
k <- 10000  					#объем генерируемой выборки
w <- 1000							#количество итераций
m <- c(2500, 3000, 3500) 		#верхние априорные пороги обрезки данных
n <- c(1, 1.5, 2, 2.5)				#множители для стандартного отклонения
p <- c(0.05)						#процент обрезки данных, поровну с каждого края

#################### variables
meanlog <-c(NULL)
sdlog <-c(NULL)
rw <- frame.data(NULL)
accTr <- c("2500", "3000", "3500")
accSD <- c("SD", "1,5SD", "2SD", "2,5SD")
accQ <- c("5%")


#################### fitting distribution
for (i in 1:3) {					#оцениваем распределение каждого массива данных времени реакции
  for (j in 1:6) {				#dtipes <-c("NO", "LOGNO", "GA", "WEI", "IG", "exGAUS")
    if (j==1) dstr <-"NO"
    if (j==2) dstr <-"LOGNO"
    if (j==3) dstr <-"GA"
    if (j==4) dstr <-"WEI"
    if (j==5) dstr <-"IG"
    if (j==6) dstr <-"exGAUS"
    aic <- gamlss(formula=empdata[,i] ~ 1, family=dstr)
    rw[i,j] <- aic$aic				#значение критtрия Акаике, сводная таблица
  }
}
colnames(rw) <-c("NO", "LOGNO", "GA", "WEI", "IG", "exGAUS")
rownames(rw) <-c("anagrams", "embedded", "mirrors")


#################### parameters estimating
for (i in 1:3) {							#оценка параметров для каждого массива данных 
  emp <- empdata[,i]						#удаляем NA из данных
  emp <- na.omit(emp)
  emp <- as.numeric(emp)
  f1 <- fitdist(emp,"lnorm",method="mle")	#подгонка параметров распределение, метод максимального правдоподобия
  b1 <- bootdist(f1)						#бустреп-оценка опараметров распределения
  c1 <- b1$CI
  meanlog[i] <- c1[1,1]
  sdlog[i] <- c1[2,1]
}
estMean <- mean(meanlog)					#средненее по трем эмпиричсеким выборкам 
estSD <- mean(sdlog)						#стандартное отклонение (для логнормального распределения) по выборке
spMean <- (estMean*1.1)						#среднее для моделирования распределения выбросов
spSD <- (estSD*1.1)							#стандартное отклонение для моделирования распределения выбросов

#################### accuracy of outliers cutoff strategies
for (i in 1:6) {						#процент выбросов в выборке, 5%-30%
  n1data <- k*(1-i*0.05) 				#генерируемое число сигналов (данных времени реакции)
  n2data <- k*i*0.05 					#генерируемое число выбросов (данных, не относящихся к времени реакции)
  f1 <- data.frame(rep(0,n1data))
  f2 <- data.frame(rep(1,n2data))
  names(f1) <- 0.05 -> names(f2)
  odds <- data.frame(rbind(f1, f2))	#группирующая переменная выбросов и сигналов
  ####Tresholds
  effTrt <- (NULL)
  for (l in m) {											#эффективность обрезки по априорным порогам
    effTr <- (NULL)
    for (j in 1:w) {
      g <- rlnorm(n1data, meanlog=estMean, sdlog=estSD) 		#генерируем "сигналы"
      spg <- rlnorm(n2data, meanlog=spMean, sdlog=spSD)		#генерируем "выбросы"
      gendataTr <- c(g,spg)
      gendataTr[gendataTr < 100] <- 1							#определяем границы обрезки, исходя из заданного априори порога обрезки данных
      gendataTr[gendataTr > l] <- 1
      gendataTr[gendataTr != 1] <- 0
      tbTr <- table(data.frame(odds, gendataTr))				#кросстабуляция группирующей и реальных данных
      NCorrOutlTr <- tbTr[4]										#здесь возможна ошибка
      NCorrSgnlTr <- tbTr[1]
      effTr <- data.frame(rbind(effTr, c(NCorrOutlTr, NCorrSgnlTr)))	#таблица точности по сгененированному массиву
    }
    perNSgnlTr <- mean(NCorrSgnlTr)/n1data
    perNOutlTr <- mean(NCorrOutlTr)/n2data
    effTrt <- data.frame(rbind(effTrt, c(perNSgnlTr, perNOutlTr)))		#общая таблица точности по методу обрезки
  }
  accTr <- data.frame(cbind(accTr, effTrt))							#общая таблица точности по разным стратегиям и с разной долей выбросов
  ####SD
  effSDt <- (NULL)
  for (l in n) {											#эффективность обрезки по стандартному отклонению
    effSD <- (NULL)
    for (j in 1:w) {
      g <- rlnorm(n1data, meanlog=estMean, sdlog=estSD) 		#генерируем "сигналы"
      spg <- rlnorm(n2data, meanlog=spMean, sdlog=spSD)		#генерируем "выбросы"
      gendataSD <- c(g,spg)
      genMx <- mean(gendataSD)+sd(gendataSD)					#задаем границы диапазона обрезки по стандартным отклонениям
      genMn <- mean(gendataSD)-sd(gendataSD)					#задаем границы диапазона обрезки по стандартным отклонениям
      genMn[genMn < 0] <- 0
      gendataSD[gendataSD < l*genMn] <- 1
      gendataSD[gendataSD > l*genMx] <- 1
      gendataSD[gendataSD != 1] <- 0
      tbSD <- table(data.frame(odds, gendataSD))				#кросстабуляция группирующей и реальных данных
      NCorrOutlSD <- tbSD[4]
      NCorrSgnlSD <- tbSD[1]
      effSD <- data.frame(rbind(effSD, c(NCorrOutlSD, NCorrSgnlSD)))	#таблица точности по сгененированному массиву
    }
    perNSgnlSD <- mean(NCorrSgnlSD)/n1data
    perNOutlSD <- mean(NCorrOutlSD)/n2data
    effSDt <- data.frame(rbind(effSDt, c(perNSgnlSD, perNOutlSD)))		#общая таблица точности по методу обрезки
  }
  accSD <- data.frame(cbind(accSD, effSDt))							#общая таблица точности по разным стратегиям и с разной долей выбросов
  ####Quantles
  effQt <- (NULL)
  for (l in p) {											#эффективность обрезки по квантилям
    effQ <- (NULL)
    for (j in 1:w) {
      g <- rlnorm(n1data, meanlog=estMean, sdlog=estSD)			#генерируем "сигналы"
      spg <- rlnorm(n2data, meanlog=spMean, sdlog=spSD)			#генерируем "выбросы"
      gendata <- c(g,spg)
      qdata <- gendata
      p975 <- quantile(qdata, .975)							#определяем квантили
      p025 <- quantile(qdata, .025)							#определяем квантили
      gendataQ <- gendata
      gendataQ[gendataQ > p975] <- 1							#маркируем данные, как принадлежащие обрезаемым квантилям
      gendataQ[gendataQ < p025] <- 1
      gendataQ[gendataQ != 1] <- 0
      tbQ <- table(data.frame(odds, gendataQ))					#кросстабуляция группирующей и реальных данных
      NCorrOutlQ <- tbQ[4]										#нередко возникает ошибка, вариант "лечения"  - ifelse (median(gendataQ)==0, NCorrOutl <- 0, tb[4]), дает подозрительные результаты
      NCorrSgnlQ <- tbQ[1]
      effQ <- data.frame(rbind(effQ, c(NCorrOutlQ, NCorrSgnlQ)))	#общая таблица точности по методу обрезки
    }
    perNSgnlQ <- mean(NCorrSgnlQ)/n1data
    perNOutlQ <- mean(NCorrOutlQ)/n2data
    effQt <- data.frame(rbind(effQt, c(perNSgnlQ, perNOutlQ))) 		#общая таблица точности по методу обрезки
  }
  accQ <- data.frame(cbind(accQ, effQt))								#общая таблица точности по разным стратегиям и с разной долей выбросов
}
colnames(accTr) <- c("Method", "CorrSgnl_5%", "CorrOut_5%", "CorrSgnl_10%", "CorrOut_10%", "CorrSgnl_15%", "CorrOut_15%","CorrSgnl_20%", "CorrOut_20%", "CorrSgnl_25%", "CorrOut_25%", "CorrSgnl_30%", "CorrOut_30%")
colnames(accSD) <- c("Method", "CorrSgnl_5%", "CorrOut_5%", "CorrSgnl_10%", "CorrOut_10%", "CorrSgnl_15%", "CorrOut_15%","CorrSgnl_20%", "CorrOut_20%", "CorrSgnl_25%", "CorrOut_25%", "CorrSgnl_30%", "CorrOut_30%")
colnames(accQ) <- c("Method", "CorrSgnl_5%", "CorrOut_5%", "CorrSgnl_10%", "CorrOut_10%", "CorrSgnl_15%", "CorrOut_15%","CorrSgnl_20%", "CorrOut_20%", "CorrSgnl_25%", "CorrOut_25%", "CorrSgnl_30%", "CorrOut_30%")
acc <- data.frame(rbind(accTr, accSD, accQ))		#финальная сводная таблица
