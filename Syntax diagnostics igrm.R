##Install Packages
library(MASS)
library(goftest)
library(actuar)
library(fitdistrplus)
library(ggplot2)
library(statmod)
library(readxl)

DataPDRB <- read_excel("C:/Users/user/Documents/SKRIPSIII/DATA/DATA PENELITIAN.xlsx")
View(DataPDRB)
y <- DataPDRB$Y

hist(DataPDRB$Y, 
     breaks = 10,                       
     col = "beige",                    
     main = "Histogram PDRB",      
     xlab = "Nilai PDRB",                       
     border = "black")

# Estimasi parameter distribusi inverse gaussian
mu <- mean(y)
sigma2 <- var(y)

# fungsi distribusi kumulatif untuk Inverse Gaussian
inverse_gaussian_cdf <- function(y, mu, sigma2) {
  
  term <- sqrt(1 / (sigma2 * y))
  part1 <- pnorm(term * (y / mu - 1))
  part2 <- exp(2 / (sigma2 * mu)) * pnorm(-term * (y / mu + 1))
  
  F_y <- part1 + part2
  
  return(F_y)
}

# Uji Kolmogorov-Smirnov
ks_test <- ks.test(y, function(x) inverse_gaussian_cdf(x, mu, sigma2))

mydata <- DataPDRB[,]
View(mydata)
summary(mydata)

y <- mydata$Y
X1 <- mydata$X1
X2 <- mydata$X2
X3 <- mydata$X3
X4 <- mydata$X4
X5 <- mydata$X5

# Scatter plot
# VISUAL PLOT X1-X5 VS 1/(Y^2) Data Full
plot(mydata$X1,1/(mydata$Y^2))
plot(mydata$X2,1/(mydata$Y^2))
plot(mydata$X3,1/(mydata$Y^2))
plot(mydata$X4,1/(mydata$Y^2))
plot(mydata$X5,1/(mydata$Y^2))

warna_titik <- "navy"
p1 <- ggplot(mydata, aes(x = X1, y = 1/(Y^2))) +
  geom_point(color = warna_titik, size = 2) +
  labs(title = "X1 vs 1/Y²", x = "X1", y = "1/Y²") +
  theme_minimal()
p2 <- ggplot(mydata, aes(x = X2, y = 1/(Y^2))) +
  geom_point(color = warna_titik, size = 2) +
  labs(title = "X2 vs 1/Y²", x = "X2", y = "1/Y²") +
  theme_minimal()
p3 <- ggplot(mydata, aes(x = X3, y = 1/(Y^2))) +
  geom_point(color = warna_titik, size = 2) +
  labs(title = "X3 vs 1/Y²", x = "X3", y = "1/Y²") +
  theme_minimal()
p4 <- ggplot(mydata, aes(x = X4, y = 1/(Y^2))) +
  geom_point(color = warna_titik, size = 2) +
  labs(title = "X4 vs 1/Y²", x = "X4", y = "1/Y²") +
  theme_minimal()
p5 <- ggplot(mydata, aes(x = X5, y = 1/(Y^2))) +
  geom_point(color = warna_titik, size = 2) +
  labs(title = "X5 vs 1/Y²", x = "X5", y = "1/Y²") +
  theme_minimal()
grid.arrange(p1, p2, p3, p4, p5, ncol = 3)

# NEWTON RAPHSON
X <- model.matrix(~ X1 + X2 + X3 + X4 + X5, data = mydata)

igrm <- function(X, y, max_iter = 100, tol = 1e-6) {
  n <- length(y)
  p <- ncol(X)
  
  # Inisialisasi beta0
  fit_init <- glm(y ~ X - 1, family = inverse.gaussian)
  betaold <- coef(fit_init)
  for (r in 1:max_iter) {
    eta_hat <- X %*% betaold                  
    mu_hat <- 1 / sqrt(eta_hat)                
    if (any(mu_hat <= 0)) stop("Nilai mu <= 0, periksa X atau 
  inisialisasi.")
    W <- diag(as.vector(mu_hat^3))             
    z <- eta_hat + (y - mu_hat) / (mu_hat^3)   
    XtW <- t(X) %*% W
    betanew <- solve(XtW %*% X) %*% XtW %*% z
    if (sum(abs(betanew - betaold)) < tol) {
      cat("Konvergen pada iterasi ke-", r, "\n")
      break
    }
    betaold <- betanew
  }
  
  return(list(beta = betanew, iterasi = r))
}
hasil <- igrm(X, y)
hasil$beta

#Model Data Full
modelIG <- glm(Y ~ X1 + X2 + X3 + X4 + X5, 
               family = inverse.gaussian(link = "1/mu^2"), 
               data = mydata)
summary(modelIG)

#Menghitung nilai prediksi mu pengamatan ke-1
(modelIG$fitted.values)[1]
(modelIG$residuals)[1]
predict(modelIG, mydata[1,2:6])
x1 <- mydata[1,2:6]
x2 <- x1+c(1,0,0,0,0)

eta1 <- predict(modelIG, x1)
eta2 <- predict(modelIG, x2)

mu1 <- 1/sqrt(eta1)
mu2 <- 1/sqrt(eta2)

# DIAGNOSTIK PENGAMATAN BERPENGARUH DATA FULL
# LEVERAGE
hii <- hatvalues(modelIG)

# Cutoff
n <- nrow(X)
p <- ncol(X)-1
cutoff_leverage <- 2 * (p + 1) / n
hii_df <- data.frame(Wilayah=mydata$Wilayah,
                     Observasi = 1:nrow(mydata),  
                     Y = mydata$Y,
                     Cook_Distance = hii,
                     Cutoff = cutoff_leverage,
                     Influential = hii> cutoff_leverage
)

hii_df[which(hii> cutoff_leverage), ]

##BUANG INDRAMAYU, PURWAKARTA, PANGANDARAN, KOTA SUKABUMI, DAN KOTA BANJAR

mydata<-DataPDRB[-c(12, 14,18, 20, 27),]
modelIG <- glm(Y ~ X1+X2+X4+X5+X3,
               family = inverse.gaussian(link = "1/mu^2"), 
               data = mydata)
summary(modelIG)

# DIAGNOSTIK PENGAMATAN BERPENGARUH 22 PENGAMATAN
# LEVERAGE
hii <- hatvalues(modelIG)

# Cutoff
n <- nrow(X)
p <- ncol(X)-1
cutoff_leverage <- 2 * (p + 1) / n
hii_df <- data.frame(Wilayah=mydata$Wilayah,
                     Observasi = 1:nrow(mydata),  
                     Y = mydata$Y,
                     Cook_Distance = hii,
                     Cutoff = cutoff_leverage,
                     Influential = hii> cutoff_leverage
)

hii_df[which(hii> cutoff_leverage), ]

# PEARSON RESIDUAL
mu_hat<-modelIG$fitted.values 
chi_i <- (modelIG$y - mu_hat) / sqrt(1.200674e-06*mu_hat^3)
pearson_df <- data.frame(Wilayah=mydata$Wilayah,
                         Observasi = 1:nrow(mydata),  
                         Y = mydata$Y,
                         mu_hat = mu_hat,
                         Chi_i = chi_i,
                         Outlier = abs(chi_i) > 2
)
print(pearson_df)
pearson_df[abs(chi_i) > 2, ]

# COOK'S DISTANCE
n <- length(y)
# Estimasi parameter dispersi (phi)
phi_hat <- sum((y - mu_hat)^2 / mu_hat^3) / (n - p - 1)

# Pearson residual terstandarisasi
chi_i_std <- chi_i / sqrt(phi_hat * (1 - hii))

cook_D <- (chi_i_std^2 / p + 1) * (hii / (1 - hii))
cutoff_cook <- 4 / (n - 1)

cook_df <- data.frame(Wilayah=mydata$Wilayah,
                      Observasi = 1:nrow(mydata),  
                      Y = mydata$Y,
                      Cook_Distance = cook_D,
                      Cutoff = cutoff_cook,
                      Influential = cook_D > cutoff_cook
)

print(cook_df)
cook_df[which(cook_D > cutoff_cook), ]

# COVARIANCE RATIO
p_prime <- p + 1  # Jumlah parameter termasuk intersep
CVR_i <- ((n - p - chi_i_std^2) / (n - p - 1))^p_prime / (1 - hii)

cutoff_CVR_upper <- 1 + (3 * (p + 1)) / n
cutoff_CVR_lower <- 1 - (3 * (p + 1)) / n
cvr_df <- data.frame(Wilayah=mydata$Wilayah,
                     Observasi = 1:nrow(mydata),  
                     Y = mydata$Y,
                     CVR = CVR_i,
                     Influential = (CVR_i > cutoff_CVR_upper) | 
                       (CVR_i    < cutoff_CVR_lower)
) 

print(cvr_df)
cvr_df[which(CVR_i > cutoff_CVR_upper| CVR_i < cutoff_CVR_lower), ]

# WELSCH'S DISTANCE
t_i <- chi_i_std * sqrt((n - p - 1) / (n - p - chi_i_std^2))

WD_i <- (n - 1) * t_i^2 * hii / (1 - hii)^2
cutoff_WD <- 3 * sqrt(p + 1)

wd_df <- data.frame(Wilayah=mydata$Wilayah,
                    Observasi = 1:nrow(mydata),  
                    Y = mydata$Y,
                    WD = WD_i,
                    Influential = WD_i > cutoff_WD
)

print(wd_df)
wd_df[which(WD_i > cutoff_WD), ]
