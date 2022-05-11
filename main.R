suppressMessages({
  library(spatstat)
  library(data.table)
  library(dbmss)})
load("./CSBIGS.Rdata")
setDT(Emergencies)
nrow(unique(Emergencies[, c(1, 2)])) # unique number of points
Emer <- Emergencies[!duplicated(Emergencies[, 1:2]), ]
#Emer <- Emer[Emer$Month == 1, ]
area.owin(Region)/(1000^2) # Area of the region

Emergencies[, X := X + rnorm(nrow(Emergencies), 0, 0.1)]

mon_full <- ppp(x = Emergencies$X, y = Emergencies$Y, window = Region)
mon_full <- affine.ppp(mon_full, diag(c(0.001, 0.001)))
unitname(mon_full) <- c("kilometre", "kilometers")
summary(mon_full)
plot(mon_full, main = "Emergencies' Locations")

quad <- quadratcount(mon_full, nx = 3, ny = 3, cex = 2)
quad2 <- quadratcount(mon_full, nx = 5, ny = 5, cex = 2)
plot(mon_full, main = "Quadrat 3 X 3")
plot(quad, add = TRUE, col = "red", cex = 2)
plot(mon_full, main = "Quadrat 5 X 5")
plot(quad2, add = TRUE, col = "red", cex = 2)

dens <- density.ppp(mon_full, sigma = 1.5, kernel = "epanechnikov")
plot(log(dens), main = "")

prob <- quantile(Emergencies$M, c(0, 0.95, 1))
cate <- cut(Emergencies$M, prob, labels = c("Lower", "Higher"),
           include.lowest = TRUE)
mon_full_wmppp <- wmppp(data.frame(X = Emergencies$X, Y = Emergencies$Y,
                                   PointType = cate), win = Region)

KdE <- KdEnvelope(X = mon_full_wmppp, r = seq(0, 10, 0.1),
                  NumberOfSimulations = 1000, ReferenceType = "Higher", Global = TRUE)
plot(KdE, main = "Envelope of KDE under Null Hypothesis")
rm(dens, KdE, mon_full_wmppp, cate, prob, quad, quad2);gc()

quad_test_res <- list()
for(i in 3:7){
  quad_test_res[[i - 2]] <- quadrat.test.ppp(mon_full, i, i)
}
quad_res <- data.frame(Chi = rep(0, times = 5), DoF = rep(0, times = 5),
                       Pval = rep(0, times = 5))
for(i in 1:nrow(quad_res)){
  quad_res[i, 1] <- quad_test_res[[i]][[1]][[1]]
  quad_res[i, 2] <- quad_test_res[[i]][[2]][[1]]
  quad_res[i, 3] <- quad_test_res[[i]][[3]]
}
quad_res

Z <- ppp(Emergencies$X, Emergencies$Y, window=Region, marks=Emergencies[,c("M", "T")])
ZM <- subset(Z, select = M)
ZC <- cut(ZM, breaks = quantile(Z$marks$M,
                                probs = c(0:8)/8),
          include.lowest = TRUE, labels = 1:8)
ZCS <- shift(rescale(ZC, 1000), origin="midpoint")

fit0 <- ppm(ZCS, ~ 1)
fit1 <- ppm(ZCS, ~ x)
fit2 <- ppm(ZCS, ~ y)
fit3 <- ppm(ZCS, ~ polynom(x, 2))
fit4 <- ppm(ZCS, ~ polynom(x, 3))
fit5 <- ppm(ZCS, ~ x + y)
fit6 <- ppm(ZCS ~ polynom(y, 2))
fit7 <- ppm(ZCS ~ polynom(y, 5))
fit8 <- ppm(ZCS ~ marks)
fit9 <- ppm(ZCS ~ marks + x + y)
fit10 <- ppm(ZCS ~ marks + polynom(x, y, 2))

anova(fit0, fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8, fit9, fit10, test="Chi")

plot(predict(fit1, type = "se"), ribside = "bottom")
plot(predict(fit2, type = "se"), ribside = "bottom")
plot(predict(fit3, type = "se"), ribside = "bottom")
plot(predict(fit4, type = "se"), ribside = "bottom")
plot(predict(fit5, type = "se"), ribside = "bottom")
plot(predict(fit10, type = "se"), ribside = "bottom")
miplot(ZCS)

rm(fit0, fit1, fit2, fit3, fit4, fit6, fit7, fit8);gc()

ZCS$window$type <- "rectangle"
miplot(ZCS)
ZCS$window$type <- "polygonal"
plot(distmap(ZCS))
plot(nndist(ZCS), ylab = "Distance")


Emer <- Emergencies[Emergencies$Month %in% 1:6, ]
Z2 <- ppp(Emer$X, Emer$Y, window=Region, marks=Emer[,c("M", "T")])
GC <- Kest(Z2, nlarge = 12000)
plot(GC, main = "K")
LC <- Lest(Z2, nlarge = 12000)
plot(LC, main = "Empirical and Theoretical L")
EG <- envelope(Z2, Kest, 39, 1)
plot(EG, main = "Pointwise Envelope")
EL <- envelope(Z2, Lest, 39, 1)
plot(EL, main = "L Envelope")
rm(GC, LC, EG, EL, ZC, ZM);gc()

ZCSEnv <- envelope(fit10, Lest, nsim = 40, global = TRUE, correction = "border")
plot(ZCSEnv)

Emer_ppp <- ppp(Emer$X, Emer$Y, window=Region)
fit <- kppm(Emer_ppp, ~ 1, "Thomas")
plot(fit)
E <- envelope(fit,Kest,nsim= 59, global = TRUE, correction = "border")
plot(E) 

lam <- predict(fit5, locations = ZCS)
Linhomo <- Linhom(ZCS, lam)
plot(Linhomo)

fit <- kppm(Emer_ppp ~ x + y, clusters = "MatClust")
plot(fit)