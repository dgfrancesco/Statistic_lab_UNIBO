### Statics lab group's poject for exam at Alma Mater Studiorum University of Bologna 

# more details of meuse poject are on books: 
# Bivand, Pebesma, Gomez-Rubio "Applied Spatial Data Analysis with R"
# P.A. Burrough, R.A. McDonnell, 1998 "Principles of Geographical Information Systems"
########################################################################
# library needed   # install.packages() 
library(gstat)
library(geoR)
library(fields)
library(akima)

# Data import & geodata object
data(meuse.all)
meuse.all$zinc <- log(meuse.all$zinc)
meuse <- as.geodata(meuse.all, coords.col = c(2,3), data.col = 7)
plot(meuse)

# Creation of some variables
zinco <- meuse$data
x <- meuse$coords[,1]
y <- meuse$coords[,2]

# Option 0,H0 we have called m0 and the option, H1 that we have called m1
# Fisher's F test & Student's t test

m0<-lm(zinco~1)
m1<-lm(zinco~x+y)
summary(m1)
anova(m0,m1)

# 3D map
xgrid <- seq(min(x), max(x), length = 20)
ygrid <- seq(min(y), max(y), length = 20)
xy.grid <- expand.grid(x = xgrid, y = ygrid)

par(mfrow=c(3,3))
data.interp.lin <- interp(x, y, zinco)
data.surf.lin <- predict(m1, newdata = xy.grid)
res.interp.lin <- interp(x, y, zinco-fitted.values(m1))
?drape.plot
drape.plot(data.interp.lin$x, data.interp.lin$y, data.interp.lin$z, zlim = range(data.interp.lin$z, na.rm = TRUE), 
           zlim2 = NULL, add.legend = F, horizontal = TRUE, theta = 30, phi = 20, breaks = NA,
           ticktype="detailed", main="data", cex.main=3)

drape.plot(xgrid, ygrid, matrix(data.surf.lin, 20, 20), zlim = range(data.surf.lin, na.rm=TRUE), 
           zlim2 = NULL, add.legend = F, horizontal = TRUE, theta = 30, phi = 20, breaks = NA,
           ticktype="detailed", main="trend", cex.main=3)

drape.plot(res.interp.lin$x, res.interp.lin$y, res.interp.lin$z, zlim = range(res.interp.lin$z, na.rm = TRUE), 
           zlim2 = NULL, add.legend = F, horizontal = TRUE, theta = 30, phi = 20, breaks = NA,
           ticktype="detailed", main="residuo", cex.main=3)

##########################################################################
drape.plot(data.interp.lin$x, data.interp.lin$y, data.interp.lin$z, zlim = range(data.interp.lin$z, na.rm = TRUE), 
           zlim2 = NULL, add.legend = F, horizontal = TRUE, theta = 120, phi = 20, breaks = NA,
           ticktype="detailed", cex.main=3)

drape.plot(xgrid, ygrid, matrix(data.surf.lin, 20, 20), zlim = range(data.surf.lin, na.rm=TRUE), 
           zlim2 = NULL, add.legend = F, horizontal = TRUE, theta = 120, phi = 20, breaks = NA,
           ticktype="detailed", cex.main=3)

drape.plot(res.interp.lin$x, res.interp.lin$y, res.interp.lin$z, zlim = range(res.interp.lin$z, na.rm = TRUE), 
           zlim2 = NULL, add.legend = F, horizontal = TRUE, theta = 120, phi = 20, breaks = NA,
           ticktype="detailed", cex.main=3)
########################################################################
drape.plot(data.interp.lin$x, data.interp.lin$y, data.interp.lin$z, zlim = range(data.interp.lin$z, na.rm = TRUE), 
           zlim2 = NULL, add.legend = T, horizontal = TRUE, theta = 210, phi = 20, breaks = NA,
           ticktype="detailed", cex.main=3)

drape.plot(xgrid, ygrid, matrix(data.surf.lin, 20, 20), zlim = range(data.surf.lin, na.rm=TRUE), 
           zlim2 = NULL, add.legend = T, horizontal = TRUE, theta = 210, phi = 20, breaks = NA,
           ticktype="detailed", cex.main=3)

drape.plot(res.interp.lin$x, res.interp.lin$y, res.interp.lin$z, zlim = range(res.interp.lin$z, na.rm = TRUE), 
           zlim2 = NULL, add.legend = T, horizontal = TRUE, theta = 210, phi = 20, breaks = NA,
           ticktype="detailed", cex.main=3)
#########################################################################
drape.plot(data.interp.lin$x, data.interp.lin$y, data.interp.lin$z, zlim = range(data.interp.lin$z, na.rm = TRUE), 
           zlim2 = NULL, add.legend = F, horizontal = TRUE, theta = 30, phi = 90, breaks = NA,
           ticktype="detailed", main="data", cex.main=3)

drape.plot(xgrid, ygrid, matrix(data.surf.lin, 20, 20), zlim = range(data.surf.lin, na.rm=TRUE), 
           zlim2 = NULL, add.legend = F, horizontal = TRUE, theta = 30, phi = 90, breaks = NA,
           ticktype="detailed", main="trend", cex.main=3)

drape.plot(res.interp.lin$x, res.interp.lin$y, res.interp.lin$z, zlim = range(res.interp.lin$z, na.rm = TRUE), 
           zlim2 = NULL, add.legend = F, horizontal = TRUE, theta = 30, phi = 90, breaks = NA,
           ticktype="detailed", main="residuo", cex.main=3)
#########################################################################
# Variogram of residues
# with m0 valid(trend = "cte")
# with m1 valid(trend = "1st") 

par(mfrow=c(2,2))
plot(variog(meuse,trend = "cte"), main="m0")
plot(variog(meuse,trend = "cte",option="cloud"), main="m0")
plot(variog(meuse,trend = "1st"), main="m1")
plot(variog(meuse,trend = "1st",option="cloud"), main="m1")

# The best variogram would appear to be the trend m1
# Vari0 is the variogam obteined by applying a constant model to the data
# Vari is the variogram obteined from applying a linear model to data

par(mfrow=c(2,2))
n.bins<-20
vari0 <- variog(meuse, trend="cte", option="bin", uvec = n.bins)
vari <- variog(meuse, trend="1st", option="bin", uvec = n.bins)
plot(vari0, main="Trend costante")
plot(vari, main="Trend lineare")

# The maximum range appears to break off at 1700

vari0 <- variog(meuse, trend="cte", option="bin", 
                max.dist = 1700, uvec = n.bins)
vari <- variog(meuse, trend="1st", option="bin", 
               max.dist = 1700,uvec = n.bins)
plot(vari0)
plot(vari)
eyefit(vari0)
eyefit(vari)
# eyefit results:
#    cov.model sigmasq phi tausq kappa kappa2  practicalRange
# 1 exponential     0.6 350  0.06  <NA>   <NA> 1048.5062957203
#  cov.model sigmasq  phi tausq kappa kappa2 practicalRange
# 1 spherical     0.4 1100  0.12  <NA>   <NA>           1100
#    cov.model sigmasq phi tausq kappa kappa2 practicalRange
# 1 exponential    0.51 500  0.07    NA     NA       1497.866

expon.trend.cte <- likfit(meuse, trend = "cte", ini=c(0.57,550),
                          nugget = 0.05, cov.model = "matern", kappa = 0.5)
spheri<-likfit(meuse, trend=~x+y, ini=c(0.37, 1300), 
               nugget=0.16, cov.model = "spherical")
expon <- likfit(meuse, trend = ~x+y, ini=c(0.51,660), nugget = 0.06, 
                cov.model = "matern", kappa = 0.5)

par(mfrow=c(2,1))
plot(vari0, ylim=c(0,1), main = "vari0")
lines(expon.trend.cte, lty=3, col=4)
plot(vari, ylim=c(0,1), main = "vari")
lines(spheri, lty=3, col=2)
lines(expon, lty=3, col=3)

# Creation of regular grid
krsfe<-pred_grid(x,y,by=30)

par(mfrow=c(1,1))
points(meuse,col="red")
points(krsfe,pch="+")

# Kriging using 3 models (interpolation)
par(mfrow=c(1,2))

# Kriging linear sferic
par(mfrow=c(1,2))
KC.sph <- krige.control(type="SK",obj.mod=spheri, 
                        trend.d="1st",trend.l="1st")
kr.sph <- krige.conv(meuse, krige=KC.sph, loc=krsfe, 
                     output=output.control(signal=TRUE))

image(kr.sph,value=kr.sph$predict,col=topo.colors(50), main="Prediction")
contour(kr.sph, add=TRUE, nlevels=10)
points(meuse,add=T)
image(kr.sph,value=sqrt(kr.sph$krige.var),col=topo.colors(50), main="St.err.")
points(meuse,add=T)

# Kriging linear exponential (matern)
KC.exp<- krige.control(type="SK",obj.mod=expon, 
                       trend.d="1st",trend.l="1st")
kr.exp <- krige.conv(meuse,krige=KC.exp, loc=krsfe, 
                     output=output.control(signal=TRUE))

image(kr.exp,value=kr.exp$predict,col=topo.colors(50))
contour(kr.exp, add=TRUE, nlevels=10)
points(meuse,add=T)
image(kr.exp,value=sqrt(kr.exp$krige.var),col=topo.colors(50))
points(meuse,add=T)

# Kriging constant exponential
KC.exp.cte<- krige.control(type="SK", obj.mod=expon.trend.cte, 
                           trend.d="cte",trend.l="cte")
kr.exp.cte <- krige.conv(meuse,krige=KC.exp.cte, loc=krsfe, 
                         output=output.control(signal=TRUE))

image(kr.exp.cte,value=kr.exp.cte$predict,col=topo.colors(50))
contour(kr.exp.cte, add=TRUE, nlevels=10)
points(meuse,add=T)
image(kr.exp.cte,value=sqrt(kr.exp.cte$krige.var),col=topo.colors(50))
points(meuse,add=T)


# Comparison of predictive goodness of the 3 models with crossvalidation
cross.spheri <- xvalid(meuse,model = spheri)
cross.expon <- xvalid(meuse,model = expon)

# Both give the matrix singularity error to be inverted in the maximum likelihood algorithm
# They return:
#xvalid: number of data locations       = 164
#xvalid: number of validation locations = 164
#xvalid: performing cross-validation at location ... 1, 
#Error in solve.default(ttivtt, crossprod(ivtt, as.vector(data))) : 
 # system is computationally singular: reciprocal condition number = 
 # Modello sferico: 7.98184e-17
 # Modello esponenziale: 1.03251e-16

cross.expon.trend.cte <- xvalid(meuse,model = expon.trend.cte)
plot(cross.expon.trend.cte)
#Returns:
#xvalid: number of data locations       = 164
#xvalid: number of validation locations = 164
#xvalid: performing cross-validation at location ... [1-164]
#xvalid: end of cross-validation

mean(cross.expon.trend.cte$error^2 / cross.expon.trend.cte$krige.var)
#[1] 0.9849874

summary(cross.expon.trend.cte$error)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-1.249000 -0.217600 -0.018820 -0.002658  0.191400  2.272000 

# Compare the 3 models based on the goodness of fit to the data by comparing AIC

spheri$AIC
#[1] 256.2527

expon$AIC
#[1] 258.3124

expon.trend.cte$AIC
#[1] 262.7259

##########################################################################
# according to the AIC criterion (which takes into account the complexity of the model)
# the best model is spheri, even if the values of spheri and expon are close
# let's say we have enough evidence that a model including a linear trend in the x and y covariates is appropriate, 
# even looking at the significance of the regression coefficients in summary (m1)
