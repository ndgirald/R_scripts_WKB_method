
#--Ejemplo calculo valor anualidad vida tasas Lognor
#  aprox WKB

D = read.csv("FIRMX.csv",header = TRUE, stringsAsFactors=FALSE)

u = D$Close

fechas = as.Date(D$Date,format="%Y-%m-%d")


np = length(u)

ejex.mes = seq(fechas[1],fechas[np], "months")
ejex.año = seq(fechas[1],fechas[np],"years")

#----------rendimientos diarios m = 360
imk = diff(log(u),1,1)

plot(fechas[-1],imk, xaxt="n", 
panel.first = grid(),lwd=1,col='gray',
type='l',ylab='tasa rendimiento diaria')
axis.Date(1, at=ejex.mes, format="%m/%y")
axis.Date(1, at=ejex.año, labels = FALSE, tcl = -0.2)
abline(h=0)

#-----------agregar desde dia a mes
library(highfrequency)
library(xts)
imk = xts(x=imk, order.by = fechas[-1])
#------------- definir una tasa promedio positiva !
library(GMCM)

plot(fechas[-1],GMCM:::cummean(imk), xaxt="n", panel.first = grid(),
type='l',ylab='valor unidad fondo FTHRX')
axis.Date(1, at=ejex.mes, format="%m/%y")
axis.Date(1, at=ejex.año, labels = FALSE, tcl = -0.2)
abline(h=0)

ts.month <- apply.monthly(imk,FUN=sum)
dmk = log(1+ts.month)
cmk = GMCM:::cummean(dmk)

ts.plot(cmk)
abline(h=0)

n1=which(cmk>0)

mimk = mean(cmk[n1])
ia = (1+mimk)^(12)-1

#-----------------parametros ym(k) ~ Normal
mu = mimk
sigma = sd(cmk[n1])



#-- ley Makeham Beard

muxt.mb = function(t,x,pars){
a = pars[1]; b = pars[2]; k=pars[3]; r =pars[4];
(k+ a*exp(b*(x+t))  )/( 1 +a*r*exp(b*(x+t)) )
}

tpx.mb = function(t,x,pars){
a = pars[1]; b = pars[2]; k=pars[3]; r =pars[4];
f=(1+a*r*exp(b*x))/(1+a*r*exp(b*(x+t)))
p=ifelse(t < 110-x,exp(-k*t)*f^((1-k*r)/(b*r)),0)
return(p)
}
#-------------parametros MB

pars.mb = c(0.00004720,0.09048063,0.00016508,0.02963878)

#------------------para simular T(x), Kmx
require(GoFKernel)
f <- function(t) 1-tpx.mb(t,x,pars.mb)

#------------parametros de la anualidad

x = 60; w = 110;

#------calculo esperanza de vida

fn.m = function(t){tpx.mb(t,x,pars.mb)}
ex.m = integrate(fn.m,0,110-x)$value
na = ceiling(ex.m)

#-----------Valor de la anualidad cierta
m = 12
q = 1
iq = 0.017
C = 2.5

Gavqmn = function(i,m,q,n,iq){
try(if(iq > i) stop("tasa inflacion invalida"))
try(if(m%%q != 0) stop("m no es divisible por q"))
t = seq(1,n*m,1)
res = (1/m)*sum((1+i)^(-t/m)*(1+iq)^(floor((t-1)*q/m)/q))
return(res)}

Cp = double(4)
names(Cp) = c("Geo.cierta","Geo.vida","Geo.vida.Lognor","Geo.vida.Lognor.sim")

(Cp[1] = C*m*Gavqmn(ia,m,q,na,iq))

#-----Valores de la anualidad de vida
Gavqxm = function(x,i,iq,m,q,pars){
try(if(iq > i) stop("tasa inflacion invalida"))
try(if(m%%q != 0) stop("m no es divisible por q"))
v = 1/(1+i)
k = seq(1,m*(110-x))
kmpx = sapply(k,function(k)tpx.mb(k/m,x,pars))
vkm = v^(k/m)
vqm = (1+iq)^(floor(k*q/m)/q)
a = (1/m)*v^(1/m)*sum(vkm*vqm*kmpx)
return(a)}

(Cp[2] = C*m*Gavqxm(x,ia,iq,m,q,pars.mb))

#-------------define los pagos rk

nm = (110-x)*m
t = seq(1,nm)
rk = C*(1+iq)^floor((t-1)/m)

#-------determina la distribucion de Kmx+1/m = techo
tm = t/m
dKmx = tpx.mb(tm,x,pars.mb)*(1-tpx.mb(1/m,x+tm,pars.mb))
(A=sum(dKmx))
dKmx = dKmx/A
plot(tm,dKmx,type='l')

#
#-------------------Metodo WKB la fda de valor presente
#                   anualidad de vida tasas LogNor
#                   
library(lognorm)


#--------------media y matriz cov

mu.w = -mu*seq(1,nm)+log(rk)

R = matrix(0,nm,nm)
for(i in 1:nm){
for(j in 1:nm){
R[i,j] = ifelse(i <= j, i,j)}}
R = matrix(sigma^2,nrow(R),ncol(R))*R
Rho = cov2cor(R)

#-----------------determina los parametros LogNormal
#                 del metodo WKB, 
muS = double(nm)
sdS = double(nm)


for(j in 1:nm){
mu.j = mu.w[1:j]
R.j = R[1:j,1:j]
Rho.j = Rho[1:j,1:j]
coefSum2 = estimateSumLognormal(mu.j,
sqrt(diag(R.j)),corr=Rho.j )
muS[j] = coefSum2[1]
sdS[j] = coefSum2[2]
}

head(muS)
head(sdS)

muS[1] = 1.23
sdS[1] = 0.011

#-------------simula Z con mezcla LogNor y Kmx

N = 5000
p <- dKmx

# install.packages("mixtools")
library(mixtools)

sim_data <- rnormmix(n=N,  lambda = p, mu = muS, sigma = sdS)

# The simulated data 'sim_data' is on the log scale.
# To get the data on the original scale (non-log normal), exponentiate:

Zs <- exp(sim_data)

Cp[3] = mean(Zs)

#-------------determina la fdp mezcla LogNor y Kmx
#-------------del valor presente de  S

fdpS = function(t,muS,sdS,dKmx){
nm = length(muS)
fdpL = double(nm)
for(j in 1:nm){
fdpL[j] = dKmx[j]*dlnorm(t,meanlog = muS[j],sdlog=sdS[j])
}
p=sum(fdpL)
return(p)}

#----------------no usar
fdpS = function(t,muS,sdS,dKmx){
nm = length(muS)
fdpL = double(nm)
for(j in 1:nm){
fdpL[j] = d.sum.of.mixtures.LNLN(t,1,p.vector=dKmx,
mu.vector=muS,sigma.vector=sdS,logdens=TRUE)
}
}
#----------------------

library(rmutil)
fn = function(t) fdpS(t,muS,sdS,dKmx)
(A = int(fn,0,950))
hn = function(t) t*fdpS(t,muS,sdS,dKmx)/A
(Cp[4] = int(hn,0,950))

#---------------------histograma y densidad
tz = seq(0,1250)
nz = length(tz)
fz = double(nz)
for(j in 1:nz){
fz[j] = fdpS(tz[j],muS,sdS,dKmx)/A}

#for(j in 1:nz){
#fz[j] = d.sum.of.mixtures.LNLN(300,1,p.vector=dKmx,
#mu.vector=muS,sigma.vector=sdS,logdens=TRUE)}


#-----modificar parte de la cola izquierda
#      para corregir discontinuidades
tz0 = seq(1,150)
lfz = fz[1:150]
mlofz = loess(lfz~tz0, degree = 2)
lfz.s = predict(mlofz, data.frame(tz0 = tz0), se = TRUE)
fz[1:150] = lfz.s$fit


hist(Zs,100,prob=TRUE)
points(Cp,rep(0,4),pch=rep(20,4),cex=rep(1.5,4),
col=c('blue','red','orange','magenta'))
lines(tz,fz,col='red')

#----------------fda y funcion cuantil

fdaS = function(t,muS,sdS,dKmx){
mn = length(muS)
fdaL = double(mn)
for(j in 1:mn){
fdaL[j] = dKmx[j]*plnorm(t,meanlog = muS[j],sdlog=sdS[j])
}
p=sum(fdaL)
return(p)}

Fs.e <- ecdf(Zs)
plot(Fs.e, main = "Empirical Cumulative Distribution Function",
     xlab = "Value", ylab = "Cumulative Probability")

ktz = knots(Fs.e)
nz = length(ktz)
fdaz = double(nz)
for(j in 1:nz){
fdaz[j] = fdaS(ktz[j],muS,sdS,dKmx)}

plot(ktz,Fs.e(ktz),type='s')
lines(ktz,fdaz,col='red')

fdaS(0,muS,sdS,dKmx)
fdaS(1000,muS,sdS,dKmx)

#---------------probabilidad de insolvencia
(PS = 1- fdaS(Cp[4],muS,sdS,dKmx))


#----------------funcion cuantil de S

fS = function(t)fdaS(t,muS,sdS,dKmx)
fS.inv <- inverse(fS,lower=0,upper=1000)
fS.inv(0.99)

#--------------------------grafica cuantil S
Qs = double(length(p))
for(j in 1:length(p)){
Qs[j] = fS.inv(p[j])}

plot(p,Qs,type='l')

#--------------------momentos de S con Q(p)
library(rmutil)

(ES = int(fS.inv,0,1))

fS2 = function(s)(fS.inv(s)-ES)^2

(VarS = int(fS2,0,1))

#-----------momentos de S con 1-F

fS = function(s)1-fdaS(s,muS,sdS,dKmx)
(ES = int(fS,0,1000))
(mean(Zs))

fS2 = function(s)2*s*fS(s)

(ES2 = int(fS2,0,1000))
(VarS = ES2-(ES)^2)
(var(Zs))


#-------fda con cota inferior comonotonicidad

# redefinir mu para no incluír log(rk)


mu.j = mu*seq(1,nm)
sigma2.j = diag(R)

a = tpx.mb(tm,x,pars.mb)*rk*exp(-mu.j+(1/2)*sigma2.j)

Var.Lambda <- t(a) %*% R %*% a

rho.j = double(length(a))

for(j in 1:length(a)){
S.j = tpx.mb(tm,x,pars.mb)*rk*sqrt(sigma2.j)*Rho[j,]*exp(-mu.j+(1/2)*sigma2.j)
rho.j[j] = sum(S.j)/sqrt(Var.Lambda)
}

#------------simular Sl, ec. (5.19)
f = function(t){
ifelse(t < 110-x, 1-tpx.mb(t,x,pars.mb),1)}


Zs.l = double(N)
for(j in 1:N){
Tx = random.function(n=1, f, lower = 0, upper = 110-x, 
kind = "cumulative")
k = ceiling(m*Tx)
mu.k = -mu.j[1:k]+(1/2)*(1-rho.j[1:k]^2)*sigma2.j[1:k]
sigma.k = sqrt(sigma2.j[1:k])*rho.j[1:k]
Zs.l[j] = sum(rk[1:k]*exp(mu.k+sigma.k*rep(rnorm(1,0,1),k)))
}

Cp = c(Cp,mean(Zs.l))
names(Cp[5])="Geo.vida.simul.Sl"

hist(Zs.l,100,prob=TRUE)
points(Cp,rep(0,5),pch=rep(20,5),cex=rep(1.5,5),
col=c('blue','red','orange','magenta','green'))
lines(tz,fz,col='red')


Fsl.e <- ecdf(Zs.l)

ktz = knots(Fs.e)
nz = length(ktz)
fdaz = double(nz)
for(j in 1:nz){
fdaz[j] = fdaS(ktz[j],muS,sdS,dKmx)}

plot(ktz,Fs.e(ktz),type='s')
lines(ktz,Fsl.e(ktz),col='blue')
lines(ktz,fdaz,col='red')

n1 = 2501
n2 = 4500
plot(ktz[n1:n2],Fs.e(ktz[n1:n2]),type='s')
lines(ktz[n1:n2],Fsl.e(ktz[n1:n2]),col='blue')
lines(ktz[n1:n2],fdaz[n1:n2],col='red')


#--------------cuantil condicional de la cota Sl
cond.quant.Sl.k = function(p,k,mu.j,rho.j,sigma2.j){
mu.k = -mu.j[1:k]+(1/2)*(1-rho.j[1:k]^2)*sigma2.j[1:k]
sigma.k = sqrt(sigma2.j[1:k])*rho.j[1:k]
s = sum(rk[1:k]*qlnorm(p,meanlog=mu.k,sdlog=sigma.k))
return(s)
}

cond.quant.Sl.k(p=0.5,k=100,mu.j,rho.j,sigma2.j)

# para cada y resolver una ecuacion no lineal

y = 227.4728
k = 100

fn = function(p){
cond.quant.Sl.k(p,k,mu.j,rho.j,sigma2.j)}

p= seq(0.01,0.99,0.01)

wn = double(length(p))
for(j in 1:length(p)){
wn[j] = fn(p[j])}

plot(p,wn,type='l')
abline(h=0)

#------
Zk.y = function(y){
fnToFindRoot = function(p) fn(p)-y
uniroot(fnToFindRoot, 
c(0, 1), tol = 0.001)$root
}

Zk.y(226)
Zk.y(230)
Zk.y(430)

#-------------la fda de la cota inferior Sl

fdaSl = function(y,mu.j,rho.j,sigma2.j,dKmx){
fdaL = double(nm)
for(k in 1:nm){

fn = function(p){
cond.quant.Sl.k(p,k,mu.j,rho.j,sigma2.j)}

Zk.y = function(y){
fnToFindRoot = function(p) fn(p)-y
uniroot(fnToFindRoot, 
c(0, 1), tol = 0.001)$root
}

fdaL[k] = dKmx[k]*Zk.y(y)
}
p=sum(fdaL)
return(p)}

fdaSl(Cp[3],mu.j,rho.j,sigma2.j,dKmx)

#----------------funcion cuantil de Sl

fSl = function(t)fdaSl(t,mu.j,rho.j,sigma2.j,dKmx)


fSl.inv <- inverse(fSl,lower=0,upper=1000)
fSl.inv(0.99)

(PS = 1- fdaS(Cp[4],muS,sdS,dKmx))



#-----------momentos de Sl con Q(p)

library(rmutil)

(ES = int(fS.inv,0,1))

fS2 = function(s)(fS.inv(s)-s)^2

(VarS = int(fS2,0,1))


#--------------------------grafica cuantil Sl
Qsl = double(length(p))
for(j in 1:length(p)){
Qsl[j] = fSl.inv(p[j])}

plot(p,Qs,type='l')
lines(p,Qsl,col='red')

#-----------momentos de Sl con 1-F
fSl = function(s)1-fdaSl(s,mu.j,rho.j,sigma2.j,dKmx)
(ESl = int(fSl,0,900))

fSl2 = function(s)2*s*fSl(s)

(ESl = int(fSl,0,900))
(ESl2 = int(fSl2,0,900))
(VarSl = ESl2-(ESl)^2)
(VarS)

#---------------------graficas de comparacion
ktz = knots(Fs.e)
nz = length(ktz)
fdaz = double(nz)
fdazl = double(nz)
for(j in 1:nz){
fdaz[j] = fdaS(ktz[j],muS,sdS,dKmx)
fdazl[j] = fdaSl(ktz[j],mu.j,rho.j,sigma2.j,dKmx)
}

plot(ktz,Fs.e(ktz),type='s')
lines(ktz,Fsl.e(ktz),col='blue')
lines(ktz,fdaz,col='red')
lines(ktz,fdazl,col='orange')

n1 = 2501
n2 = 4500
plot(ktz[n1:n2],Fs.e(ktz[n1:n2]),type='s')
lines(ktz[n1:n2],Fsl.e(ktz[n1:n2]),col='blue')
lines(ktz[n1:n2],fdaz[n1:n2],col='red')
lines(ktz[n1:n2],fdazl[n1:n2],col='orange')


#-----------invertir la funcion cuantil 
library(Qest)

invQ

#------------------------------------
library(pracma)
shat.x = function(x){
fn = function(t) Kx.1(t)-x
bisect(fn,-0.049, 0.049)$root}
shat = shat.x(300)

install.packages("nleqslv")
library(nleqslv)
y.p = nleqslv(F0, fn, method="Newton", 
global="cline", control=list(trace=1,stepmax=5))
(y.p$x)
#----------------------------
  install.packages("rootSolve")
    library(rootSolve)
  solution <- multiRoot(f = fn, start = F0)
    print(solution)



#--------------probabilidades de ruina 

#--------------usa las fda de S por WKB y S cota inferior
#              fdaS(ktz[j],muS,sdS,dKmx)
#              fdaSl(ktz[j],mu.j,rho.j,sigma2.j,dKmx)

#-------------el caso de WKB
#-------------hay que definir las fda-wkb de cada S(k)
#-------------distribucion fdp de T
#-------------requiere un valor F(0) que es uno en Cp

fdpT = function(muS,sdS,F0){
n = length(muS)
p = double(n)
p[1] = 1-plnorm(F0,meanlog = muS[1],sdlog=sdS[1])
for(k in 2:n){
p[k] = plnorm(F0,meanlog = muS[k-1],sdlog=sdS[k-1])-
plnorm(F0,meanlog = muS[k],sdlog=sdS[k])}
return(p)}

T = fdpT(muS,sdS,F0=Cp[4])
nT = which( T > 0.0001)

plot(tm,dKmx,type='l')
lines(tm[nT],T[nT]*max(dKmx)/max(T[nT]),type='h',col='gray')

(p.ruina = sum(tpx.mb(tm,x,pars.mb)*T))

T.ruina = tpx.mb(tm,x,pars.mb)*T/p.ruina


#-------------el caso de Comonotonicidad


fdpTl = function(mu.j,rho.j,sigma2.j,F0){
p = double(nm)
for(k in 2:nm){

fn.k = function(p){
cond.quant.Sl.k(p,k,mu.j,rho.j,sigma2.j)}

fn.k1 = function(p){
cond.quant.Sl.k(p,k-1,mu.j,rho.j,sigma2.j)}

Zk.y = function(y){
fnToFindRoot = function(p) fn.k(p)-y
uniroot(fnToFindRoot, 
c(0, 1), tol = 0.001)$root
}

Zk1.y = function(y){
fnToFindRoot = function(p) fn.k1(p)-y
uniroot(fnToFindRoot, 
c(0, 1), tol = 0.001)$root
}
p[k] = Zk1.y(F0)-Zk.y(F0)
}
return(p)}


Tl = fdpTl(mu.j,rho.j,sigma2.j,F0=Cp[4])
nTl = which( Tl > 0.0001)



plot(tm,dKmx,type='l')
lines(tm[nT],T[nT]*max(dKmx)/max(T[nT]),type='h',col='gray')

lines(tm[nTl],Tl[nTl]*max(dKmx)/max(Tl[nTl]),type='h',col='lightblue')


(p.ruina.l = sum(tpx.mb(tm,x,pars.mb)*Tl))

