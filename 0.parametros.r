
#--parametros mortalidad y tasas Lognor
#  

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

#-----------------parametros ym(k) ~ Normal

 
mu = mean(cmk[n1])
sigma = sd(cmk[n1])

ia = (1+mu)^(12)-1
sigma.a = 12*sigma

sigmas <- c(0.001, 0.005, 0.01, 0.05, 0.1)*12

#-----------parametros mortalidad

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
m = 12
q = 1
iq = 0.017
C = 2.5

#------calculo esperanza de vida

fn.m = function(t){tpx.mb(t,x,pars.mb)}
ex.m = integrate(fn.m,0,110-x)$value
na = ceiling(ex.m)

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



#-----------Valor de la anualidad cierta

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










