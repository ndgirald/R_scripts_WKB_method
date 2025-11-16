
D = list()

country = c("ESP","GBR_NP")

#country = c("Spain","Great Britain","Australia",
# "Italy","France","Canada")


aprox = c("WKB","CMLB")

modelo = c("G-M","M-B", "K-M")

edad = seq(60,70)
n.edad = length(edad)
n.tau = 10
n.country = length(country)

D = list(aprox=list(WKB=list(age = double(n.edad),
     modelo = list(G-M=list(male(
       parameters = double(3),
       weights = double(600),
       lognor.param = double(2)),
       female(
       parameters = double(3),
       weights = double(600),
       lognor.param = double(2)),    
     CMLB=list(age = double(n.edad),
     modelo = list(G-M=list(male(
       parameters = double(3),
       weights = double(600),
       lognor.param = double(2)),
       female(
       parameters = double(3),
       weights = double(600),
       lognor.param = double(2))))

print(D)
str(D)

library(listviewer)
jsonedit(D)

library(rlist)
list.save(D, 'D.rdata')

library(rlist)
D = list.load('D.rdata')

G = list.load('G.rdata')

#------------ grafica lambda  x cohorte x edad

tau = matrix(0,n.tau,2)
tau[,1] = seq(1962,1991)
tau[,2] = seq(1991,2020)

#lambda = matrix(0,n.tau,n.edad)
nx=3

par(mfrow=c(1,2))
plot(tau[,2],D$genre$male$lambda[nx,,1],ylim=c(-0.6,0.6),
type='l',col='blue',lwd=2,ylab='lambda by age',xlab='cohort year')
for(j in 2:n.edad){
lines(tau[,2],D$genre$male$lambda[nx,,j],col='gray')}
lines(tau[,2],D$genre$male$lambda[nx,,j],col='red',lwd=2)

plot(tau[,2],D$genre$female$lambda[nx,,1],ylim=c(-0.6,0.2),
type='l',col='blue',lwd=2,ylab='lambda by age',xlab='cohort year')
for(j in 2:n.edad){
lines(tau[,2],D$genre$female$lambda[nx,,j],col='gray')}
lines(tau[,2],D$genre$female$lambda[nx,,j],col='red',lwd=2)


#------------ grafica 3D lambda  x cohorte x edad

library(plot3D)

pdf(file = "lambda.3D.GBR.female.male.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5) # The height of the plot in inches

par(mfrow=c(1,2))

persp3D(z=D$genre$female$lambda[nx,,], 
phi = 30, theta = 60, main=paste('females',country[nx]),
ylab='age',
xlab='tau',
zlab='lambda',border="black", lwd=0.3,
col = rgb(1,1,1,0.5))

persp3D(z=D$genre$male$lambda[nx,,], 
phi = 30, theta = 60, main=paste('males',country[nx]),
ylab='age',
xlab='tau',
zlab='lambda',border="black", lwd=0.3,
col = rgb(1,1,1,0.5))

# Step 3: Run dev.off() to create the file!
dev.off()

#-------graficas lambda x = 60, ajustada vs observada
#             por pais, mujeres

library(lattice)

tau = D$genre$female$tau

cp1 = numeric(0)
cp2 = cp1


for(nx in 1:n.country){
cp1 = c(cp1,D$genre$female$lambda[nx,,1])
cp2 = c(cp2,D$genre$female$lambda.hat[nx,])
}

freq = rep(tau,n.country)
countries = rep(country,length(tau))

np = order(freq,cp1,cp2,countries)

pdf(file = "adjuste.vs.observed.females.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5) # The height of the plot in inches


xyplot(cp1[np]+cp2[np]~freq[np]|countries, type='b', pch=c(1,19),
lty = rep(1,7),lwd=rep(1,7),
scales = list(y = list(relation = 'free')),
   main="Adjusted vs observed lambda, Females, age=60",
   ylab="lambda hat vs observed", xlab="years")



# Step 3: Run dev.off() to create the file!
dev.off()
#------------ lambda x = 60 ajustada vs observada
#             por pais, hombres


cp1 = numeric(0)
cp2 = cp1

for(nx in 1:n.country){
cp1 = c(cp1,D$genre$male$lambda[nx,,1])
cp2 = c(cp2,D$genre$male$lambda.hat[nx,])}

tau = D$genre$female$tau

freq = rep(tau,n.country)
countries = rep(country,length(tau))

np = order(freq,cp1,cp2)

pdf(file = "adjuste.vs.observed.males.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 5) # The height of the plot in inches


xyplot(cp1[np]+cp2[np]~freq[np]|countries, type=c('b'), pch=c(1,19),
lty = rep(1,7),lwd=rep(1,7),
scales = list(y = list(relation = 'free')),
   main="Adjusted vs observed lambda, Males age=60",
   ylab="lambda hat vs observed", xlab="years")


# Step 3: Run dev.off() to create the file!
dev.off()


#------------grafica tpx proyectada a 2050 para x = 60 
str(G)

country = c("Spain","Great Britain","Australia",
"Italy","France","USA","Canada")
n.edad = 99-60+1
n.tau = 2020-1991+1
n.country = length(country)
genre = c("male","female")
n.genre = 2

#-----------------utilizand lattice
cp1 = numeric(0)
cp2 = cp1
for(nx in 1:n.country){
cp1 = c(cp1,G$male$Gqcax[nx,],G$female$Gqcax[nx,])
cp2 = c(cp2,G$male$Gqcax.0[nx,],G$female$Gqcax.0[nx,])}

cpcountry = rep(country[1],1000)
for(j in 2:7){
cpcountry = c(cpcountry,rep(country[j],1000))}

cpgenre = c(rep(genre[1],500),rep(genre[2],500))
for(j in 2:7){
cpgenre = c(cpgenre,c(rep(genre[1],500),rep(genre[2],500)))}



base = data.frame(
price0 = cp2,
price = cp1,
country = cpcountry,
genre = cpgenre)

#Gqcax=readRDS("Gqcax.rda")


library(lattice)
library(latticeExtra)
library(ggplot2)


my_theme <- trellis.par.get()
my_theme$strip.background$col <- "grey80"
my_theme$plot.symbol$pch <- 16
my_theme$plot.symbol$col <- "grey60"
my_theme$plot.polygon$col <- "grey90"

dens_lattice <- densityplot(~ price+price0 | cpgenre+cpcountry, 
data = base,
as.table = TRUE,
par.settings = my_theme,
plot.points = FALSE,
between = list(x = 0.2, y = 0.2),
scales = list(x = list(rot = 45)))

l_den <- useOuterStrips(dens_lattice)

pdf(file = "present.value.projected.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height = 5) # The height of the plot in inches


print(l_den)

# Step 3: Run dev.off() to create the file!
dev.off()



par(mfrow=c(3,3))

for(nx in 1:7){
plot(density(G$male$Gqcax.0[nx,]),main=country[nx],
lty=1,lwd=2,xlim=c(0,450),ylim=c(0,0.008))
lines(density(G$male$Gqcax[nx,]),lty=5,lwd=2)
}

legend("topleft", c(country), 
lty=c(8,seq(1,6)),lwd=c(rep(2,7)))



plot(density(G$male$Gqcax.0[1,]),lty=1,lwd=2,xlim=c(-50,400),ylim=c(0,0.0081))
for(nx in 2:7){
lines(density(G$male$Gqcax.0[nx,]),lty=nx,lwd=2)}

legend("topleft", c(country), 
lty=c(8,seq(1,6)),lwd=c(rep(2,7)))



legend("topleft", c(country), 
pch=c(seq(21,24)),cex=rep(1.5,4),
pt.bg=c(rep('gray',3),'red',rep('gray',1)),
lty=c(rep("solid",4)),lwd=c(rep(1,4)),
col=c(rep('black',4)))







