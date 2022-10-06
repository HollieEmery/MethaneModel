# CH4-flux process model ####
# Written by: Hollie Emery
# Date: 2013-12-13
# Adapted to R: 2018-02-23

# based on the diffusion/advection/transport process model dscribed in 
# Walter & Heimann 2000 (Global Biogeochemical Cycles), this script adds 
# tidal forcing to the model, thus making the freshwater wetland model
# appropriate for salt marshes and other tidal wetlands

# DEFINE DOMAIN #####

#temporal values
tidal = 60*60*12.4         # high tide every 12.4 hours (s)
day = 60*60*24             # day is 24 hours (s)
year = day * 365           # year is 365 days (s)
spring = day * 90          # growing season begins April 1 (DOY 90)
fall = day * 304           # growing season ends October 31 (DOY 304) 

#model time
dt = 60*60                 # 1 hr time step (s)
tmax = year                # max time of 1 yr (s)
nt = tmax/dt+1

#model space
zmin = 0                   # position (cm)
zmax = 300                 # position (cm) - 3m tidal range
platform = 225             # marsh platform height (cm)
root = platform-50         # plant rooting depth (cm)
sed = 100                  # depth of marsh sediments
dz = 1                     # increment (cm)
z = sed:platform           # position array
rz = which(z>root)         # root zone index
nz = length(z)


# DEFINE CONSTANT PARAMS #####
  
#Diffusion
fcoarse = 0.5              # relative volume of coarse pores
Di = 0.2e-4                # methane diffusivity in water/saturated sediments (cm^2/s)
Dii = 0.2                  # methane diffusivity in air/dry sediments (cm^2/s)
Dw = Di * .66 * fcoarse
Da = Dii * .66 * fcoarse

#production 
R0 = 1                     # production rate tuning parameter
Q10p = 6                   # Methane production Q10
Tbar = 12                  # Mean annual temperature

#oxidation 
Vmax = 30/(60*60)          # Michealis-Menton max CH4 oxidation rate (µM/hr)
Km = 5                     # Michealis-Menton parameter (µM)
Q10o = 2                   # Methane oxidation Q10

#plant mediated transport
kp = 0.1/(60*60)           # rate constant (per hour)
Tveg = 10                  # quality of transport (unitless)
Pox = 0.5                  # plants oxidize half the methane they pick up

#ebullition
Cthresh = 500              # concentration of methane for bubble formation
ke = 1/(60*60)             # rate constant (per hour)


# DEFINE PROFILE PARAMS ####

#root biomass in root zone
froot=2*(z-root)/(platform-root)
froot[-rz]=0

# organic matter profile
forg=exp(-abs(z-root)/10)
forg[-rz]=0

# MODEL FORCING ####
time = seq(0,tmax,dt)

#tide is modelled with 12.4 hr and 14.7 day harmonics
t1 = cos(time*2*pi/tidal)
t2 = 3+cos(time*2*pi/(day*14.7))
tide = 150+50*t1*t2

#air temperature is modelled to vary seasonally between -5 and 29∞C
airT = 12+15*sin(time*2*pi/year+10.8)+2*sin(time*2*pi/day+4.3)

#soil temperature is interpolated between airT and constant temp at depth
Tdeep = 4                  
T = matrix(0,nt,nz) #zeros(nt,nz)
for(i in 1:nt){
  T[i,] = Tdeep + (airT[i]-Tdeep) * z/platform
}

#binary function for freezing temperatures
fT = matrix(1,nrow(T),ncol(T))#ones(size(T))
fT[which(T<=0)] = 0

# DERIVED FORCINGS ####

#plant growth
Tgr = 7
Tmat = Tgr+10
Lmin=0
L=4
Lmax=Lmin+L
fgrow=rep(0,nt)
for(i in 1:nt){
  if(T[i,50] < Tgr){
    fgrow[i] = Lmin
  }else{
    if(T[i,50] > Tmat){
      fgrow[i] = Lmax
    }else{
      fgrow[i] = Lmin + L*(1-((Tmat-T[i,50])/(Tmat-Tgr))^2)
    }
  }
}

#Evapotranspiration modelled from air temp and plant growth
ETf = .1 + 0.001 * airT^2 + 0.1 * fgrow

# water table (wt) height in cm
wt=rep(0,nt)
wt[1]=platform
for(i in 2:nt){
  if(tide[i] >= platform){
    wt[i] = platform
  }else{
    wt[i] = wt[i-1] - ETf[i]
  }
}

# BOUNDARY CONDITIONS ####

Cdeep = 0                  # soil concentration of methane at depth in µM
Catm = 0.076               # atmospheric concentration of methane in µM

# INITIAL CONDITIONS ####

C = prod=oxid=ebull=plant=H=matrix(0,nt,nz)        # initialize storage
Tot_ebull=Febull=Fplant=Fdiff=Ftot=rep(0,nt)

# Initial concentration profile built based on spin-up
C[1,] = 60000*rev(dweibull(1:nz,2,50))


# SCHEME COEFFICIENTS ####

Bw = Dw * dt / (dz^2)      # wet beta value
Ba = Da * dt / (dz^2)      # dry beta value

Aw = matrix(0,nz,nz)
# implicit differencing coefficient matrix (defaults to wet diffusion)
#for i=2:nz-1
#   Aw(i,i-1) = -Bw;
#  Aw(i,i)=(1+2*Bw);
# Aw(i,i+1) = -Bw;
#end

# Crank-Nicolson coefficient matrix (defaults to wet diffusion)
for(i in 2:(nz-1)){
  Aw[i,i-1] = -Bw
  Aw[i,i] = 2+2*Bw
  Aw[i,i+1] = -Bw
}

Aw[1,1] = 1
Aw[nz,nz] = 1

# MAIN LOOP ####

# initialize counters
time = 0                   # initialize time counter [s]
i = 1                      # initialize counter for loop
flag = 1                   # toggle to escape while loop
t = rep(0,nt)
t[1] = time                # initial time

while(flag==1){
i = i + 1              # update step counter
time = time + dt       # update time counter
t[i] = time            # update time

# diffusion
# build coefficient matrix based on water table height
A = Aw
wtable = floor(wt[i]-sed)
for(k in wtable:(nz-1)){
  A[k,k-1] = -Ba
  A[k,k] = (1 + 2*Ba)
  A[k,k+1] = -Ba
}
A[1,1] = 1
A[nz,nz] = 1

# make right hand side (known) vector
rhs = rep(0,nz)
#rhs(2:nz-1) = C(i-1,2:nz-1);    # implicit 

# (crank-nicolson style)
for(m in 2:(nz-1)){
  watertable = wt[i] - sed
  b = ifelse(m > watertable, Ba, Bw)
  rhs[m] = b * (C[i-1, m-1] + C[i-1, m+1]) + 2 * (1-b) * C[i-1, m]
}    

rhs[1] = Cdeep
rhs[nz] = Catm

# calculate unknown vector
c.unk = solve(A,rhs) #c=A\rhs

# loop through depth to calculate source/sink term
for(j in 1:nz){
  
# production and consumption
if(j < (wt[i]-sed)){
  prod[i,j] = R0 * forg[j] * fT[i,j] * Q10p^((T[i,j] - Tbar) /10)  # ch4 production in the layer
  oxid[i,j] = 0 # ch4 is not consumed below the water table
}else{
  prod[i,j] = 0 # ch4 is not produced above the water table
  oxid[i,j] = -Vmax * c.unk[j] / (Km + c.unk[j]) * Q10o^((T[i,j] - Tbar) /10) # ch4 consumption
}

# ebullition
if(c.unk[j] >= Cthresh){
  ebull[i,j] = -ke * (c.unk[j] - Cthresh)
}else{
  ebull[i,j]=0
}

# plant mediated transport
plant[i,j] = -kp * Tveg * froot[j] * fgrow[i] * c.unk[j]

# total
H[i,j] = prod[i,j] + oxid[i,j] + ebull[i,j] + plant[i,j]

}

# update concentration data
C[i,] = c.unk + H[i,]

# calculate fluxes
Tot_ebull[i] = -sum(ebull[i,])    # total removed by ebullition

if(wt[i] < platform){                 # if soil surface is exposed
# add ebullulated methane to lowest dry layer, not total flux
C[i,wtable+1] = C[i,wtable+1] + Tot_ebull[i]
Febull[i] = 0
# diffusive flux is through air-filled sediment pores
Fdiff[i] = Da * C[i,nz-1]
}else{ # if surface is inundated
# ebullition bubbles rise to the surface and join flux
Febull[i] = Tot_ebull[i]
# diffusive flux is through water-filled sediment pores
Fdiff[i] = Dw * C[i,nz-1]
}

Fplant[i] = -sum(plant[i,]) * (1-Pox)     # methane removed from the entire the root zone  
Ftot[i] = Febull[i] + Fdiff[i] + Fplant[i]  # total of flux

# have we gone past max time?
if(time >= tmax){
  flag = 0               # set flag = 0 if t >= tmax
}
}

plt = t/(60*60*24)

#figures ####
plot(plt,tide,
     type="l",
     xlim=c(0, 28),
     xlab='Time (days)',
     ylab='Water height (cm)')

#figure
matplot(1:366,T[seq(1,nt,24), seq(1,nz,25)],
        type="l",lty=1,
        xlab='Time (DOY)',
        ylab='Soil Tmperature (∞C)',
        main='Soil Temperature at 6 depths below the marsh platform')

legend("topleft",
       c('125 cm','100 cm','75 cm','50 cm','25 cm','0 cm'),
       lty=1,col=1:6)

#figure
plot(plt,T[,50],
     type="l",
     xlim = c(150,160),
     ylim = c(10,30),
     xlab='Time (DOY)',
     ylab='Soil Tmperature (∞C)')
lines(plt,T[,nz],lwd=2)
legend("bottomright",c('Surface','50 cm'),lwd=2:1)

#figure
plot(plt,wt,
     type="o",pch=19,cex=.3,
     xlim=c(0,4),
     ylim=c(214,226),
     xlab='Time (DOY)',
     ylab='Water table height (cm)',
     main='Winter')

#figure
plot(plt,wt,
     type="o",pch=19,cex=.3,
     xlim=c(200,204),
     ylim=c(214,226),
     xlab='Time (DOY)',
     ylab='Water table height (cm)',
     main='Summer')

#figure
matplot(plt,C[,seq(1,nz,10)],
        type="l",lty=1,
xlab='Time',
ylab='Concentration at 10 cm intervals')

#figure
matplot(t(C[seq(3,nt,30*24),]),z,lwd=2,
        xlim=c(0,4000),
        lty=1,type="l",
xlab='Methane Concentration',
ylab='Depth',
ylim=c(sed, platform))
legend("bottomright",c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),col=1:12,lty=1)

#figure
plot(plt,Febull,pch=".",col="magenta",
     xlim=c(0,365),ylim=c(-2,16),
     xlab='Time (DOY)',
     ylab='Flux',
     main='Relative contributions of transport pathways')
points(plt,Fdiff,col="blue",pch=".")
points(plt,Fplant,col="green",pch=".")
legend("topleft",c('Ebullition','Diffusion','Plant-Mediated'),col=c("magenta","blue","green"),lty=1)

#figure
plot(plt,Febull,col="magenta",lwd=2,type="l",
     xlim=c(210,215),ylim=c(-2,16),
     main='Relative contributions of transport pathways',
     xlab='Time (DOY)',
     ylab='Flux')
lines(plt,Fdiff,col="blue",lwd=2)
lines(plt,Fplant,col="green",lwd=2)
legend("topright",c('Ebullition','Diffusion','Plant-Mediated'),col=c("magenta","blue","green"),lty=1)

