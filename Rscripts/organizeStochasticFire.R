#title: "organizeStochasticFire"
#author: "ErinJHanan"

# FireSizes columns are:
# 1. the number of patches burned
# 2. year
# 3. month
# 4. wind_direction
# 5. wind_speed
# 6. number of ignitions during that month.

# read in results from batched runs
rs.fs1 = read.table("real.suppressed/stochastic/FireSizes.txt", head=F)
ru.fs1 = read.table("real.unsuppressed/stochastic/FireSizes.txt", head=F)
ws.fs1 = read.table("withoutHuman.suppressed/stochastic/FireSizes.txt", head=F)
wu.fs1 = read.table("withoutHuman.unsuppressed/stochastic/FireSizes.txt", head=F)

rs.fs2 = read.table("real.suppressed/stochastic2/FireSizes.txt", head=F)
ru.fs2 = read.table("real.unsuppressed/stochastic2/FireSizes.txt", head=F)
ws.fs2 = read.table("withoutHuman.suppressed/stochastic2/FireSizes.txt", head=F)
wu.fs2 = read.table("withoutHuman.unsuppressed/stochastic2/FireSizes.txt", head=F)

rs.fs3 = read.table("real.suppressed/stochastic3/FireSizes.txt", head=F)
ru.fs3 = read.table("real.unsuppressed/stochastic3/FireSizes.txt", head=F)
ws.fs3 = read.table("withoutHuman.suppressed/stochastic3/FireSizes.txt", head=F)
wu.fs3 = read.table("withoutHuman.unsuppressed/stochastic3/FireSizes.txt", head=F)

rs.fs4 = read.table("real.suppressed/stochastic4/FireSizes.txt", head=F)
ru.fs4 = read.table("real.unsuppressed/stochastic4/FireSizes.txt", head=F)
ws.fs4 = read.table("withoutHuman.suppressed/stochastic4/FireSizes.txt", head=F)
wu.fs4 = read.table("withoutHuman.unsuppressed/stochastic4/FireSizes.txt", head=F)

# combine results from batched runs
rs.fs = rbind(rs.fs1, rs.fs2)
rs.fs = rbind(rs.fs, rs.fs3)
rs.fs = rbind(rs.fs, rs.fs4)

ru.fs = rbind(ru.fs1, ru.fs2)
ru.fs = rbind(ru.fs, ru.fs3)
ru.fs = rbind(ru.fs, ru.fs4)

ws.fs = rbind(ws.fs1, ws.fs2)
ws.fs = rbind(ws.fs, ws.fs3)
ws.fs = rbind(ws.fs, ws.fs4)

wu.fs = rbind(wu.fs1, wu.fs2)
wu.fs = rbind(wu.fs, wu.fs3)
wu.fs = rbind(wu.fs, wu.fs4)

# subset records where a fire occurred
rs.fs.fire = subset(rs.fs, rs.fs$V1>0)
ru.fs.fire = subset(ru.fs, ru.fs$V1>0)
ws.fs.fire = subset(ws.fs, ws.fs$V1>0)
wu.fs.fire = subset(wu.fs, wu.fs$V1>0)

# combine outputs for each method
a=rs.fs.fire$V1
b=ru.fs.fire$V1
c=ws.fs.fire$V1
d=wu.fs.fire$V1
w=rep(1,length(a))
x=rep(2,length(b))
y=rep(3,length(c))
z=rep(4,length(d))

FireSize=append(a, b, after = length(a))
FireSize=append(FireSize, c, after = length(FireSize))
FireSize=append(FireSize, d, after = length(FireSize))

method=append(w, x, after = length(w))
method=append(method, y, after = length(method))
method=append(method, z, after = length(method))

e=rep(1,(length(a)+length(b)))
f=rep(0,(length(c)+length(d)))
ClimateChange=append(e, f, after = length(e))

g=rep(1,length(a))
h=rep(0,length(b))
i=rep(1,length(c))
j=rep(0,length(d))
Suppression=append(g, h, after = length(g))
Suppression=append(Suppression, i, after = length(Suppression))
Suppression=append(Suppression, j, after = length(Suppression))

FireSizes.all = as.data.frame(cbind(method,FireSize,ClimateChange,Suppression))
FireSizes.rs = subset(FireSizes.all, FireSizes.all$method==1)
FireSizes.ru = subset(FireSizes.all, FireSizes.all$method==2)
FireSizes.ws = subset(FireSizes.all, FireSizes.all$method==3)
FireSizes.wu = subset(FireSizes.all, FireSizes.all$method==4)

# Determine number of ignitions
numSims = 100
simYears = 37
totalYears = numSims * simYears
count(method)

# sub in correct number for each method below
rs.firesPerSim = 3794/numSims
ru.firesPerSim = 3862/numSims
ws.firesPerSim = 3237/numSims
wu.firesPerSim = 3233/numSims

rs.firesPerSim
ru.firesPerSim
ws.firesPerSim
wu.firesPerSim

# determine largest fire in each simulation (assumes >500 patches)
rs.largest = subset(rs.fs1, rs.fs1$V1>500)
rs.largest
rs.sim.num=8364/(simYears*12)
rs.sim.num

ru.largest = subset(ru.fs1, ru.fs1$V1>500)
ru.largest
ru.sim.num=6144/(simYears*12)
ru.sim.num

ws.largest = subset(ws.fs1, ws.fs1$V1>500)
ws.largest
ws.sim.num=7032/(simYears*12)
ws.sim.num

wu.largest = subset(wu.fs1, wu.fs1$V1>500)
wu.largest
wu.sim.num=11028/(simYears*12)
wu.sim.num





