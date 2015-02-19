wolff <- read.table("temperature_sweep.dat",header=T)
attach(wolff) #Attaches wolff to R search path; wolff$Value becomes Value


sus.max = which.max(Susceptibility)
cv.max = which.max(Cv)
print(paste("Susceptibility max occurs at:",Temperature[sus.max]) ,sep = " ")
print(paste("Cv max occurs at:",Temperature[cv.max]) ,sep = " ")

par(mfrow = c(2,2), oma = c(0,0,2,0), mar = c(5,5,2,2) - .3)
#opens graphics device with 2x2 matrix for plots, 2 line outer margin

plot(Temperature,Magnetization,type="p",
    xlab="Temperature",ylab="<Magnetization>")
abline(v = Temperature[sus.max],col = "red",lty = 3)
abline(v = Temperature[cv.max],col = "green",lty = 3)

plot(Temperature,Susceptibility,type="p",
    xlab="Temperature",ylab="Susceptibility")
abline(v = Temperature[sus.max],col = "red",lty = 3)
abline(v = Temperature[cv.max],col = "green",lty = 3)

plot(Temperature,Energy,type="p",
    xlab="Temperature",ylab="<Energy>")
abline(v = Temperature[sus.max],col = "red",lty = 3)
abline(v = Temperature[cv.max],col = "green",lty = 3)

plot(Temperature,Cv,type="p",
    xlab="Temperature",ylab="<Heat Capacity>")
abline(v = Temperature[sus.max],col = "red",lty = 3)
abline(v = Temperature[cv.max],col = "green",lty = 3)


mtext("Monte Carlo simulation of Ising model", line = .5, outer = T)


detach(wolff)

dev.copy2eps(file="system_plots.eps")
