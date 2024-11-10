# Then, we explored the knots parameter to enhance the fitting effect
# when knots = 10
lPriors_c1 <- list(scale.mu=log(mu1_c1) - (sigma1_c1^2)/2, scale.sd=sigma1_c1, bl.smooth=20^11, bl.knots=10,
                   beta.mu=mu2_c1, beta.sd=sigma2_c1 , noise.sd=0.1, noise.nu=4)

tm_c1 <- system.time(result_c1 <- fitSpectraSMC(shift_c1, raw_spectra_c1, peakLocations_c1, lPriors_c1))
bl.est_c1 <- result_c1$basis %*% result_c1$alpha[,1,] 
b_c1 <- apply( bl.est_c1, 1, mean) 
b_c1 <- as.matrix(b_c1) 
expFn_c1 <- apply(result_c1$expFn,2,mean)

raw_spectra_mean_c1 <- apply(raw_spectra_c1, 2, mean)
plot(shift_c1, raw_spectra_mean_c1, type='l',xlab=expression(paste("Raman shift (cm"^{-1}, ")")), ylab="Intensity (a.u.)", ylim = c(0,1.2), main = 'the Spectrum for Cyclohexane-1 when knots=10' )
lines(shift_c1, b_c1, col="red") #baseline
lines(shift_c1, b_c1 + expFn_c1,col="blue") #spectral signature
legend("topright", inset = .05, c("baseline","baseline + spectral signature"),lty = c(1,1), col = c("red", "blue"))
abline(v = peakLocations_c1, col = "green", lty = 3)

intensity_shift_c1 <- data.frame( raman_shift = round(shift_c1,0), intensity =  b_c1 + expFn_c1 )
for ( i in 1:length(peakLocations_c1)){
  j = peakLocations_c1[i]
  text(j, intensity_shift_c1$intensity[which(intensity_shift_c1$raman_shift == j)], j)
}
intensity_shift_c1 <- data.frame( raman_shift = ceiling(shift_c1), intensity =  b_c1 + expFn_c1 )

#when knots = 50

lPriors_c1_50 <- list(scale.mu=log(mu1_c1) - (sigma1_c1^2)/2, scale.sd=sigma1_c1, bl.smooth=10^11, bl.knots=50,
                      beta.mu=mu2_c1, beta.sd=sigma2_c1 , noise.sd=0.1, noise.nu=4)

tm_c1_50 <- system.time(result_c1_50 <- fitSpectraSMC(shift_c1, raw_spectra_c1, peakLocations_c1, lPriors_c1_50))
bl.est_c1_50 <- result_c1_50$basis %*% result_c1_50$alpha[,1,] 
b_c1_50 <- apply( bl.est_c1_50, 1, mean) 
b_c1_50 <- as.matrix(b_c1_50) 
expFn_c1_50 <- apply(result_c1_50$expFn,2,mean)

raw_spectra_mean_c1 <- apply(raw_spectra_c1, 2, mean)
plot(shift_c1, raw_spectra_mean_c1, type='l',xlab=expression(paste("Raman shift (cm"^{-1}, ")")), ylab="Intensity (a.u.)", ylim = c(0,1.2), main = 'the Spectrum for Cyclohexane-1 when knots=10' )
lines(shift_c1, b_c1_50, col="red") #baseline
lines(shift_c1, b_c1_50 + expFn_c1_50,col="blue") #spectral signature
legend("topright", inset = .05, c("baseline","baseline + spectral signature"),lty = c(1,1), col = c("red", "blue"))
abline(v = peakLocations_c1, col = "green", lty = 3)

intensity_shift_c1 <- data.frame( raman_shift = round(shift_c1,0), intensity =  b_c1_50 + expFn_c1_50 )
for ( i in 1:length(peakLocations_c1)){
  j = peakLocations_c1[i]
  text(j, intensity_shift_c1$intensity[which(intensity_shift_c1$raman_shift == j)], j)
}

intensity_shift_c1 <- data.frame( raman_shift = ceiling(shift_c1), intensity =  b_c1_50 + expFn_c1_50 )


#when knots = 100

lPriors_c1_100 <- list(scale.mu=log(mu1_c1) - (sigma1_c1^2)/2, scale.sd=sigma1_c1, bl.smooth=10^11, bl.knots=100,
                       beta.mu=mu2_c1, beta.sd=sigma2_c1 , noise.sd=0.1, noise.nu=4)

tm_c1_100 <- system.time(result_c1_100 <- fitSpectraSMC(shift_c1, raw_spectra_c1, peakLocations_c1, lPriors_c1_100))
bl.est_c1_100 <- result_c1_100$basis %*% result_c1_100$alpha[,1,] 
b_c1_100 <- apply( bl.est_c1_100, 1, mean) 
b_c1_100 <- as.matrix(b_c1_100) 
expFn_c1_100 <- apply(result_c1_100$expFn,2,mean)