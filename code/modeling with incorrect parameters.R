# Next, we explored the flexibility of the model by offering the inaccurate peak locations (intervals=5,10,30,50,80,100)

peakLocations_c1_50 <- seq(800, 1800, 50 )

mu1_c1 <- median(peakScale_c1)
sigma1_c1 <- sd(log(peakScale_c1))
mu2_c1 <- median(peakAmplitude_c1)
sigma2_c1 <- sd(peakAmplitude_c1)

lPriors_c1_50 <- list(scale.mu=log(mu1_c1) - (sigma1_c1^2)/2, scale.sd=sigma1_c1, bl.smooth=10^11, bl.knots=30,
                      beta.mu=mu2_c1, beta.sd=sigma2_c1 , noise.sd=0.1, noise.nu=4)

tm_c1_50 <- system.time(result_c1_50 <- fitSpectraSMC(shift_c1, raw_spectra_c1, peakLocations_c1_50, lPriors_c1_50))
bl.est_c1_50 <- result_c1_50$basis %*% result_c1_50$alpha[,1,] 
b_c1_50 <- apply( bl.est_c1_50, 1, mean) 
b_c1_50 <- as.matrix(b_c1_50) 
expFn_c1_50 <- apply(result_c1_50$expFn,2,mean)

raw_spectra_mean_c1 <- apply(raw_spectra_c1, 2, mean)
plot(shift_c1, raw_spectra_mean_c1, type='l',xlab=expression(paste("Raman shift (cm"^{-1}, ")")), ylab="Intensity (a.u.)", ylim = c(-1,1.2), main = 'the Spectrum for Cyclohexane-1 when knots=10' )
lines(shift_c1, b_c1_50, col="red") #baseline
lines(shift_c1, b_c1_50 + expFn_c1_50,col="blue") #spectral signature
legend("bottomright", inset = .05, c("baseline","baseline + spectral signature"),lty = c(1,1), col = c("red", "blue"))
abline(v = peakLocations_c1_50, col = "green", lty = 3)


ci_hdi_c1_50 <- ci(result_c1_50$beta, method = "HDI")
raw_spectra_mean_c1 <- apply(raw_spectra_c1, 2, mean)
raw_data <- data.frame(shift_c1, raw = raw_spectra_mean_c1)
baseline_c1 <- data.frame(shift_c1, b_c1_50)
signature_c1 <- data.frame( shift_c1, sig = b_c1_50 + expFn_c1_50)


p_c1_50 <- ggplot()+
  geom_line(data = raw_data, aes(x = shift_c1, y = raw))+
  geom_ribbon(data = raw_data, aes(ymin = raw-ci_hdi_c1_50$CI_low , ymax = raw + ci_hdi_c1_50$CI_low, x = shift_c1, alpha = 0.01), col = "grey", fill = "grey")+
  geom_line(data = signature_c1, aes(x = shift_c1, y = sig), color = "royalblue")+
  geom_ribbon(data = signature_c1, aes(ymin = sig - ci_hdi_c1_50$CI_low, ymax = sig + ci_hdi_c1_50$CI_low, x = shift_c1, alpha = 0.01), col = "royalblue", fill = "royalblue")+  
  labs(x = expression(paste("Raman shift (cm"^{-1}, ")")), y = "Intensity (a.u.)")

p_c1_50+axis_theme
