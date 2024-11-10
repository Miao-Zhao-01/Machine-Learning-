# first, we employ the Bayesian model for our collected data: Cyclohexane, Sesame Oil, and Lung Tissue
# Lung Tissue
setwd("/Users/zhaomiao/Downloads/dissertation/rawdata/WMRS_785nm_10measurements_1s/Tissue-1")

for(j in 11:20){
  i <- formatC(j, width = 3, flag = 0) #turn 1 into 001
  path <- paste0("Raman_QEP001961_",i, ".txt") 
  d <- read.csv(file = path, header = TRUE) 
  d = d[-c(1:12),]
  delimiter <- "\t" # Replace with the appropriate delimiter used in file
  d <- strsplit(d, delimiter)
  raman <- data.frame(matrix(unlist(d), nrow = length(d), byrow = TRUE))
  raman[,2] <- as.numeric(raman[,2])
  raman_spectra_t1 <- matrix(0, 1,1044)
  raman_spectra_t1[1,] <- as.matrix(raman[,2])
  max_spectra <- max(raman_spectra_t1)
  raman_spectra_t1 <- raman_spectra_t1 / max_spectra
  wavelength_t1 <- as.numeric(unlist(raman[,1]))
  
  if (j==11){
    raw_spectra_t1 = raman_spectra_t1
  }
  else{
    raw_spectra_t1 = rbind(raw_spectra_t1, raman_spectra_t1)
  }
}

raman_2rows <- matrix(0, 2,1044)
raman_2rows[1,] <- apply(raw_spectra_t1, 2, sum)
max<- max(raman_2rows[1,])
raman_2rows[1,] <- raman_2rows[1,]/max

# read another data
setwd("/Users/zhaomiao/Downloads/dissertation/rawdata/WMRS_785nm_1measurement_10s/Tissue-1")
j <- 1
i <- formatC(j, width = 3, flag = 0) 
path <- paste0("Raman_QEP001961_",i, ".txt") 
d <- read.csv(file = path, header = TRUE) 
d = d[-c(1:12),]
delimiter <- "\t" # Replace with the appropriate delimiter used in file
d <- strsplit(d, delimiter)
raman <- data.frame(matrix(unlist(d), nrow = length(d), byrow = TRUE))
raman[,2] <- as.numeric(raman[,2])
raman_2rows[2,] <- raman[,2]
max <- max(raman_2rows[2,])
raman_2rows[2,] <- raman_2rows[2,]/max

exc_wavelength <- 785.423297
shift_t1 <- ( ( 1 / exc_wavelength) - (1/ wavelength_t1)) * 1e7

peakLocations_t1 <- c(1131, 1302, 1449,1555, 1624, 1733, 2333)
peakScale_t1 <- c( 9, 64, 36, 12, 17, 4, 10)
peakAmplitude_t1 <- c(0.017, 0.06 , 0.076, 1,0.095, 0.017, 0.2)

mu1_t1 <- median(peakScale_t1)
sigma1_t1 <- sd(log(peakScale_t1))
mu2_t1 <- median(peakAmplitude_t1)
sigma2_t1 <- sd(peakAmplitude_t1)

lPriors_t1 <- list(scale.mu=log(mu1_t1) - (sigma1_t1^2)/2, scale.sd=sigma1_t1, bl.smooth=10^11, bl.knots=30,
                   beta.mu=mu2_t1, beta.sd=sigma2_t1 , noise.sd=0.1, noise.nu=4)

tm_t1_2rows <- system.time(result_t1_2rows <- fitSpectraSMC(shift_t1, raman_2rows, peakLocations_t1, lPriors_t1))
bl.est_t1_2rows <- result_t1_2rows$basis %*% result_t1_2rows$alpha[,1,] 
b_t1_2rows <- apply( bl.est_t1_2rows, 1, mean) 
b_t1_2rows <- as.matrix(b_t1_2rows) 
expFn_t1_2rows <- apply(result_t1_2rows$expFn,2,mean)

plot(shift_t1, raman_2rows[1,], col='blue',type='l', xlab=expression(paste("Raman shift (cm"^{-1}, ")")), ylab="Intensity (a.u.)", cex.axis=1.5,cex.lab=1.5 )
lines(shift_t1, b_t1_2rows, col="red") #baseline
lines(shift_t1, b_t1_2rows + expFn_t1_2rows,col="blue") #spectral signature
legend("topright", cex = 0.8,inset = .01, c("Raman Spectrum","peak locations"),lty = c(1,3), col = c("blue", "green"))
abline(v = peakLocations_t1, col = "green", lty = 3)

intensity_shift_t1 <- data.frame( raman_shift = round(shift_t1,0), intensity =  b_t1_2rows + expFn_t1_2rows )
for ( i in 1:length(peakLocations_t1)){
  j = peakLocations_t1[i]
  text(j, intensity_shift_t1$intensity[which(intensity_shift_t1$raman_shift == j)], j)
}

intensity_shift_t1 <- data.frame( raman_shift = ceiling(shift_t1), intensity =  b_t1_2rows+ expFn_t1_2rows )
intensity_shift_t1 <- data.frame( raman_shift = floor(shift_t1), intensity =  b_t1_2rows+ expFn_t1_2rows )


ci_hdi_t1 <- ci(result_t1_2rows$beta, method = "HDI")

raw_data <- data.frame(shift_t1, raw = raman_2rows[1,])
baseline_t1 <- data.frame(shift_t1, b_t1_2rows)
signature_t1 <- data.frame( shift_t1, sig = b_t1_2rows + expFn_t1_2rows)
sigma_range_t1 <- sd(signature_t1$sig)

axis_theme <- theme(
  axis.title = element_text(
    size = 20
  )
)

p_t1 <- ggplot()+
  geom_line(data = raw_data, aes(x = shift_t1, y = raw))+
  geom_ribbon(data = raw_data, aes(ymin = raw-ci_hdi_t1$CI_low , ymax = raw + ci_hdi_t1$CI_low, x = shift_t1, alpha = 0.01), col = "grey", fill = "grey")+
  geom_line(data = signature_t1, aes(x = shift_t1, y = sig), color = "royalblue")+
  geom_ribbon(data = signature_t1, aes(ymin = sig - ci_hdi_t1$CI_low, ymax = sig + ci_hdi_t1$CI_low, x = shift_t1, alpha = 0.01), col = "blue", fill = "blue")+  
  labs(x = expression(paste("Raman shift (cm"^{-1}, ")")), y = "Intensity (a.u.)")

p_t1+axis_theme
#----------------------------------------------------------------------------------
# Sesame Oil
setwd("/Users/zhaomiao/Downloads/dissertation/rawdata/WMRS_785nm_10measurements_1s/Sesame_oil-1")

for(j in 11:20){
  i <- formatC(j, width = 3, flag = 0) #turn 1 into 001
  path <- paste0("Raman_QEP001961_",i, ".txt") 
  d <- read.csv(file = path, header = TRUE) 
  d = d[-c(1:12),]
  delimiter <- "\t" # Replace with the appropriate delimiter used in file
  d <- strsplit(d, delimiter)
  raman <- data.frame(matrix(unlist(d), nrow = length(d), byrow = TRUE))
  raman[,2] <- as.numeric(raman[,2])
  raman_spectra_s1 <- matrix(0, 1,1044)
  raman_spectra_s1[1,] <- as.matrix(raman[,2])
  max_spectra <- max(raman_spectra_s1)
  raman_spectra_s1 <- raman_spectra_s1 / max_spectra
  wavelength_s1 <- as.numeric(unlist(raman[,1]))
  
  if (j==11){
    raw_spectra_s1 = raman_spectra_s1
  }
  else{
    raw_spectra_s1 = rbind(raw_spectra_s1, raman_spectra_s1)
  }
}

raman_2rows_s1 <- matrix(0, 2,1044)
raman_2rows_s1[1,] <- apply(raw_spectra_s1, 2, sum)
max<- max(raman_2rows_s1[1,])
raman_2rows_s1[1,] <- raman_2rows_s1[1,]/max

# read another data
setwd("/Users/zhaomiao/Downloads/dissertation/rawdata/WMRS_785nm_1measurement_10s/Sesame_oil-1")
j <- 1
i <- formatC(j, width = 3, flag = 0) 
path <- paste0("Raman_QEP001961_",i, ".txt") 
d <- read.csv(file = path, header = TRUE) 
d = d[-c(1:12),]
delimiter <- "\t" # Replace with the appropriate delimiter used in file
d <- strsplit(d, delimiter)
raman <- data.frame(matrix(unlist(d), nrow = length(d), byrow = TRUE))
raman[,2] <- as.numeric(raman[,2])
raman_2rows_s1[2,] <- raman[,2]
max <- max(raman_2rows_s1[2,])
raman_2rows_s1[2,] <- raman_2rows_s1[2,]/max

exc_wavelength <- 785.423297
shift_s1 <- ( ( 1 / exc_wavelength) - (1/ wavelength_t1)) * 1e7

peakLocations_s1 <- c(1131, 1302, 1449,1555, 1624, 1733, 2333)
peakScale_s1 <- c( 9, 64, 36, 12, 17, 4, 10)
peakAmplitude_s1 <- c(0.017, 0.06 , 0.076, 1,0.095, 0.017, 0.2)


mu1_s1 <- median(peakScale_s1)
sigma1_s1 <- sd(log(peakScale_s1))
mu2_s1 <- median(peakAmplitude_s1)
sigma2_s1 <- sd(peakAmplitude_s1)

lPriors_s1 <- list(scale.mu=log(mu1_s1) - (sigma1_s1^2)/2, scale.sd=sigma1_s1, bl.smooth=0.5^11, bl.knots=30,
                   beta.mu=mu2_s1, beta.sd=sigma2_s1 , noise.sd=0.1, noise.nu=4)

tm_s1_2rows <- system.time(result_s1_2rows <- fitSpectraSMC(shift_s1, raman_2rows_s1, peakLocations_s1, lPriors_s1))
bl.est_s1_2rows <- result_s1_2rows$basis %*% result_s1_2rows$alpha[,1,] 
b_s1_2rows <- apply( bl.est_s1_2rows, 1, mean) 
b_s1_2rows <- as.matrix(b_s1_2rows) 
expFn_s1_2rows <- apply(result_s1_2rows$expFn,2,mean)

plot(shift_s1, raman_2rows_s1[1,], col='blue',type='l', xlab=expression(paste("Raman shift (cm"^{-1}, ")")),ylim=c(0,1.1), ylab="Intensity (a.u.)", cex.axis=1.5,cex.lab=1.5 )
lines(shift_s1, b_s1_2rows, col="red") #baseline
lines(shift_s1, b_s1_2rows + expFn_s1_2rows,col="blue") #spectral signature
legend("topleft", cex = 0.8,inset = .01, c("Raman Spectrum","peak locations"),lty = c(1,3), col = c("blue", "green"))
abline(v = peakLocations_s1, col = "green", lty = 3)

intensity_shift_s1 <- data.frame( raman_shift = round(shift_s1,0), intensity =  b_s1_2rows + expFn_s1_2rows )
for ( i in 1:length(peakLocations_s1)){
  j = peakLocations_s1[i]
  text(j, intensity_shift_s1$intensity[which(intensity_shift_s1$raman_shift == j)], j)
}

intensity_shift_s1 <- data.frame( raman_shift = ceiling(shift_s1), intensity =  b_s1_2rows+ expFn_s1_2rows )
intensity_shift_s1 <- data.frame( raman_shift = floor(shift_s1), intensity =  b_s1_2rows+ expFn_s1_2rows )


ci_hdi_s1 <- ci(result_s1_2rows$beta, method = "HDI")

raw_data <- data.frame(shift_s1, raw = raman_2rows_s1[1,])
baseline_s1 <- data.frame(shift_s1, b_s1_2rows)
signature_s1 <- data.frame( shift_s1, sig = b_s1_2rows + expFn_s1_2rows)


p_s1 <- ggplot()+
  geom_line(data = raw_data, aes(x = shift_s1, y = raw))+
  geom_ribbon(data = raw_data, aes(ymin = raw-ci_hdi_s1$CI_low , ymax = raw + ci_hdi_s1$CI_low, x = shift_s1, alpha = 0.01), col = "grey", fill = "grey")+
  geom_line(data = signature_s1, aes(x = shift_s1, y = sig), color = "royalblue")+
  geom_ribbon(data = signature_s1, aes(ymin = sig - ci_hdi_s1$CI_low, ymax = sig + ci_hdi_s1$CI_low, x = shift_s1, alpha = 0.01), col = "royalblue", fill = "royalblue")+  
  labs(x = expression(paste("Raman shift (cm"^{-1}, ")")), y = "Intensity (a.u.)")

p_s1+axis_theme

#----------------------------------------------------------------------------------
# Cyclohexane
setwd("/Users/zhaomiao/Downloads/dissertation/rawdata/WMRS_785nm_10measurements_1s/Cyclohexane-1")

for(j in 11:20){
  i <- formatC(j, width = 3, flag = 0) #turn 1 into 001
  path <- paste0("Raman_QEP001961_",i, ".txt") 
  d <- read.csv(file = path, header = TRUE) 
  d = d[-c(1:12),]
  delimiter <- "\t" # Replace with the appropriate delimiter used in file
  d <- strsplit(d, delimiter)
  raman <- data.frame(matrix(unlist(d), nrow = length(d), byrow = TRUE))
  raman[,2] <- as.numeric(raman[,2])
  raman_spectra_c1 <- matrix(0, 1,1044)
  raman_spectra_c1[1,] <- as.matrix(raman[,2])
  max_spectra <- max(raman_spectra_c1)
  raman_spectra_c1 <- raman_spectra_c1 / max_spectra
  wavelength_c1 <- as.numeric(unlist(raman[,1]))
  
  if (j==11){
    raw_spectra_c1 = raman_spectra_c1
  }
  else{
    raw_spectra_c1 = rbind(raw_spectra_c1, raman_spectra_c1)
  }
}

raman_2rows_c1 <- matrix(0, 2,1044)
raman_2rows_c1[1,] <- apply(raw_spectra_c1, 2, sum)
max<- max(raman_2rows_c1[1,])
raman_2rows_c1[1,] <- raman_2rows_c1[1,]/max

# read another data
setwd("/Users/zhaomiao/Downloads/dissertation/rawdata/WMRS_785nm_1measurement_10s/Cyclohexane-1")
j <- 1
i <- formatC(j, width = 3, flag = 0) 
path <- paste0("Raman_QEP001961_",i, ".txt") 
d <- read.csv(file = path, header = TRUE) 
d = d[-c(1:12),]
delimiter <- "\t" # Replace with the appropriate delimiter used in file
d <- strsplit(d, delimiter)
raman <- data.frame(matrix(unlist(d), nrow = length(d), byrow = TRUE))
raman[,2] <- as.numeric(raman[,2])
raman_2rows_c1[2,] <- raman[,2]
max <- max(raman_2rows_c1[2,])
raman_2rows_c1[2,] <- raman_2rows_c1[2,]/max

exc_wavelength <- 785.423297
shift_c1 <- ( ( 1 / exc_wavelength) - (1/ wavelength_c1)) * 1e7

peakLocations_c1 <- c( 1030, 1159, 1268, 1350, 1444,1555,1622,2333)
peakAmplitude_c1 <- c(0.434, 0.094, 0.454, 0.075, 0.536, 0.501, 0.048, 0.2)
peakScale_c1 <- c(17, 12, 17, 12, 18, 9, 17,10)

mu1_c1 <- median(peakScale_c1)
sigma1_c1 <- sd(log(peakScale_c1))
mu2_c1 <- median(peakAmplitude_c1)
sigma2_c1 <- sd(peakAmplitude_c1)

lPriors_c1 <- list(scale.mu=log(mu1_c1) - (sigma1_c1^2)/2, scale.sd=sigma1_c1, bl.smooth=10^11, bl.knots=30,
                   beta.mu=mu2_c1, beta.sd=sigma2_c1 , noise.sd=0.1, noise.nu=4)

tm_c1_2rows <- system.time(result_c1_2rows <- fitSpectraSMC(shift_c1, raman_2rows_c1, peakLocations_c1, lPriors_c1))
bl.est_c1_2rows <- result_c1_2rows$basis %*% result_c1_2rows$alpha[,1,] 
b_c1_2rows <- apply( bl.est_c1_2rows, 1, mean) 
b_c1_2rows <- as.matrix(b_c1_2rows) 
expFn_c1_2rows <- apply(result_c1_2rows$expFn,2,mean)

plot(shift_c1, raman_2rows_c1[1,], col='blue',type='l', xlab=expression(paste("Raman shift (cm"^{-1}, ")")), ylab="Intensity (a.u.)", cex.axis=1.5,cex.lab=1.5 )
lines(shift_c1, b_c1_2rows, col="red") #baseline
lines(shift_c1, b_c1_2rows + expFn_c1_2rows,col="blue") #spectral signature
legend("topright", cex = 0.8,inset = .01, c("Raman Spectrum","peak locations"),lty = c(1,3), col = c("blue", "green"))
abline(v = peakLocations_c1, col = "green", lty = 3)


legend("topleft", cex = 0.8,text.font = 0.4,inset = .01, c("Raman Spectrum"),lty = 1, col="blue")

intensity_shift_c1 <- data.frame( raman_shift = round(shift_c1,0), intensity =  b_c1_2rows + expFn_c1_2rows )
for ( i in 1:length(peakLocations_c1)){
  j = peakLocations_c1[i]
  text(j, intensity_shift_c1$intensity[which(intensity_shift_c1$raman_shift == j)], j)
}

intensity_shift_c1 <- data.frame( raman_shift = ceiling(shift_c1), intensity =  b_c1_2rows + expFn_c1_2rows )
intensity_shift_c1 <- data.frame( raman_shift = floor(shift_c1), intensity =  b_c1_2rows + expFn_c1_2rows )


ci_hdi_c1 <- ci(result_c1_2rows$beta, method = "HDI")

raw_data <- data.frame(shift_c1, raw = raman_2rows_c1[1,])
baseline_c1 <- data.frame(shift_c1, b_c1_2rows)
signature_c1 <- data.frame( shift_c1, sig = b_c1_2rows + expFn_c1_2rows)


p_c1 <- ggplot()+
  geom_line(data = raw_data, aes(x = shift_c1, y = raw))+
  geom_ribbon(data = raw_data, aes(ymin = raw-ci_hdi_c1$CI_low , ymax = raw + ci_hdi_c1$CI_low, x = shift_c1, alpha = 0.01), col = "grey", fill = "grey")+
  geom_line(data = signature_c1, aes(x = shift_c1, y = sig), color = "royalblue")+
  geom_ribbon(data = signature_c1, aes(ymin = sig - ci_hdi_c1$CI_low, ymax = sig + ci_hdi_c1$CI_low, x = shift_c1, alpha = 0.01), col = "royalblue", fill = "royalblue")+  
  labs(x = expression(paste("Raman shift (cm"^{-1}, ")")), y = "Intensity (a.u.)")

p_c1+axis_theme