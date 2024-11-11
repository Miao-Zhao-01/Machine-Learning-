# Introduction
The project titled "Machine Learning for Spectroscopy in Biological Applications", is my dissertation work. It focuses on mitigating the influence of fluorescence intensity by using an open-source ssersBayes package within the R programming. Leveraging a Bayesian hierarchical model, we employ it to analyse our collected data. Additionally, experiments are conducted to evaluate the model’s fitting performance under variations in prior information. The processing  workflow of the model is illustrated below:
<div align=center>
<img width="596" alt="image" src="https://github.com/user-attachments/assets/782261ad-86aa-4d37-adc2-e8edf2bca52d">
</div>

# Data Set
1) Cyclohexane: Cyclohexane was chosen owing to its well-known Raman spectrum and absence of fluorescence background. 
2) Sesame oil: The Raman spectrum of sesame oil exhibits prominent, well-known peaks, but it is also characterized by a high level of fluorescence background. 
3) Lung tissue: This unique data was obtained from the Royal Infirmary of Edinburgh (NHS Lothian BioResource, reference 15/ES/0094). An ex vivo human lung tissue sample was collected from a patient diagnosed with suspected or confirmed lung cancer.

# Results
Through experiments, we successfully tackled the research issues and attained our stated goals. Notably, after modifying the model parameters, such as smoothness and the number of knots, our models show increased fitting impact.

Moreover, our models exhibit superior fitting capabilities for Cyclohexane and Lung Tissue compared to Sesame Oil, one example for lung tissue is shown in the following figure. This discrepancy in performance can be attributed to the fact that Sesame Oil data is noisier and more contaminated, which consequently impedes the model’s ability to accurately capture the peak amplitudes and scales.

A noteworthy aspect to highlight is that we only provided detailed peak locations, while amplitudes and scales are specified only by their median and standard derivation value. This also contributes to the challenge of fitting the signal signature. Despite the less-than-ideal conclusion for Sesame Oil, the outcomes obtained from Cyclohexane and Lung Tissue demonstrate how effective our models are at capturing peaks’ characteristics.

Additionally, we thoroughly investigated the influence of prior information on the model’s performance. We discovered that when peculiar values, such as identical values, zeros, or extracted values at regular intervals, were introduced as prior information, the fitting effectiveness of the models diminished significantly.

Conversely, we made an intriguing observation: when we supplied incorrect but plausible prior information to the models, they exhibited notably better fitting results. This finding highlights the models’ capacity to adapt to variations in prior data and underscores the importance of providing realistic and reasonable prior information for optimal performance.

<div align=center>
<img width="638" alt="image" src="https://github.com/user-attachments/assets/5910b51b-9d25-4fc0-ae6f-e6e52c258b2f">
</div>

# Main references
[1] Charles J Geyer. Introduction to markov chain monte carlo. Handbook of markov chain monte carlo, 20116022:45, 2011.  
[2] Nia C Jenkins, Katjana Ehrlich, Andr´as Kufcs´ak, Stephanos Yerolatsitis, Susan Fernandes, Irene Young, Katie Hamilton, Harry AC Wood, Tom Quinn, Vikki Young, et al. Computational fluorescence suppression in shifted excitation raman spectroscopy. IEEE Transactions on Biomedical Engineering, 2023.  
[3] Roger M Jarvis, Alan Brooker, and Royston Goodacre. Surface-enhanced raman spectroscopy for bacterial discrimination utilizing a scanning electron microscope with a raman spectroscopy interface. Analytical Chemistry, 76(17):5198–5202, 2004.  
[4] Andrzej Kudelski. Analytical applications of raman spectroscopy. Talanta, 76(1):1–8, 2008.  
[5] Jun S Liu and Jun S Liu. Monte Carlo strategies in scientific computing, volume 75. Springer, 2001.  
[6] Chris Sherlock, Paul Fearnhead, and Gareth O Roberts. The random walk metropolis: linking theory and practice through a case study. 2010.
