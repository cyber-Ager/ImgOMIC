# ImgOMIC

[EU] R-ko pakete hau nire doktoretzako proiektuan erabiliko ditudan funtzioak gordetzeko sortua da. Proiektu honek irudietatik ateratako ezaugarriak zientzia omikoekin (genomika, trankriptomika, proteomika etab.) aztertzea du helburu, hortik ImgOMIC izena.

[EN] This package was created with the intention of recopilating the neccessary functions to analyse medical imaging phenotypes with OMIC data. The package contains the neccessary function to perform these analyses from the data preprocessing to associative model building.

## Precision Calculation

In most cases when working with microarrays for the quantificatyion of the presence of miRNA, RNA or proteins duplicates are used. The duplicates ensure that we reduce the variability of the measure by computing the mean between the measures. But this also allows to quantify that variability and define the precision of the machine that we are using. Quantifying the precision allows us to know the variability of each measurement and it is important to take it in to account in order to ensure that the associations that we are identifying are real.

### Measure precision

This function is based on the 2009 article _"Estimating Precision Using Duplicate Measurements"_ by **Nicole Pauly Hyslop and Warren H. White** (https://www.tandfonline.com/doi/abs/10.3155/1047-3289.59.9.1032). 

The aim of this function is to quantify the precision of each of the measuring "sensors" used to detect the presence of the particle of interest. For example, in the case of miRNA detection by qPCR in a microarray, the "sensor" would correspond to the different miRNAs. Even if the same sample is introduced into each well, each well will emit fluorescent light based on the presence of a specific miRNA, as the "sensor mechanisms" (in this case, the stem-loop primers) differ for each miRNA. In this context, the variability of the measurements for each miRNA may vary, but this variability can be assessed using duplicates.

The function `measure_precision` measures the precission based on Root Mean Square (RMS) and Mean Absolute Difference (MAD). The equation used to measure this precission are the next one: 

$$
D_1 = \frac{(C_{i1}-C_{i2})/\sqrt{2}}{\bar{C_{i}}}
$$

$$
\text{RMS Precision} = \sqrt{\frac{1}{n}\sum_{i = 1}^{n}D_i^2}\times100\%
$$

$$
\text{MAD Precision} = \sqrt{\frac{\pi}{2}}\frac{1}{n}\sum_{i = 1}^{n}|D_i|\times100\%
$$
