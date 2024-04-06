# bayesian-size-spectrum
A Bayesian approach to parameterizing size spectrum models using electrofishing stream fish data (Richter et al., In Review). This approach is modified from Edwards et al. (2017) and focuses on using fish data to fit a Pareto distribution with binned size data to estimate the slope parameter of the size spectrum model. As this approach is primarily focused on stream fish assemblages, we incorporate a correction to the capture probability of stream fish due to the size selectivity that is associated with stream fish electrofishing samples (Richter et al. 2022). The associated script for this model was written in R and uses the JAGS program (https://sourceforge.net/projects/mcmc-jags/) to implement the Bayesian sampler. We provide simulated fish catch data that can be used to parameterize a size spectrum model using the script. 

Model steps:
- Assign all fish to size classes (log2-scale)
- Determine the lower and upper limits of each size bin along with the size bin midpoints
- Calculate the total catch for each size bin
- Use the size selectivity correction model to calculate the capture probability for each size bin
- Define input data for the model
- Input parameters for JAGS program
- Create model file (sizespectrum.model.jags)
- Run the JAGS model
- Save model output

## Data
This repository includes model parameters for the size selectivity correction (b0.summary.csv, beta.summary.csv) and simulated stream fish data (samplefishdata.csv). The fish data includes body size measurements, both length and weight, of all fish captured at the site. 

## Model structure
The fish data is used to assign each captured fish to size bins (log~2~-scale). The model assumes that the size bin catches $\bf{c}$ follow a multinomial distribution ($\bf{c} \sim \text{Multinomial}(n, \bf{p})$) where $n$ represents the total catch and number of independent trials and $\bf{p}$ is the probability of observing a fish in each size bin. The probability of sampling a fish from size bin $b$ is the capture probability for that size bin ($q_{b}$) multiplied by the Pareto probability density function integrated across the size bin. The sampling probability is then standardized across all values to ensure their sum equalled 1 for the multinomial distribution. 
Pb formula
This fitted Pareto distribution is essentially a truncated power law distribution and the focus of this model is to parameterize the shape parameter ($\alpha$) of the Pareto distribution, which is directly related to the slope of the size spectrum $\lambda$ ($\lambda$ = ($\alpha + 1)$). The location parameter was assigned the value of the minimum body size used for analysis (1 g).

For more details about the structure of the model, please refer to Richter et al. (In Review).

## References
Edwards, A. M., Robinson, J. P., Plank, M. J., Baum, J. K., & Blanchard, J. L. (2017). Testing and recommending methods for fitting size spectra to data. Methods in Ecology and Evolution, 8(1), 57-67. https://doi.org/10.1111/2041-210X.12641  
Richter, I. A., Giacomini, H. C., De Kerckhove, D. T., Jackson, D. A., & Jones, N. E. (2022). Correcting for size bias in electrofishing removal samples. Ecological Modelling, 467, 109929. https://doi.org/10.1016/j.ecolmodel.2022.109929  
Richter, I. A., Giacomini, H. C., De Kerckhove, D. T., Jones, N. E., & Jackson, D. A. In Review. From individual sites to the entire watershed: variability in size spectrum models among stream fish communities.![image](https://github.com/ianarichter/bayesian-size-spectrum/assets/52668791/72da9db0-e728-4171-bcb5-7dafc04d42ea)
