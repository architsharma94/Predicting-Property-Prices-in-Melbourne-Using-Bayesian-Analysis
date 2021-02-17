# Predicting-Property-Prices-in-Melbourne-Using-Bayesian-Analysis

## Repository guide:

* Data => Assignment2PropertyPrices.csv
* Analysis and Reporting => Analysis and Prediction Report.docx
* Code => code.R


## Introduction

Before COVID-19, Australia’s real estate was presumed as a constantly booming industry, for a significantly long period of time. As part of the development plans for its major cities, assessment and prediction of volatility in property prices has been a hot topic amongst data scientists ever since. The analysts have been putting in tremendous efforts to collect and examine property prices in the city of Melbourne, which is a major contributor in Australia’s economy.
The core procedure executed for this analysis is Markov Chain Monte Carlo (MCMC) Bayesian analysis, on a dataset comprising of area of property, number of bedrooms, bathrooms, carparks and type, along with corresponding property prices. The assessment was carried out to understand the statistical properties of the data collected and build a multiple linear model for predicting property prices in Melbourne.

## About Data

The dataset comprises of 10,000 observations for 5 predictors and 1 response variable (Y), i.e. sale price of Melbourne’s properties in 100,000 of AUD. The predictors are as follows:
1.	Area (X1): Land size in m2 of the sold property
2.	Bedrooms (X2): The number of bedrooms
3.	Bathrooms (X3): The number of bathrooms
4.	CarParks (X4): The number of car parks
5.	PropertyType (X5): The type of the property (0: House, 1: Unit) 

## Methodology

The analysis was carried out using Just Another Gibbs Sampler(JAGS) and R programming software. After obtaining the descriptive statistics of the data, the next steps involve the following:
* Creation of JAGS model diagram for the setting
  * Assess the descriptive statistics and graphical representations of variables
  * Determine the distribution of likelihood and prior
  * Visualize the distributions into a model diagram
* Specification of prior distributions as per the following expert knowledge provided in the problem statement
  * Area (X1): Every m2 increase in land size increases the sales price by 90 AUD. This is a very strong expert knowledge.
  * Bedrooms (X2): Every additional bedroom increases the sales price by 100,000 AUD. This is a weak expert knowledge.
  * Bathrooms (X3): There is no expert knowledge on the number of bathrooms.
  * CarParks (X4): Every additional car space increases the sales price by 120,000 AUD. This is a strong expert knowledge.
  * PropertyType (X5): If the property is a unit, the sale price will be 150,000 AUD less than that of a house on the average. This is a very strong expert knowledge.
* Compile the model using R programming and JAGS, and perform MCMC diagnostics to assess the effect of:
  * Number of chains
  * Burn-in period
  * Thinning
  * Number of saved steps
* After obtaining suitable diagnostics, visualize the posterior distribution for each of the variables and draw relevant inferences
* Obtain prediction model using the Bayesian point estimates
* Find predictions of sale prices for the information provided in problem statement

## Conclusion

Based on the data and expert knowledge provided for each of the 5 predictors, a multiple linear model was successfully simulated to predict property prices in Melbourne. The best diagnostics for the MCMC model were obtained using the following parameters:
* Number of Chains = 3
* Burn-in Period = 5000
* Thinning = 30
* Number of Saved Steps = 5000
* Adapt Steps = 1500

With the abovementioned configuration, satisfactory diagnostics were achieved for all beta priors, except for a few autocorrelations in intercept, type of property, number of bedrooms and bathrooms.
The β prior of area and type of property were found to be insignificant but were still considered for formulating the predictive model because their point estimates were extremely small. The final predictive model is given below:

Y = 4.97 + 2.54x10-6 xi1 + 0.11 xi2+ 0.42 xi3 + 0.109 xi4 + 0.0391 xi5

The goodness-of-fit of predictive model was examined by plotting a density plot of observed vs. predicted values. It was observed that the model failed to capture almost 50% of the high-density peak of observed values. This issue could be attributed to the following:
* Presence of possible outliers in data, which led to distortion of trend and subsequent misinterpretation of relationship between response and predictors
* Correlation between predictors, even though small, might have led to issues in modelling
* Considering beta priors of area and property type as significant, which could be handled by completely dropping them from the model
* Presence of a few autocorrelations in some beta priors’ diagnostics
