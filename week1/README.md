## Week 1 Homework: Linear Regression

# - Step 2.2

There is a detectable but weak positive relationship between maternal age and de Novo mutation (DNM) count, with a beta1 of 0.3776. The p-value of the model is 6.88e-24, meaning it is statistically significant. These conclusions are corroborated by the graph generated in Step 2.1 (ex2_b.png), as a gentle upward trend is visible across the data points, with a good deal of variation in de Novo mutation count per age. The r-squared value of the relationship is 0.228, which indicates maternal age only explains a small proportion (about 23%) of the variation seen in de Novo mutation frequency among sampled mothers. This too is consistent with understandings of biology, as  many factors contribute to mutation rate beyond age, such as environmental exposures, family history, and to some extent, simple random chance.

## Step 2.3 

For paternal de Novo mutations, the positive relationship between age and DNMs is stronger than in the maternal data, being 1.3538. The paternal model has a higher r-squared value (0.619), meaning the model explains about 62% of the variation in observed paternal de Novo mutations. Based on the p-value, the results are statistically significant with a p-value of 1.55e-84. When comparing graphs of the paternal and maternal dataset, the stronger positive relationship is evident as the slope are both greater in magnitude. Likewise, points are clearly less varied relative to the trend and to themselves. 

## Step 2.4

beta0 = y - beta1(x)

beta0 + beta1(x) = y 

10.3263 + 1.3538(50.5) = dnm

de Novo Mutations at 50 years, 6 months = 78.69 (79 dnms)

## Step 2.6

I used a t-test for related samples of scores in paternal and maternal datasets. I chose this test because the data are related in that each observation of de Novo Mutations occurs at the same age in the other dataset. From this test, I got a t-statistic of 61.61 and a p-value of 1.12e-204, meaning the data are statistically significantly different from one another. This means the number of DNMs in paternal samples is statistically different from the observed frequencies of maternal DNMs. This is evident from the data and the graphs, as the DNM counts are visually different. However, statistical testing is necessary to say with greater confidence that the data are truly statistically different. 

