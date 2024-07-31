# GaussianCdfPdfProduct
This simulation helps to visualize posterior probability distributions that are proportional to the product of a Gaussian PDF and CDF:

$$p\left( x \right) \propto \sigma_1^{-1} \ \varphi\left( \frac{x - \mu_1}{\sigma_1} \right) \ \Phi \left( \frac{x - \mu_2}{\sigma_2} \right)$$

![](Misc/animation.gif)

We approximate the mean and variance by performing Newton's method once from both means, picking whichever evaluates to larger posterior density.  
