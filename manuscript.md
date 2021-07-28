---
bibliography: [references.bib]
---

# Introduction



As we enter a new generation of biodiversity monitoring, many forms of data are
increasingly available---remote sensing, . Ecological interactions are hard to
find [@Jordano2016ChaEco].

Yet sampling of ecological interactions detection requires in-situ observation.
This is subject to many biases:

Interactions vary in space and time [@Poisot2015SpeWhy], we are more likely to
observe interactions between species with high relative abundance [@cite]. As a
result of these biases, the data we collect is noisy and likely contains  many
false-negatives [@Poisot2021ImpMam]. This has many practical consequences for
answering questions about species interactions and how human activity is
effecting them.

In this manuscript we seek to answer:
How "wrong" are the measurements of network structure
(connectance, nestedness, modularity) as a function of false-negative probability?
How wrong are predictions about interactions?


![The probability an observed interaction is a "true negative" (y-axis) given
how many times it has been sampled as a non-interaction. Another thing is true,
which is that this function will never reach
1. So many assumptions here about probability.
](./figures/bernoulli.png){#fig:bernoulli}




# False-negatives as a product of relative abundance

Does a false negative rate of 0.9 seem unrealistic? Consider how the probability
of observation occurs as a function of abundance.

In this section we demonstrate the realised false-negative rate (FNR) is
different for high abundance vs low abundance species.

Consider a probability of false negative detection per unit biomass.
In this model every observation is drawn from the distribution of the biomass
distribution at a particular place and time. If we assume that this distribution
is the same everywhere (again unlikely)

Seeing two low biomass species interacting requires two relatively low prob
events, which is detecting each species of low biomass.

What if there is a strength of association? Covariance of biomass of i and
biomass of j due to cooccurence because this interaction is "important" for each
species.

  --> This implies that interactions that are variable/opportunistic are subject
  to ever higher false-negative rate.


# Effects of false-negatives on network properties

Here we simulate a bunch of food webs using generative models.
We then simulate the effects of observation error to generate
false negatives in the sample and compare the computed network
properties of the "true" networks to the computed properties on
the observed network in order to see how much false negatives
effect our quantification of network structure.


# Effects of false negatives on ability to make predictions

Use the same model and data as [@Strydom2021]. Seed the training
data with false negatives at a rate $p_{fn}$. Don't do anything to
the test data. Make ROC-PR AUC plots for 3 levels of $p_{fn}$. Same
model, same data, different levels of predictive capacity.


# Conclusion

How does this influence our models of interaction prediction?

How does this effect how we design samples of interactions?

How can we correct for this bias in existing data?

# References
