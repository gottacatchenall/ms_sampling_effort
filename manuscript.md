---
bibliography: [references.bib]
---

# Introduction


Ecological interactions are hard to sample [@Jordano2016SamNet]. Still,
collecting data about species interactions is imperative to measure and mitigate
the effects of human activity on Earth's biodiversity [@Andys90AuthorPaper;
@Jordano2016ChaEco], and to predict potential spillover of zoonotic diseases
[@cite]. Over the past decade biodiversity data has become increasingly
available, both due to improved sensing technology [@Stephenson2020TecAdv] but
also growth of open and FAIR databases. Remote sensing has enabled data on
spatial scales previously unsampleable, and in-situ observations in the form of
both cameras and environmental sensors have greatly improved the resolution of
data. Yet sampling of ecological interactions detection often requires human
sampling as coexistance is not necessarily indicative of interaction
[@Blanchet2020CooNot].

This induces constraints on sampling of interactions based on the spatial and
temporal scales feasible to human sampling. These sampling constraints go on to
bias species interaction data: we only observe but a small fraction of the
variance in species interactions in space and time, and these observations
themselves reflect the distribution of abundance within communities.
[@Poisot2015SpeWhy]. Further sampling of species interactions is geographically
biased toward the usual suspects [@Poisot2021GloKno]. These biases the data we.
collect is noisy and likely contains many false-negatives. This has many
practical consequences for answering questions about species interactions and
how human activity is effecting them [@deAguiar2019RevBia].

In this manuscript we seek to determine how false negatives impact analysis and
prediction of ecological networks, and how understanding the relationship
between sampling effort and probability of a true negative can guide how we
design surveys of ecological interactions [@Jordano2016SamNet]. The fundemental
questions we seek to answer are: 1) How many times do you have to observe a
non-interaction between two species to be confident in saying that is a true
negative? 2) How "wrong" are the measurements of network structure modularity as
a function of false-negative probability? 3) How do false-negatives impact our
ability to make reliable predictions about interactions?
We conclude by arguing for using null models such as those used  

# How many observations of a non-interaction do we need to classify it as a true negative?

A naive model of interaction detection would assume that every true interaction
between two species is incorrectly observed as a non-interaction with an
independent and fixed probability, which we denote $p_{fn}$ and subsequently
refer to as the False-Negative Rate (FNR). In this model, if we observe the same
species not-interacting $O$ times the probability of a true negative, $p_{tn}$,
is given by $p_{tn} = 1 - (p_{fn})^O$. This relation is shown in
@fig:negativebinom for varying values of the false negative rate $p_{fn}$. This
illustrates a fundamental link between our ability to reliably say an
interaction doesn't exist---$p_{tn}$---and our sampling effort $O$.

![The probability an observed interaction is a "true negative" (y-axis) given
how many times it has been sampled as a non-interaction. Another thing is true,
which is that this function will never reach 1. So many assumptions here about
probability. It's the birthday paradox, but backwards.
](./figures/negativebinom.png){#fig:negativebinom}

From @fig:negativebinom it is evident that the more times we see two species
_present_ but _not interacting_, the more likely the interaction is a true
negative. But what should this threshold of number of observations be?





Observations of species occur according to their relative abundance, and
this can lead to high realized values of $p_{fn}$ for species with low relative
abundance.

# False-negatives as a product of relative abundance

In this section we demonstrate the realized probability of false-negative
changes drastically with sampling effort simply as a function of the
distribution of species abundances within a community.
We do this by simulating the observation process on both NN empirical food webs
from the Mangal database [@Banville2021ManJl] and food-webs generated using
the niche model.

A simple model of observation assumes each observed species is drawn from this
abundance distribution. Seeing two low biomass species interacting requires two
relatively low probability events, which is detecting each species of low
biomass. Generally across communities, the shape of this abundance distribution
can be reasonably-well described by a log-normal distribution
[@Volkov2003NeuThe]. Controversies around theory of species abundance
distributions and neutral theory aside,

For simplicity we simulate abundances from $N_S$ independent draws from a
standard-log-normal distribution. For an ecological network $A$ with $N_S$
species, for each true interaction $A_{ij} = 1$ we estimate the probability of
observing both species $i$ and $j$ at given place and time by simulating a
distribution of $O$ observations, where the species observed at the
$1,2,\dots,O$-th observation is drawn from the abundance distribution. If both
$i$ and $j$ are present in the $O$ observations, the observation is computed as
a true-negative, and if not, as a false-negative. The results for applying this
can be found in @fig:samplingeffort, where the left panel is
this model applied NUM empirical food-webs from the Mangal database [@] , and
on the right the niche model [@].

![False negative rate as a function of sampling effort and network size,
computed using the method described above. Left panel:  in blue. Right empirical
food webs from Mangal database in teal. The outlier on panel B is a 714 species
network. ](./figures/samplingdist.png){#fig:samplingeffort}

Empirical data on interactions, limited by the practical realities of funding
and human-work hours, tends to fall on the order on 100s [@JordanoTable1].
Yet species richness clearly effects this and should be taken into account when
designing samples.


This simple model doesn't consider the possibility that there are positive or
negative associations between observing two species together based their
interaction. Here we assume each individual observation of a given single
species $i$ within a species pool occurs according to the distribution the
abundances of the species in that species pool .

However, we can demonstrate that the probability of observing a false negative
is _higher_ if there is some positive association between occurrence of species
$A$ and $B$. We can express the probability that we observe an existing
interaction between as $P(AB)$. If there is no correlation between probability
of observing $A$ and observing $B$, then the probability of observing the
interaction $P(AB) = P(A)P(B)$. In this case, the probability of observing both
$A$ and $B$, which we denote $P(AB)$, is not equal to $P(A)P(B)$ as $P(A)$ and
$P(B)$ are not independent. If there some positive strength of association
between observing both $A$ and $B$ because this interaction is "important" for
each species, then

$$P(AB) > P(A)P(B)$$

In this case, the probability of observing a false negative is $p_{fn} = 1 -
P(AB)$ which due to the above inequality due to positive associated implies

$$p_{fn} \geq 1 - P(A)P(B)$$


Caveats: this doesn't consider variation in abundance in space and time which is
kind of a problem. In this model every observation is drawn from the
distribution of the biomass distribution at a particular place and time. We
assume that this distribution is the same everywhere (again unlikely).


We now transition toward assessing the effects of false negatives in our
data on the properties derived from these measurements, and for use as data for
predicting interactions in the future. What levels of false negatives are
acceptable to infer network properties and predict interactions.


# Effects of false-negatives on network properties

Here we simulate the effects of observation error to generate false negatives in
the samples of ecological networks and compare the computed network properties
of the original network to the computed properties on the observed network in
order to see how false negatives effect our quantification of network structure.

![fig. 1$\sigma$ in first grad, 2$\sigma$ in second ](./figures/properties_error.png){#fig:properties}

# Effects of false negatives on ability to make predictions

In this section, we assess the effect of false negatives in data on our ability
to make predictions about interactions.

We use the predictive model and dataset as in @Strydom2021RoaPre to predict
interactions between species never observed at the same place and time.

This dataset from @Hadfield2014TalTwo describes host-parasite interaction
networks sampled across 51 sites. We partition the data into 80-20 training-test
split, and then seed the training data with false negatives varying rates, but
crucially do nothing to the test data.

The model---a neural-network with 3 layers to predict outputs based on features
extracted from cooccurence, see @Strydom2021RoaPre for more details).

In @fig:rocpr, we show receiving-operating-characteristic (ROC) and
precision-recall (PR) curves for the model with varying levels of
false-negatives added to the data.

![fig](./figures/rocpr_falsenegatives.png){#fig:rocpr}

Big takeaway here is false-negatives have more effect on PR space,
unsurprisingly. Sadly this is also where the potential application of is
greatest. Still, performance doesn't matter with many added false-negatives!
Good evidence in favor of this type of model. Same caveat as previous section
that this is data that _already_ is likely to have many false-negatives. So, the
effects of adding more, as we do in this illustration, might be mitigated.

# Conclusion


***The primary takeaways from this paper***
In this paper we have demonstrated that false negatives are likely purely due
to the distribution of.


***The primary recommendations for study design that this paper provides***
Take species richness and relative abundance into account. A model similar
to that which we show here can be used to provide a neutral expectation of
true-negative probability given a number of observations of individuals at
a given place and time.

***What does the future hold for this research***
A brief note on false positives.
How does this influence our understanding of the structure of ecological
networks, and how we infer other things based on that?
How does this influence our models of interaction prediction?
How does this effect how we design samples of interactions?
How can we correct for this bias in existing data?

# References
