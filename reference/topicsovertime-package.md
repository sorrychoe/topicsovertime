# topicsovertime: Topics over Time: A Continuous-Time Topic Model

Fits the Topics over Time model of Wang and McCallum (2006)
[doi:10.1145/1150402.1150450](https://doi.org/10.1145/1150402.1150450) ,
a latent Dirichlet allocation style topic model in which each topic also
has a Beta distribution over continuous time, so topics are shaped by
both word co-occurrence and document timestamps without discretizing
time or assuming Markov dynamics. Inference uses collapsed Gibbs
sampling implemented in 'C++'. Includes prediction of topic proportions
and timestamps for new documents, topic trends over time, and simulation
of synthetic corpora.

## See also

Useful links:

- <https://sorrychoe.github.io/topicsovertime/>

- <https://github.com/sorrychoe/topicsovertime>

- Report bugs at <https://github.com/sorrychoe/topicsovertime/issues>

## Author

**Maintainer**: Jaeseong Choe <cjssoote@gmail.com> \[copyright holder\]

Authors:

- Jaeseong Choe <cjssoote@gmail.com> \[copyright holder\]
