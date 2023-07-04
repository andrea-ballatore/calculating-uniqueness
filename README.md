# A uniqueness index for multivariate observations

**Abstract**: The concept of uniqueness can play an important role when the assessment of an observation’s distinctiveness is essential. This article introduces a distance-based uniqueness measure that quantifies the relative rarity or commonness of a multi-variate observation within a dataset. Unique observations exhibit rare combinations of values, and not necessarily extreme values. Taking a cognitive psychological perspective, our measure defines uniqueness as the sum of distances between a target observation and all other observations. After presenting the measure u and its corresponding standardised version uz, we propose a method to calculate a p value through a probability density function. We then demonstrate the measure’s behaviour through a case study on the uniqueness of Greater London boroughs, based on real-world socioeconomic variables. This initial investigation indicates that u can support exploratory data analysis.

**Associated publication**: Ballatore, Andrea and Cavazzi, Stefano (2023) Why is Greenwich so common? Quantifying the uniqueness of multivariate observations. The 12th International Conference on Geographic Information Science (GIScience), Short Papers, Leeds, UK.

## Contents

The `uniqueness.Rmd` notebook contains the development code and analyis of the uniqueness index.

The core uniqueness functions are in file `uniq_functions.R`, which can be easily imported.