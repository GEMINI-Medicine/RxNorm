# RxNorm
The GEMINI-RxNorm code automates the use of multiple RxNorm tools in tandem with other datasets to identify drug concepts from pharmacy orders. 

Please note rxnorm **match_score** in the rxnorm API **has changed** (<https://lhncbc.nlm.nih.gov/RxNav/news/RxNormApproxMatch.html>) since the writing of the manuscript and has not been validated by GEMINI. GEMINI recommends manual review of all rxnorm results without filtering based on match score threshold.

# Installation

The script should be saved as a local file, and `source` within your R script. The current function uses the GEMINI data holdings. In order for this function to work with other data holdings, the database connection, and database name should be changed to customize the given database. 

# Our Paper on GEMINI-RxNorm tool

<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10409892/>
