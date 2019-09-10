sigan
-------
A slightly modified version of the Broad's SignatureAnalyzer tool described in [10.1038/ng.3557](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4936490/).

## SignatureAnalyzer
SignatureAnalyzer is a program for detecting somatic mutational signatures.
It uses a Bayesian (Gibbs-sampling) NMF variant. It also provides handling for
hypermutator samples and both L1 and L2 regularization.


The original SignatureAnalyzer code is provided as-is in SignatureAnalyzer.R.

## sigan

sigan is identical to SignatureAnalyzer except:  
- it relies on the COSMIC v3 signatures for signature comparison
- It takes an output name for naming samples files.

Certain options are still hardcoded. At some point, it would be useful to make 
L1/L2 regularization and hypermutator handling available at the CLI. Right now, these
have to be edited manually in the script.

## Basic sage:

```
Rscript --vanilla ./sigan.R <MAF FILE> <"STUDY_NAME">
```

## Credit
Jaegil Kim wrote the original SignatureAnalyzer and maintains the copyright.
Eric Dawson modified the code to produce sigan.
