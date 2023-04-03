# ML4G

## Conda environment
- project1_base.yml: default conda environment
- project1.yml: custom R conda environment (to be submitted)

## Engineered features
We extracted the "scores" of the peaks from different experimentation methods. From DNase we extracted the 5th column as the score and for everything else we extracted the 7th column and further transformed these scores to create features to train our model. We also defined the regions "near" to the gene as the region 5kb before and after the TSS start site and end site. We defined the "far" regions as the region 20kb before and after. Below we have listed the list of different features that we have engineered for each gene that we are trying to predict the expression of. 

- gene length 
- sum of peak scores in the near and far region of the TSS 
- average of peak scores 
- max of peak scores
- number of peaks in the near and far region of the TSS
- number of peaks in the far region that are not near a TSS 
- distance to closest peak (signed and unsigned)
- 