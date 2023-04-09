# ML4G

## Group members
- Gudmundur Björgvin Magnusson 
- Hui Jeong Jung 

## Conda environment
- project1_base.yml: default conda environment
- project1.yml: custom R conda environment (to be submitted)

## Engineered features
We extracted the "scores" of the peaks from different experimentation methods. From DNase we extracted the 5th column as the score and for everything else we extracted the 7th column and further transformed these scores to create features to train our model. We also defined the regions "near" to the gene as the region 5kb before and after the TSS start site and end site. We defined the "far" regions as the region 20kb before and after. Below we have listed the list of different features that we have engineered for each gene that we are trying to predict the expression of. 

- gene length 
- sum of peak scores in the near and far region of the TSS 
- average of peak scores in the near and far region of the TSS 
- max of peak scores 
- number of peaks in the near and far region of the TSS
- number of peaks in the far region that are not near a TSS 
- max score of a peak in the far region that is not near a TSS 
- mean score of a peak in the far region that is not near a TSS
- distance to closest peak (signed and unsigned)
- score of the closest peak 
- score of the furthest peak 
- sum of the scores of peaks that are in the 10 "bins", or one subsection of the near region split into 10 folds 
- number of genes in the far neighborhood
- which strand the gene is located on (+ or - binarized)

