This repository contains the files needed for deploying the shiny version of the Manhattan++ plot.

To setup please perform the following steps.
1) Install the latest version of R (https://www.r-project.org/) and R-studio (https://posit.co/download/rstudio-desktop/). You also need to configure Git to work with Rstudio.
2) Create a new Git Version Control project in R-studio with this repository.
3) Use the script genfiles/genmatrix.R and genfiles/run.R to generate the files required by the shinyApp. You need to pass a GWAS results file with the following columns: maf, chr, pvalue, pos, conseq.
4) Modify the gwas_table.txt to add a GWAS to the App.
5) Test the shiny App by selecting run on the App.R files.
6) Publish the app to your shiny account.

If you use this repository please cite this paper:
Grace et al
Manhattan++: displaying genome-wide association summary statistics with multiple annotation layers
BMC Bioinformatics 2019; 20(1):610
PubMed: 31775616
