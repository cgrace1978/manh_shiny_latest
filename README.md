This repository contains the files needed for deploying the shiny version of the Manhattan++ plot.

To setup please perform the following steps.
1) Install the latest version of R (https://www.r-project.org/), R-studio (https://posit.co/download/rstudio-desktop/) and Git (https://git-scm.com/downloads).
2) Configure R Studio to use Git
3) Create a new Git Version Control project in R-studio with this repository (https://github.com/cgrace1978/manh_shiny_latest).
4) Use the scripts **genfiles/genmatrix.R** and **genfiles/run.R** to generate the files required by the Manh++ shinyApp. You need create a GWAS results file which has the following columns:
   * **snp** - Identifier for the variant (rsid etc.)
   * **chr** - Chromosome
   * **pos** - Position
   * **pvalue** - P-value for the variant
   * **conseq** - High Impact Flag - mark as 1 if the variant is annotated as High Impact else 0.
   * **maf** - Minor allele frequency of the variant.
   * **genes** - Annotation regarding genes proximal to the variant.
6) Modify the **gwas_table.txt** file to add a GWAS to the App.
7) Test the shiny App by selecting run on the **App.R** files.
8) Publish the app to your shiny account.

If you use this repository please cite this paper:
Grace et al
Manhattan++: displaying genome-wide association summary statistics with multiple annotation layers
BMC Bioinformatics 2019; 20(1):610
PubMed: 31775616
