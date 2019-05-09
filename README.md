# SalTy-paper
This repo contains analysis and figure scripts for my paper on the phylogenetics of Salmonella Typhimurium.
Each folder contains the scripts and data necessary to produce the basic plots.

###Working Title:
**The population structure of Salmonella enterica serovar Typhimurium reveals signatures of anthropogenic selection and niche adaptation**

## Main Figures
### Figure 1
![Figure 1](Figure_1/Figure_1.pdf)
**The population structure of *Salmonella* Typhimurium.**

This figure also has a series of boxplot to show different measures of genome degredation.
![Figure 1 boxplots](Figure_1/boxplots/boxplots.pdf)
These are the number of matched HACS sequences. The mean deltabitscore and the invasivness index.
However, I was not happy with the matched HACS sequences or how the DBS plot looked so I used the mean DBS using SL1344 as a reference, and HACS as number of predicted loss of function genes. This led to the following:

![Figure 1 boxplots2](Figure_1/boxplots/boxplots-SL1344.pdf)

In this case the bottom row has 3rd level group 1 removed which leads to the ifference between the 1st level groups being significant.