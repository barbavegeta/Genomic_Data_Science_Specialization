\# Statistics for Genomic Data Science (Course 6)



This directory contains my R scripts for the four module projects from \*\*Course 6 – Statistics for Genomic Data Science\*\* in the Johns Hopkins Genomic Data Science Specialization.



These scripts focus on applying core statistical methods to genomic data using \*\*R\*\*: exploratory analysis, hypothesis testing, simple models, and visualisation on teaching datasets.



---



\## Files



\- `module\_1\_project.R`  

&nbsp; R script for the \*\*module 1 project\*\*.  

&nbsp; Covers basic exploratory data analysis and simple summary statistics on genomic-style datasets.



\- `module\_2\_project.R`  

&nbsp; Script for the \*\*module 2 project\*\*.  

&nbsp; Implements statistical comparisons (e.g. hypothesis tests, confidence intervals) relevant to genomic measurements.



\- `module\_3\_project.R`  

&nbsp; Script for the \*\*module 3 project\*\*.  

&nbsp; Focuses on modelling approaches (e.g. simple regression or related models) applied to genomic or high-throughput data.



\- `module\_4\_project.R`  

&nbsp; Script for the \*\*module 4 project\*\*.  

&nbsp; Brings together previous topics into a small end-to-end analysis: data import, cleaning, statistical testing, and result visualisation.



> Note: each script is aligned with a specific course project. They are small, targeted analyses rather than a general-purpose statistics package.



---



\## How to run



1\. Open R or RStudio in this folder:



&nbsp;  ```r

&nbsp;  setwd("6 - Statistics for Genomic Data Science")

&nbsp;  ```

2\. Make sure you have R and the required CRAN packages installed. From R:



   ```r

&nbsp;  install.packages(c("tidyverse", "data.table"))

&nbsp;  # add other packages here if a specific script reports that one is missing

&nbsp;  ```



3\. Run a project script either by sourcing it:



&nbsp;  ```r

&nbsp;  source("module\_1\_project.R")

&nbsp;  ```



