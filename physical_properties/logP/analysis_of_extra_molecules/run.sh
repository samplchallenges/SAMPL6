#!/bin/bash

# Activate environment
source activate sampl6_logP

# Update requirements
conda list --export > requirements.txt

# Run first step of analysis and create collection file
python logP_analysis.py
python logP_analysis2.py

# Compile LaTeX statistic table twice for better rendering
pdflatex ./analysis_outputs_withrefs/StatisticsTables/statisticsLaTex/statistics.tex
pdflatex ./analysis_outputs_withrefs/StatisticsTables/statisticsLaTex/statistics.tex
rm statistics.log 
rm statistics.aux 
rm texput.log
mv statistics.pdf ./analysis_outputs_withrefs/StatisticsTables/  

