#! /bin/bash

grep SSQ best_fit.txt | awk '{print $3}' > SSQ
wc -l SSQ
xmgrace SSQ 

