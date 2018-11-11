# FCOP

This package is developed with C
OS: Linux 2.6.18-194.26.1.el5 
Compiler: gcc version 4.4.7 20120313 (Red Hat 4.4.7-3) (GCC)

## PACKAGE CONTENT ##
Source Code:   
               FCOP.cpp                    source code
               fcop                        executable file

Example Data:  
               mapping_GM12878             the mapping matrix among enhancers and TF binding peaks
               mapping_GM12878_random      a randomized mapping matrix (simplely permutate each column of mapping_GM12878)
               ANNoutput_GM12878           probability of each enhancer (this is optional)
               GM12878_TFlist              list of the TFs

Sample Output: 
               Cuts_GM12878                Pattern Sepcific Cutoff, output of Step1
               output                      FCOPs, output of Step2


## GENERAL INSTRUCTIONS ##
For running instructions, just run the program without any inputs
./fcop
   
Usage: ./fcop [other options]
==================== Inputes ==========================
-f     -ffilename      File for input data
-X     -X5000         No. of rows of input data (genes/enhancers)
-Y     -Y62            No. of columns of input data (TFs)

-C     -CCuts_GM12878  Pattern specific cutoff file
-s     -s5            Pattern length
-m     -m80            Least support
-p     -p0.8          Minimal frequentness probability
-------------------- Optional Inputes ------------------
-t     -tTFlist    File for the TFs
-h     -hEnhancerProbabilities  Files for the enhancer confidence levels
==================== Examples ==========================
Step1: Generate Pattern Specific Cutoff at pvalue = 0.001
 -fmapping_GM12878_random -X5000 -Y62 -s5 -p0.001  -hEnhancerP  > Cuts_GM12878

Step2: Search with Frequentness Probability on Minimum Support and Pattern Specific Cutoff
 -fmapping_GM12878 -X5000 -Y62 -CCuts_GM12878 -m50 -p0.8 
 -fmapping_GM12878 -X5000 -Y62 -CCuts_GM12878 -m50 -p0.8 -tTFfile -hEhancerP
========================================================

## STEP BY STEP INSTRUCTIONS ##
Run the following 3 commonds and get it working!
Complie source code and get executable program,
./g++ -o fcop  FCOP.cpp

Step1. Generate pattern specific cutoff for all patterns with length up to 5
./fcop -fmapping_GM12878_rand -X5000 -Y62 -s5 -p0.01 -hANNoutput_GM12878  > Cuts_GM12878

Step2. Generate all FCOPs with frequentness probability larger than 0.8 at pattern specific cutoff or least support (whichever is larger)
./fcop -fmapping_GM12878 -X5000 -Y62 -CCuts_GM12878  -m80 -p0.8   > output
or
./fcop -fmapping_GM12878 -X5000 -Y62 -CCuts_GM12878  -m80 -p0.8 -tGM12878_TFlist >output

## OUTPUT FORMAT ##
columns of Cuts_GM12877
1.    pattern specific cutoff 
2...  index of the TFs in the pattern (start from 1)

columns of output
1.    number of support of the pattern 
2.    frequentness probability
3...  index of the TFs in the FCOP (start from 1)


## Notice ##
Theoretically, we can calculated pattern specific cutoff for all possible patterns. However, it is not economic to do that. 
For real data, longer patterns have much fewer support (in our analysis nearly all patterns longer than 5 have pattern specific cutof close to 1). 
To speed up the process, when generate the pattern sepcific cutoff (Step1) we implemented a maximum pattern length. In the example in 'step by step instructions', we only calculate pattern specific cutoff for patterns with length up to 5 (-s5).
Then, we implemented a 'least support' (-m80) when generating all FCOPs (Step2). The frequentness probability of a pattern is calculated based on the pattern specific cutoff or the least support, whichever is larger. For example, pattern candidate [TF1, TF2, TF3] has a pattern specific cutoff of 150, then its frequentness probability will be calculated based on 150, while pattern candidate [TF3 TF5 TF9] has a pattern specific cutoff of 50, then its frequentness probability will be calculated based on 80.
By doing so the non-frequent patterns are eliminated much more efficiently.





