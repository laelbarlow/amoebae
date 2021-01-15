
## Additional background for new bioinformaticians


1. Like most bioinformatics software, AMOEBAE is run using the command line.
   So, to run AMOEBAE effectively you will definitely need to possess some
   basic linux/macOS command line skills. These might be good places to start
   (also, just look it up on YouTube):

    https://www.techspot.com/guides/835-linux-command-line-basics/

    https://datacarpentry.org/shell-genomics/

    https://github.com/LeCoupa/awesome-cheatsheets/blob/master/languages/bash.sh


2. With the exception of some very small analyses, AMOEBAE is best suited for
   running on High Performance Computing (HPC) Clusters. If you have access to
   such a cluster, your system administrators will provide documentation
   specific to the particular cluster. Make sure that you are familiar with the
   process of running software by submitting jobs appropriately to the job
   scheduler that the cluster uses to manage resources. Importantly, when you
   run AMOEBAE on a cluster, jobs will be submitted for you automatically as
   necessary to complete each step of your analysis (this works because of
   Snakemake, see below). This has been tested on clusters including [Compute
   Canada](https://www.computecanada.ca/). However, AMOEBAE may need to be
   customized to run on some clusters, depending on the job scheduler used. If
   you have trouble with this, please submit an issue through the GitHub [issue
   tracker for AMOEBAE](https://github.com/laelbarlow/amoebae/issues).


3. If you are a graduate student in bioinformatics or a related field and you
   have not yet learned to code in the widely used Python programming language
   and also learned to use the popular Conda package and environment manager,
   then you need to do that now. Among other things, this will allow you to
   customize AMOEBAE for your project to a greater degree than you could
   otherwise. These are good places to start:

    https://wiki.python.org/moin/BeginnersGuide

    https://www.anaconda.com/

    https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html

    https://docs.conda.io/projects/conda/en/latest/user-guide/cheatsheet.html    


4. AMOEBAE is a Snakemake workflow. Snakemake is a highly popular workflow and
   environment management system based on the Python programming language. To
   run a Snakemake workflow and to customize parameters (such as E-value
   thresholds) you do not need to know how to write Python code. However, at a
   minimum, you do need to know how to find code for specific steps (or
   "rules") in the workflow definition file (the Snakefile file) and understand
   their basic syntax. This will allow you to modify the linux shell commands
   contained in the rules. The good news is that this will be incompariably
   easier and more effective than entering commands manually one at a time to
   run your similarity search analysis. Here are some places to start
   developing your understanding of Snakemake:

    https://snakemake.readthedocs.io/en/stable/index.html
    
    https://hpc-carpentry.github.io/hpc-python/11-snakemake-intro/

    https://www.youtube.com/watch?v=hPrXcUUp70Y

    https://zenodo.org/record/4063463


