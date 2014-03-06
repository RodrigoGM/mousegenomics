mouseomics
==========

This repository contains code to help with the setup of a Genomics cluster without
 root privilages.  The information detailed here, as the name implies, is oriented
 towards Mouse Genomics.  Thus, most scripts will be downloading mouse data from 
 Ensembl, UCSC or NCBI.

Software set up should be generic, though check to make sure you are downloading
  appropriate reference files for you own species of interest (not everyone is 
  interested in humans).  Pay special attention to the $PATH where things are
  located and that you're not deleting something important.

Analysis pipelines will most likely need to be adapted to your own style and 
  data structure. Though... take the ideas and best of luck with your work !