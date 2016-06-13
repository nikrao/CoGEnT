
FOCUSS-CNDL Dictionary Learning Algorithms for Sparse Representations

    Matlab Code Documentation

    Contact:
        Joseph Murray
        jfmurray@ucsd.edu
        http://dsp.ucsd.edu/~jfmurray

        Kenneth Kreutz-Delgado
        kreutz@ece.ucsd.edu
        http://dsp.ucsd.edu/~kreutz

    Copyright 2005 Joseph F. Murray
    
---- Overview ----

The Matlab code in this directory implements the FOCUSS-CNDL learning
algorithm from the paper,

  @ARTICLE{Kreutz:2003,
    author =       {Kenneth Kreutz-Delgado and Joseph F. Murray and Bhaskar D. Rao and Kjersti Engan and Te-Won Lee and Terrence J. Sejnowski},
    title =        {Dictionary Learning Algorithms for Sparse Representation},
    journal =      {Neural Computation},
    year =         {2003},
    volume =       {15},
    number =       {2},
    pages =        {349-396},
    month =        {February}
  }

This code is still in a rough, experimental form.  Please contact us (see
above) if you have any questions or difficulties.

The software is distributed according to the license described in the
license.txt file included with this distribution.

---- Learning Image Dictionaries ----

Start by creating a data set of patches drawn from the images, using
createimagedata.m.  This loads all the images in a specified directory
draws small patches at random, and saves them to a .mat file.

The dictionaries are learned based on this .mat data by runtrials_image.m,
which calls the actual learning algorithm, found in trainrd_focusscndl.m
or trainrd_focusscndl2.m (both should give similar results).

Once a dictionary has been learned, results can be plotted with
plotresultsimage.m, which also calculates the entropy of the learned
coefficients using a method of similar to Lewicki:1999.

A stand-alone version of the FOCUSS algorithm, focuss.m, can be used to
find the coding of a new data set once the dictionary is learned.

---- Synthetic Dictionaries ----

The 20x30 synthetic dictionary experiment is run using
runtrials_focusscndl.m, which generates random dictionaries.

Results are plotted with plotresults.m, which compares the learned
dictionary to the true generating dictionary.
