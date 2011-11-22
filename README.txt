This directory contains the accompanying code for the paper, "Unsupervised Discovery of Rhyme Schemes", Sravana Reddy and Kevin Knight.

Our implementation is in Python using Numpy, and contains:

1. findschemes.py  The learning algorithm described in Section 3.

2. evaluate.py  Code to evaluate output of above, as well as various parsing functions.

Right now, our implementation of the algorithm for stanza dependencies is very slow. We would like to hold off on optimizing the code before releasing it. That said, it is just a simple modification of findschemes.py using forward-backward. (Unfortunately, even findschemes.py is reasonably slow for large datasets; we apologize.)

The file allschemes.pickle is a serialization of the complete list of rhyme schemes, with overall frequencies (used in estimating the naive baseline). This should be placed in the same directory as the code.

To learn rhyme schemes, use the command 

findschemes.py <goldfile> <initialization-type> <output>

<goldfile> is the gold standard data like the files in english_gold and french_gold. The code only reads the stanzas, and obviously, makes no use of the annotations.

<initialization-type> is a character that specifies whether to initialize theta uniformly (u), with the orthographic similarity measure (o), or using CELEX pronunciations and definition of rhyme (p). The last requires you to have CELEX on your machine.

<output> is simply the name of the desired output file. The program writes stanzas and annotations in a format similar to the gold standard.

To evaluate, do

evaluate.py <goldfile> <output>

For example, to learn rhyme schemes from Kipling's poety with uniform initialization,

findschemes.py kipling.pgold u kipling.out

or poetry from 1450-1550 with orthographic initialization:

findschemes.py 1415.pgold o 1415.out

To evaluate the above runs:

evaluate.py kipling.pgold kipling.out
evaluate.py 1415.pgold 1415.out

E-mail sravana@cs.uchicago.edu with any questions or bug-fixes.

