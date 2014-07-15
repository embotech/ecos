readme.txt for pdco4

The software for PDCO is provided by SOL, Stanford University
under the terms of the OSI Common Public License (CPL)
   http://www.opensource.org/licenses/cpl1.0.php
or the BSD License
   http://www.opensource.org/licenses/bsd-license.php

17 Oct 2002: First set of files available for download from SOL.
24 Nov 2013: Files reorganized as pdco4/*.

Development history for pdco.m is included as comments within code/pdco.m.

Please send comments to
             Michael Saunders, SOL, Stanford University
             saunders@stanford.edu  650-723-1875
-----------------------------------------------------------------------------

The PDCO files are organized into the following folders:

  code  Solver files: pdco.m, pdcoSet.m, lsmr.m, minres.m
        Test files:   pdcotest*.m, LPnetlib*.m
  data  ENTROPY.*     used by pdcotestENTROPY.m
        lp*.mat       used by LPnetlib.m
  doc   notes07-PDinterior.pdf: class notes about PDCO and interior methods
        project-Method22.pdf:   report by Matt Zahr implementing Method=22

To test PDCO:
   cd matlab/pdco4 (or wherever you put pdco4)
   addpath code
   addpath data
Type any of
   help pdcotestBPDN
   help pdcotestENTROPY
   help pdcotestLP
   help pdcotestLS
   help pdcotestQP
and copy the first line of the "help" documentation into Matlab.
