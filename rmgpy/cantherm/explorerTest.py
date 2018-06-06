#!/usr/bin/env python
# -*- coding: utf-8 -*-

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2018 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

import unittest
import numpy
import os
from nose.plugins.attrib import attr

import rmgpy
from rmgpy.cantherm import CanTherm, input
from input import jobList
import rmgpy.constants as constants
from rmgpy.cantherm.statmech import InputError
from rmgpy.cantherm.explorer import ExplorerJob
################################################################################

@attr('functional')
class testExplorerJob(unittest.TestCase):
    """
    Contains tests for ExplorerJob class execute method
    """
    def setUp(self):

        cantherm = CanTherm()
        
        self.jobList = cantherm.loadInputFile(os.path.join(os.path.dirname(os.path.abspath(__file__)),'data','methoxy_explore.py'))
        for job in self.jobList:
            if not isinstance(job,ExplorerJob):
                job.execute(outputFile=None, plot=None)
            else:
                thermoLibrary,kineticsLibrary,speciesList = cantherm.getLibraries()
                job.execute(outputFile=None, plot=None, speciesList=speciesList, thermoLibrary=thermoLibrary, kineticsLibrary=kineticsLibrary)
        
        self.thermoLibrary = thermoLibrary
        self.kineticsLibrary = kineticsLibrary
        self.explorerjob = self.jobList[-1]
        self.pdepjob = self.jobList[-2]
    
    def test_reactions(self):
        """
        test that the right number of reactions are in output network
        """
        self.assertEqual(len(self.explorerjob.network.pathReactions),4)
    
    def test_isomers(self):
        """
        test that the right number of isomers are in the output network
        """
        self.assertEqual(len(self.explorerjob.network.isomers),2)
    
    def test_job_rxns(self):
        """
        test that in this case all the reactions in the job
        ended up in the final network
        """
        for rxn in self.explorerjob.jobRxns:
            self.assertIn(rxn,self.explorerjob.network.pathReactions)


    
        
if __name__ == '__main__':
    unittest.main(testRunner=unittest.TextTestRunner(verbosity=2))
