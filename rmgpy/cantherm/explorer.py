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

import os
import numpy as np
import logging
import shutil

from rmgpy.rmg.model import CoreEdgeReactionModel
from rmgpy.rmg.pdep import PDepNetwork
from rmgpy.molecule.molecule import Molecule
from rmgpy import settings
from rmgpy.data.rmg import getDB
from rmgpy.species import Species

class ExplorerJob(object):
    def __init__(self, source, pdepjob, explore_tol, energy_tol=np.inf, flux_tol=0.0):
        self.source = source
        self.explore_tol = explore_tol
        self.energy_tol = energy_tol
        self.flux_tol = flux_tol
        
        self.pdepjob = pdepjob
        
        if not hasattr(self.pdepjob,'outputFile'):
            self.pdepjob.outputFile = None

    def copy(self):
        """
        Return a copy of the pressure dependence job.
        """
        return ExplorerJob(
            source=self.source,
            pdepjob=self.pdepjob,
            explore_tol=self.explore_tol,
            energy_tol=self.energy_tol,
            flux_tol=self.flux_tol
        )
        
    def execute(self, outputFile, plot, format='pdf', print_summary=True, speciesList=None, thermoLibrary=None, kineticsLibrary=None):
        
        logging.info('Exploring network...')
        
        reactionModel = CoreEdgeReactionModel()
        
        reactionModel.pressureDependence = self.pdepjob
        
        reactionModel.pressureDependence.rmgmode = True
        
        if outputFile:
            reactionModel.pressureDependence.outputFile = os.path.dirname(outputFile)
        
        kineticsDatabase = getDB('kinetics')
        thermoDatabase = getDB('thermo')
        
        thermoDatabase.libraries['thermojobs'] = thermoLibrary
        thermoDatabase.libraryOrder.insert(0,'thermojobs')
        
        kineticsDatabase.libraries['kineticsjobs'] = kineticsLibrary
        kineticsDatabase.libraryOrder.insert(0,('kineticsjobs','Reaction Library'))
        
        reactionModel.addSeedMechanismToCore('kineticsjobs')
        
        jobRxns = [rxn for rxn in reactionModel.core.reactions]
        
        self.jobRxns = jobRxns
        
        for lib in kineticsDatabase.libraryOrder:
            if lib[0] != 'kineticsjobs':
                reactionModel.addReactionLibraryToEdge(lib[0])

        
        
        
        