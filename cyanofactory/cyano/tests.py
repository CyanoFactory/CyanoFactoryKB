'''
Whole-cell knowledge base tests

Author: Jonathan Karr, jkarr@stanford.edu
Affiliation: Covert Lab, Department of Bioengineering, Stanford University
Last updated: 2012-07-17
'''

from django.db import models as dj_models
from django.test import TestCase

from public import helpers, models

'''
TODO
- import/excel works correctly
'''

class CyanoTest(TestCase):
	# Permission list
	fixtures = ['metadata.json']
	
	def setUp(self):
		# Create users and groups
		from management.commands.autocreateinitial import Command
		Command()
	
	def testTestWorks(self):
		self.assertFalse(True == False, "Not False")
