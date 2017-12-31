import logging
import os
from krocus.PubmlstGetter import PubmlstGetter

class KrocusDatabaseDownloader:
	def __init__(self,options):
		self.logger = logging.getLogger(__name__)
		self.list_species = options.list_species
		self.species = options.species
		self.output_directory = options.output_directory
		self.verbose = options.verbose
		
		if self.verbose:
			self.logger.setLevel(logging.DEBUG)
		else:
			self.logger.setLevel(logging.ERROR)
			
	def run(self):
		pubmlst  = PubmlstGetter(verbose = self.verbose)
		
		if self.list_species:
			pubmlst.print_available_species()
		elif self.species and self.output_directory:
			pubmlst.get_species_files(self.species, self.output_directory)
		else:
			self.logger.error("Please check the input parameters")
