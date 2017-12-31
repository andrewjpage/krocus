import re
class Gene:
	def __init__(self,name, kmers_with_coverage, kmers_without_coverage):
		self.name = name
		self.kmers_with_coverage = kmers_with_coverage
		self.kmers_without_coverage = kmers_without_coverage
		
	def __str__(self):
		s = self.allele_name()+'('+ str(self.allele_number())+')'
		if not self.is_full_coverage():
			s += '*'
		return s
		
	def is_full_coverage(self):
		if self.kmers_without_coverage  == 0:
			return True
		else:
			return False
			
	def allele_number(self):
		return int(re.split('\.|-|_',self.name)[-1])
				
	def allele_name(self):
		regex = r"(.+)[\.\-_]" + str(self.allele_number() )
		
		m = re.search(regex, self.name)
		 
		if m.group:
			return m.group(1)
		else:
			return ''
		