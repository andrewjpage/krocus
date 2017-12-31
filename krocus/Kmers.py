'''Given a string of nucleotides and k, return all kmers'''

class Kmers:
	def __init__(self, sequence, k):
		self.sequence = sequence
		self.k = k

	def get_all_kmers(self):
		kmers = {}
		
		end = len(self.sequence) - self.k
		for i in range(0,end):
			seq = self.sequence[i:i+self.k]
			kmers[seq] = 0
		return kmers
		
	def get_all_kmers_array(self):
		kmers = []
		
		end = len(self.sequence) - self.k
		for i in range(0,end):
			kmers.append( self.sequence[i:i+self.k] )
		return kmers
