import sublime, sublime_plugin
from Bio.Seq import Seq
from Bio import SeqIO
import io
import re

class GenbankToFastaCommand(sublime_plugin.TextCommand):

	def run(self,edit):

		for region in self.view.sel():  
			seq_str = self.view.substr(region)
			seq_str = seq_str.strip()

			# Check that the selection begins as expected
			startmatch = re.match( r'^LOCUS', seq_str )
			# It turns out that SeqIO can handle Genbank format that
			# does not end in '//' so there is no need to check for this

			if not startmatch:
				print("Selected text must begin with 'LOCUS'")
			else:
				# Read from a string and write to a string
				seqout = io.StringIO()

				with io.StringIO(seq_str) as seqin:
					SeqIO.convert(seqin, 'genbank', seqout, 'fasta')
				seqin.close()

				# Write the fasta string to a new window at position 0			
				self.view.window().new_file().insert(edit, 0, seqout.getvalue())
