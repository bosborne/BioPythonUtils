import sublime, sublime_plugin
import re
import io
import sys
import os

class TranslateCommand(sublime_plugin.TextCommand):
	# Valid chars
	valid_bases = ['A','T','G','C','U']

	def run(self,edit):

		from Bio import SeqIO
		from Bio.Seq import Seq
		from Bio.SeqRecord import SeqRecord
		from Bio.Alphabet import IUPAC

		for region in self.view.sel():  
			seq_str = self.view.substr(region)
		
			# Fasta header pattern
			patt = re.compile('^>\s*\S+')

			# If the selection looks like Fasta
			if patt.match(seq_str):
				# Read from a string and write to a string
				seqout = io.StringIO()

				with io.StringIO(seq_str) as seqin:
					for nt_seq_record in SeqIO.parse(seqin, "fasta"):

						if len( str(nt_seq_record.seq) ) < 3:
							print("Sequence is too short to translate: " + str(nt_seq_record.seq) )
							continue

						# translate() returns a string
						aa_seq = nt_seq_record.seq.translate()
						# Copy data to the protein SeqRecord from the starting nucleotide sequence
						aa_seq_record = SeqRecord(aa_seq, 
												  id=nt_seq_record.id, 
												  description=nt_seq_record.description)
						# Collect the records
						SeqIO.write(aa_seq_record, seqout, "fasta")

					seqin.close()

				# Write the fasta string to a new window at position 0			
				self.view.window().new_file().insert(edit, 0, seqout.getvalue())

			else:
				# If it's not fasta then remove non-alphabetic ...
				seq_str = re.sub(r'[^A-Za-z]+', r'', seq_str)
				# and check that it's at least one codon long ...
				if len(seq_str) > 2:
					nt_seq = Seq(seq_str, IUPAC.unambiguous_dna)

					invalid_chars = self.validate(str(nt_seq))
					# and that it's all valid chars
					if len(invalid_chars) > 0:
						sublime.error_message("Invalid characters in sequence: " + ' '.join(invalid_chars))
					else:
						aa_seq = nt_seq.translate()
						# Write the translated sequence to a new window
						self.view.window().new_file().insert(edit, 0, str(aa_seq))
				else:
					sublime.error_message("Selection is too short to translate: " + seq_str)

	# Checks that a sequence only contains values from an alphabet
	def validate(self, seq):
		seq_arr = list(seq.upper())
		invalid = list(set(seq_arr) - set(self.valid_bases))
		return invalid

# def plugin_loaded():

# 	settings = sublime.load_settings('BioPythonUtils.sublime-settings')
# 	biopython_location = settings.get('package_directory')

# 	if biopython_location:
# 		# This approach works if the specified path has a trailing slash or not
# 		if os.path.exists( os.path.join(os.path.sep, biopython_location, 'Bio') ):
# 			sys.path.append(biopython_location)
# 		else:
# 			sublime.error_message("'Bio' directory not found in directory '" + biopython_location + "'")
# 	else:
# 		sublime.error_message("Enter package directory in BioPythonUtils -> Settings - User")
