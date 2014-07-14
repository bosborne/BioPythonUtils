import sublime, sublime_plugin
import io
import re

#
# "Download Sequence"
#
class DownloadSequenceCommand(sublime_plugin.TextCommand):

	def run(self,edit):

		from Bio import Entrez

		settings = sublime.load_settings('BioPythonUtils.sublime-settings')
		email_for_eutils = settings.get('email_for_eutils')

		if not email_for_eutils:
			sublime.error_message("Enter email address for EUtils in BioPythonUtils -> Settings - User")
			return
		else:
			Entrez.email = email_for_eutils

		for region in self.view.sel():  

			id_str = self.view.substr(region)
			id_str = id_str.strip()

			if not id_str:
				sublime.error_message("No identifiers in selection")
				continue

			ids = re.split( '[\n\s,]+', id_str )
			seq_txt = ''

			for id in ids:
				try:
					handle = Entrez.efetch(db="nucleotide", id=id, rettype="gb", retmode="text")
				except (IOError) as exception:
					print(str(exception))
					sublime.error_message("Error retrieving sequence using id '" + id + "'")

				seq_txt = seq_txt + handle.read()

			# Write the fasta string to a new window at position 0			
			self.view.window().new_file().insert(edit, 0, seq_txt)

#
# "Download Taxon"
#
class DownloadTaxonCommand(sublime_plugin.TextCommand):

	def run(self,edit):

		from Bio import Entrez

		settings = sublime.load_settings('BioPythonUtils.sublime-settings')
		email_for_eutils = settings.get('email_for_eutils')

		if not email_for_eutils:
			sublime.error_message("Enter email address for EUtils in BioPythonUtils -> Settings - User")
			return
		else:
			Entrez.email = email_for_eutils

		for region in self.view.sel():  

			id_str = self.view.substr(region)
			id_str = id_str.strip()

			if not id_str:
				sublime.error_message("No identifiers in selection")
				continue

			taxids = re.split( '[\n\s,]+', id_str )
			seq_txt = ''

			for taxid in taxids:
				try:
					links = Entrez.read(Entrez.elink(dbfrom="taxonomy", db="nucleotide", id=taxid))
				except (IOError) as exception:
					print(str(exception))
					sublime.error_message("Error retrieving sequence ids using id '" + taxid + "'")

				for link in links[0]["LinkSetDb"][0]["Link"]:
					try:
						handle = Entrez.efetch(db="nucleotide", id=link['Id'], rettype="gb", retmode="text")
					except (IOError) as exception:
						print(str(exception))
						sublime.error_message("Error retrieving sequence using id '" + link['Id'] + "'")

					seq_txt = seq_txt + handle.read()

			# Write the fasta string to a new window at position 0			
			self.view.window().new_file().insert(edit, 0, seq_txt)

#
# "Translate"
#
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
						return
					else:
						aa_seq = nt_seq.translate()
						# Write the translated sequence to a new window
						self.view.window().new_file().insert(edit, 0, str(aa_seq))
				else:
					sublime.error_message("Selection is too short to translate: " + seq_str)
					return

	# Checks that a sequence only contains values from an alphabet
	def validate(self, seq):
		seq_arr = list(seq.upper())
		invalid = list(set(seq_arr) - set(self.valid_bases))
		return invalid

#
# "Genbank To Fasta"
#
class GenbankToFastaCommand(sublime_plugin.TextCommand):

	def run(self,edit):

		from Bio import SeqIO

		for region in self.view.sel():

			seq_str = self.view.substr(region)
			seq_str = seq_str.strip()

			if seq_str:
				# Check that the selection begins as expected
				startmatch = re.match( r'^LOCUS', seq_str )
				# It turns out that SeqIO can handle Genbank format that
				# does not end in '//' so there is no need to check for this

				if not startmatch:
					sublime.error_message("Selected text does not look like Genbank: no 'LOCUS'")
				else:
					# Read from a string and write to a string
					seqout = io.StringIO()

					with io.StringIO(seq_str) as seqin:
						SeqIO.convert(seqin, 'genbank', seqout, 'fasta')
					seqin.close()

					# Write the fasta string to a new window at position 0			
					self.view.window().new_file().insert(edit, 0, seqout.getvalue())
			else:
				sublime.error_message("No selected text")
