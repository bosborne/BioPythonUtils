import sublime, sublime_plugin
import io
import re
import os
import sys

class GenbankToFastaCommand(sublime_plugin.TextCommand):

	def run(self,edit):

		from Bio import SeqIO

		for region in self.view.sel():

			seq_str = self.view.substr(region)
			seq_str = seq_str.strip()

			if not seq_str:
				sublime.error_message("No selected text")
			else:
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
