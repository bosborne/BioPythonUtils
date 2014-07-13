import sys
import sublime, sublime_plugin
import os

settings = sublime.load_settings('BioPythonUtils.sublime-settings')
biopython_location = settings.get('package_directory')

if biopython_location:
	# This approach works if the specified path has a trailing slash or not
	if os.path.exists( os.path.join(os.path.sep, biopython_location, 'Bio') ):
		sys.path.append(biopython_location)
	else:
		sublime.error_message("'Bio' directory not found in directory '" + biopython_location + "'")
else:
	sublime.error_message("Enter BioPython location in BioPythonUtils -> Settings - User")
