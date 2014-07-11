import sys

biopython_location = view.settings().get('biopython_location')

if biopython_location:
	sys.path.append(biopython_location)
else:
	print("No BioPython location specified in BioPython Settings - User")