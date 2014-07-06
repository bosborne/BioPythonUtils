BioPythonUtils
==============

BioPython utilities for Sublime Text 3

### Install BioPython for Python 3
- Do "sudo easy_install -f http://biopython.org/DIST/ biopython"
- Or see http://biopython.org/wiki/Download

### Configure Plugin.py with the BioPython location
- Sublime Text comes with its own embedded Python 3 interpreter
- This interpreter needs to know where BioPython is installed so ...
- Enter the name of directory containing BioPython into SublimeBio/Plugin.py:
- Preferences -> Packages Settings -> SublimeBio -> BioPython Location - User  
- For example:

sys.path.append('/usr/local/lib/python3.4/site-packages')

### Commands

- "Translate" translates the selected text, which can be Fasta or plain text. For example:
~~~~
     >2
     atgctatcaatcgcgattctgcttctgctaatagcagagggctcctctcaaaattacaca
     ggaaatcctgtgatatgcctggggcaccatgctgtgtccaatgggacaatggtgaaaacc
     >1
     atgctatcaatcacgattctgttcctgctcatagcagagggctcctctcagaattacaca
     gggaatcctgtgatatgcctgggacatcatgctgtatccaatgggacaatggtgaaaacc
~~~~
or:
~~~~
    atgctatcaatcacgattctgttcctgctcatagcagagggctcctctcagaattacaca
    gggaatcctgtgatatgcctgggacatcatgctgtatccaatgggacaatggtgaaaacc
~~~~
- "Genbank To Fasta" converts the selection (GenBank) to Fasta format
