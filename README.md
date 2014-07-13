BioPythonUtils
==============

BioPython utilities for Sublime Text 3

### Install BioPython for Python 3
- Do "sudo easy_install -f http://biopython.org/DIST/ biopython"
- Or see http://biopython.org/wiki/Download

### Install BioPythonUtils using Package Control
- If you don't have Package control see https://sublime.wbond.net

### Configure with the BioPython location and email address
- Sublime Text comes with its own embedded Python 3 interpreter
- This interpreter needs to know where BioPython is installed so ...
- Enter the name of directory containing BioPython:
- Preferences -> Packages Settings -> BioPythonUtils -> Settings - User  
- For example:
~~~~
{
    "package_directory": "/usr/local/lib/python3.4/site-packages",
    "email_for_eutils": "bio@bioteam.net"
}
~~~~
- The email address is required if you want to download from NCBI using EUtils

### Commands

- "Translate" translates the selected text, which can be one or more entries in Fasta format or 1 entry of plain text. For example:
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
- "Genbank To Fasta" converts the selected entries (GenBank) to Fasta format

- "Download Nucleotide" downloads sequence from NCBI using the selected ids. Ids can be delimited by commas, returns, or space. For example:
~~~~
KC781785 2
~~~~
or:
~~~~
2
KC781786
~~~~
