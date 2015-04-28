BioPythonUtils
==============

[BioPython](http://biopython.org) Utilities for [Sublime Text 3](http://www.sublimetext.com/3)

### Install [BioPython](http://biopython.org) for Python 3

* Do "easy_install -f http://biopython.org/DIST/ biopython"
* Or do "pip3 install biopython"
* Or see http://biopython.org/wiki/Download

### Install BioPythonUtils using [Package Control](https://sublime.wbond.net)

If you don't have Package Control see https://sublime.wbond.net.

### Configure with the [BioPython](http://biopython.org) location, email, or [BLAST](http://blast.ncbi.nlm.nih.gov/Blast.cgi) defaults

* Sublime Text comes with its own embedded Python 3 interpreter
* This interpreter needs to know where [BioPython](http://biopython.org) is installed so ...
* Enter the name of the directory containing BioPython:
* Preferences -> Packages Setting -> BioPythonUtils -> Settings - User  
* Add the "package directory", for example:
~~~~
{
    "package_directory": "/usr/local/lib/python3.4/site-packages",
    "email_for_eutils": "bio@bioteam.net",
    "remote_blast_app": "blastp",
    "remote_blast_db": "nr",
    "remote_blast_format": "Text"
}
~~~~

Add an email address to your "Settings - User" file ("email_for_eutils") if you want to download from [NCBI](http://www.ncbi.nlm.nih.gov) using [EUtils](http://www.ncbi.nlm.nih.gov/books/NBK25500).

### Commands

* "Translate"
* "Download Sequence"
* "Download Taxon"
* "Remote BLAST"
* "Genbank To Fasta"

#### "Translate"

Translates the selected text, which can be 1 or more entries in Fasta format or 1 or more entries of plain text. For example:
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

atgctgtcaatcacgattctgttggtgctcatagcagagggctcctctcagaattacacg
gggagtcctgtgatatgcctgggacatcatgctgtatccaatgggacaatggtgaaaacg
~~~~

#### "Download Sequence" 

Downloads sequence from [NCBI](http://www.ncbi.nlm.nih.gov) using the selected ids. Ids can be delimited by commas, returns, or space. For example:
~~~~
KC781785 2
~~~~
or:
~~~~
2
KC781786
~~~~
or:
~~~~
284218, 203807
~~~~

#### "Download Taxon" 

Downloads nucleotide sequence from [NCBI](http://www.ncbi.nlm.nih.gov) using the selected [NCBI Taxonomy](http://www.ncbi.nlm.nih.gov/taxonomy) ids. Ids can be delimited by commas, returns, or space.

#### "Remote BLAST"

Sends the selected Fasta format or "plain" sequence(s) to the [BLAST server at NCBI](http://blast.ncbi.nlm.nih.gov/Blast.cgi) and retrieves the results. Set the application, database, and result format using the Command Palette. You can also set default values for these in your "Settings - User" file ("remote_blast_app", "remote_blast_db", "remote_blast_format").

#### "Genbank To Fasta"

Converts the selection, 1 or more Genbank entries, to Fasta format.

### Issues

The interaction between the plugin and various services at NCBI are  synchronized, so Sublime Text is essentially unusable while the queries ("Download ...", Remote BLAST) are running. Not suitable for large-scale work, need to fix this.
