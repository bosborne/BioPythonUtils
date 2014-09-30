BioPythonUtils
==============

BioPython utilities for Sublime Text 3

### Install BioPython for Python 3

* Do "easy_install -f http://biopython.org/DIST/ biopython"
* Or do "pip3 install biopython"
* Or see http://biopython.org/wiki/Download

### Install BioPythonUtils using Package Control

If you don't have Package Control see https://sublime.wbond.net.

### Configure with the BioPython location, email, or BLAST defaults

* Sublime Text comes with its own embedded Python 3 interpreter
* This interpreter needs to know where BioPython is installed so ...
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

### Commands

* "Translate"
* "Download Sequence"
* "Download Taxon"
* "Genbank To Fasta"
* "Remote BLAST"

#### "Translate"

Translates the selected text, which can be one or more entries in Fasta format or 1 or more entries of plain text. For example:
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

Downloads sequence from NCBI using the selected ids. Ids can be delimited by commas, returns, or space. For example:
~~~~
KC781785 2
~~~~
or:
~~~~
2
KC781786
~~~~

Add an email address to your "Settings - User" file ("email_for_eutils") if you want to download from NCBI using EUtils.

#### "Genbank To Fasta"

Converts the selected Genbank entries to Fasta format.

#### "Download Taxon" 

Downloads nucleotide sequence from NCBI using the selected NCBI Taxonomy ids. Ids can be delimited by commas, returns, or space.

For example:
~~~~
284218, 203807
~~~~

Add an email address to your "Settings - User" file ("email_for_eutils") if you want to download from NCBI using EUtils.

#### "Remote BLAST"

Sends the selected Fasta format or "plain" sequence(s) to the BLAST server at NCBI and retrieves the results. Set the application, database, and result format using the Command Palette. You can also set default values for these in your "Settings - User" file ("remote_blast_app", "remote_blast_db", "remote_blast_format").

### Issues

The interaction between the plugin and various services at NCBI are  synchronized, so Sublime Text is essentially unusable while the queries ("Download ...", Remote BLAST) are running. Not suitable for large-scale work, need to fix this.
