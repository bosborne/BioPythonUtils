BioPythonUtils
==============

[BioPython](http://biopython.org) Utilities for [Sublime Text 3](http://www.sublimetext.com/3)

### Install BioPythonUtils using [Package Control](https://sublime.wbond.net)

If you don't have Package Control see https://sublime.wbond.net.

[BioPython](http://biopython.org) 1.68 is bundled with this package, minus the Tests and
documentation, you do not need to install it.

### Configure with your email address, and [BLAST](http://blast.ncbi.nlm.nih.gov/Blast.cgi) defaults

~~~~
{
    "email_for_eutils": "bio@bioteam.net",
    "entrez_retmax": 1000,
    "remote_blast_app": "blastp",
    "remote_blast_db": "nr",
    "remote_blast_format": "Text"
}
~~~~

The email address is required if you want to download from [NCBI](http://www.ncbi.nlm.nih.gov)
using [EUtils](http://www.ncbi.nlm.nih.gov/books/NBK25500).

### Commands

First select the relevant text, then use a command in the Tools -> BioPythonUtils menu.


* "Translate"
* "Get Sequences by Search"
* "Get Sequences by Id"
* "Get Sequences by Taxon"
* "Remote BLAST"
* "Genbank To Fasta"

#### "Translate"

Translates the selected text, which can be 1 or more entries in Fasta format or 1 or more
entries of plain text. For example:
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
Translation starts at the first codon and continues to the last, regardless of stop codons.

#### "Get Sequences by Search"

Queries [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nucleotide/) using [Entrez](https://www.ncbi.nlm.nih.gov/books/NBK184582/)
search service and the selected text and downloads the sequence results in GenBank format.
For example:
~~~~
human[ORGN] AND AKT1[GENE]
~~~~

The default maximum number of records downloaded is 20, if you want more than 20 set
the value in your "Settings - User" file ("entrez_retmax").

#### "Get Sequences by Id"

Downloads sequence from [NCBI](http://www.ncbi.nlm.nih.gov) using the selected ids. Ids can be delimited by commas,
returns, or space. For example:
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

#### "Get Sequences by Taxon"

Downloads a taxon as GenBank format entries from [NCBI](http://www.ncbi.nlm.nih.gov) using the selected
[NCBI Taxonomy](http://www.ncbi.nlm.nih.gov/taxonomy) ids. Ids can be delimited by commas, returns, or space.

#### "Remote BLAST"

Sends the selected Fasta format or "plain" sequence(s) to the [BLAST server at NCBI](http://blast.ncbi.nlm.nih.gov/Blast.cgi) and retrieves the results. You can set the application, database, and result format using the Command Palette. You can also set some default values in your "Settings - User" file ("remote_blast_app", "remote_blast_format"). Note that the available databases changes depending on the BLAST application.

#### "Genbank To Fasta"

Creates Fasta format records from the selection, 1 or more GenBank format records.

### Issues

The interaction between the plugin and various services at NCBI are  synchronized, so Sublime Text is
unusable while the queries ("Get Sequences ...", "Remote BLAST") are running. Not suitable for large-scale work.
