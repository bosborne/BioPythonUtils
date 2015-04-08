import io
import re
import sublime
import sublime_plugin
import threading
import time
import Bio
from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIWWW

# Used by "Remote BLAST"
blast_app = None
blast_db = None
blast_format = None


# "Download Sequence"
class DownloadSequenceCommand(sublime_plugin.TextCommand):

    # The underscore is used to distinguish the method name from 
    # the variable name, this is necessary to prevent recursion
    _sequences = ''

    def sequences(self):
        return self._sequences

    def add_to_sequences(self, seq):
        self._sequences = self._sequences + seq

    def run(self, edit):

        email_for_eutils = sublime.load_settings(
            'BioPythonUtils.sublime-settings').get('email_for_eutils')

        if not email_for_eutils:
            sublime.error_message(
                "Enter email address for EUtils in BioPythonUtils -> Settings - User")
            return

        for region in self.view.sel():

            id_str = self.view.substr(region)
            id_str = id_str.strip()

            if not id_str:
                sublime.error_message("No identifiers in selection")
                continue

            ids = re.split('[\n\s,]+', id_str)
            threads = []
            results = ''

            # Start jobs
            for id in ids:
                thread = EutilsCall(id=id, email=email_for_eutils)
                # AY818147time.sleep(1)
                thread.start()
                threads.append(thread)

            results = handle_threads(threads)
            self.view.window().new_file().insert(edit, 0, results )


def handle_threads(threads, results=''):
    next_threads = []

    for thread in threads:
        if thread.is_alive():
            next_threads.append(thread)
            # Required, or else this loop runs so many times that
            # the code exits, too much "recursion"
            time.sleep(.01)
        elif thread.result != False:
            print("Downloaded " + thread.id)
            results = results + str(thread.result)

    if len(next_threads) > 0:
        # Since this is recursive the "return" is required!
        return handle_threads(next_threads, results)
        # sublime.set_timeout(lambda: self.handle_threads(next_threads), 100)
    elif len(next_threads) == 0:
        return results


class EutilsCall(threading.Thread):
    # Subclass so that the results can be accessed through the instance variable
    def __init__(self, id=None, result=None, email=None):
        self.id = id
        self.email = email
        self.result = result
        threading.Thread.__init__(self)

    def run(self):
        try:
            handle = Entrez.efetch(db="nucleotide",
                                   id=self.id,
                                   email=self.email,
                                   rettype="gb",
                                   retmode="text")
            self.result = handle.read()
            return

        except (IOError) as exception:
            print(str(exception))
            sublime.error_message(
                "Error retrieving sequence using id '" + self.id + "'")

        self.result = False


# "Download Taxon"
class DownloadTaxonCommand(sublime_plugin.TextCommand):

    def run(self, edit):

        email_for_eutils = sublime.load_settings(
            'BioPythonUtils.sublime-settings').get('email_for_eutils')

        if email_for_eutils:
            Entrez.email = email_for_eutils
        else:
            sublime.error_message(
                "Enter email address for EUtils in BioPythonUtils -> Settings - User")
            return

        for region in self.view.sel():

            id_str = self.view.substr(region)
            id_str = id_str.strip()

            if not id_str:
                sublime.error_message("No identifiers in selection")
                continue

            taxids = re.split('[\n\s,]+', id_str)
            seq_txt = ''

            for taxid in taxids:
                # Check that they're numeric
                nummatch = re.match(r'^\d+$', taxid)

                if not nummatch:
                    sublime.error_message(
                        "String '" + taxid + "' is not an NCBI taxon id")
                    return

                try:
                    links = Entrez.read(
                        Entrez.elink(dbfrom="taxonomy",
                                     db="nucleotide",
                                     id=taxid))
                except (IOError) as exception:
                    print(str(exception))
                    sublime.error_message(
                        "Error retrieving sequence ids using id '" + taxid + "'")

                if len(links[0]["LinkSetDb"]) == 0:
                    sublime.error_message(
                        "No sequences retrieved with id " + taxid)
                    return

                for link in links[0]["LinkSetDb"][0]["Link"]:
                    try:
                        handle = Entrez.efetch(db="nucleotide",
                                               id=link['Id'],
                                               rettype="gb",
                                               retmode="text")
                    except (IOError) as exception:
                        print(str(exception))
                        sublime.error_message("Error retrieving sequence using id '" +
                                              link['Id'] + "':" + str(exception))

                    seq_txt = seq_txt + handle.read()

            # Write the fasta string to a new window at position 0
            self.view.window().new_file().insert(edit, 0, seq_txt)


# "Translate"
class TranslateCommand(sublime_plugin.TextCommand):

    def run(self, edit):

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

                        if len(str(nt_seq_record.seq)) < 3:
                            sublime.error_message(
                                "Sequence is too short to translate: " +
                                str(nt_seq_record.seq))
                            return

                        try:
                            # translate() returns a string
                            aa_seq = nt_seq_record.seq.translate()
                        except (BaseException) as exception:
                            print(str(exception))
                            sublime.error_message(str(exception))

                        # Copy data to the protein SeqRecord from the starting
                        # nucleotide sequence
                        aa_seq_record = SeqRecord(aa_seq,
                                                  id=nt_seq_record.id,
                                                  description=nt_seq_record.description)
                        # Collect the records
                        SeqIO.write(aa_seq_record, seqout, "fasta")

                    seqin.close()

                # Write the fasta string to a new window at position 0
                self.view.window().new_file().insert(
                    edit, 0, seqout.getvalue())

            # If selection is not Fasta
            else:
                seqout = []
                seq_num = 1
                # Could be more than one sequence, "MULTILINE" required
                for nt_str in re.split('^\s*\n', seq_str, 0, re.MULTILINE):

                    # Remove non-alphabetic ...
                    nt_str = re.sub(r'[^A-Za-z]+', r'', nt_str)
                    # and check that it's at least one codon long
                    if len(nt_str) < 3:
                        sublime.error_message(
                            "Selection is too short to translate: " + nt_str)
                        return

                    nt_seq = Seq(nt_str, IUPAC.unambiguous_dna)
                    # Check that it's all valid nucleotide
                    invalid_chars = validate_nt(str(nt_seq))

                    if len(invalid_chars) == 0:
                        try:
                            aa_seq = nt_seq.translate()
                        except (BaseException) as exception:
                            print(str(exception))
                            sublime.error_message(str(exception))

                        seqout.append(str(aa_seq))
                        seq_num += 1
                    else:
                        sublime.error_message("Invalid characters in sequence "
                                              + str(seq_num) + ": " + ''.join(invalid_chars))
                        return

                # Separate the translations with an empty line
                self.view.window().new_file().insert(
                    edit, 0, "\n\n".join(seqout))


# "Genbank To Fasta"
class GenbankToFastaCommand(sublime_plugin.TextCommand):

    def run(self, edit):

        for region in self.view.sel():

            seq_str = self.view.substr(region)
            seq_str = seq_str.strip()

            if not seq_str:
                sublime.error_message("No selected text")
                return

            # Check that the selection begins as expected
            startmatch = re.match(r'^LOCUS', seq_str)
            # It turns out that SeqIO can handle Genbank format that
            # does not end in '//' so there is no need to check for this

            if startmatch:
                # Read from a string and write to a string
                seqout = io.StringIO()

                with io.StringIO(seq_str) as seqin:
                    SeqIO.convert(seqin, 'genbank', seqout, 'fasta')
                seqin.close()

                # Write the fasta string to a new window at position 0
                self.view.window().new_file().insert(
                    edit, 0, seqout.getvalue())
            else:
                sublime.error_message(
                    "Selected text does not look like Genbank: no 'LOCUS'")
                return


# "Remote Blast"
class RemoteBlastCommand(sublime_plugin.TextCommand):

    blast_apps = ['blastp', 'blastn', 'blastx', 'tblastn', 'tblastx']
    blast_dbs = ['nr', 'refseq', 'swissprot', 'pat', 'month', 'pdb', 'env_nr']
    blast_formats = ['HTML', 'Text', 'ASN.1', 'XML']

    def run(self, edit):

        global blast_app, blast_db, blast_format

        if blast_app is None:
            blast_app = sublime.load_settings(
                'BioPythonUtils.sublime-settings').get('remote_blast_app')

        if blast_db is None:
            blast_db = sublime.load_settings(
                'BioPythonUtils.sublime-settings').get('remote_blast_db')

        if blast_format is None:
            blast_format = sublime.load_settings(
                'BioPythonUtils.sublime-settings').get('remote_blast_format')

        if not blast_db:
            sublime.error_message("No BLAST database specified")
            return

        if not blast_app:
            sublime.error_message("No BLAST application specified")
            return

        if not blast_format:
            sublime.error_message("No BLAST format specified")
            return

        for region in self.view.sel():
            seq_str = self.view.substr(region)

            # Fasta header pattern
            patt = re.compile('^>\s*\S+')

            # If the selection looks like Fasta
            if patt.match(seq_str):
                # Read from a string and write to a string
                seqout = io.StringIO()

                with io.StringIO(seq_str) as seqin:
                    for seq_record in SeqIO.parse(seqin, "fasta"):
                        result = NCBIWWW.qblast(blast_app, blast_db, seq_record.format('fasta'),
                                                format_type=blast_format)
                        # Write the result to a new window at position 0
                        self.view.window().new_file().insert(
                            edit, 0, result.read())

                    seqin.close()
            else:
                # Assume it's plain sequence and use incrementing number as an id
                seq_id = 1
                for seq_str in re.split('^\s*\n', seq_str, 0, re.MULTILINE):
                    seq_str = re.sub("[^a-zA-Z]","", seq_str)
                    seq_record = SeqRecord(Seq(seq_str), id=str(seq_id))
                    result = NCBIWWW.qblast(blast_app, blast_db, seq_record.format('fasta'),
                                            format_type=blast_format)
                    seq_id += 1
                    # Write the result to a new window at position 0
                    self.view.window().new_file().insert(
                        edit, 0, result.read())


class SelectBlastDatabase(sublime_plugin.WindowCommand):

    def run(self):
        cmd = RemoteBlastCommand
        sublime.active_window().show_quick_panel(cmd.blast_dbs, setDatabase)


class SelectBlastApplication(sublime_plugin.WindowCommand):

    def run(self):
        cmd = RemoteBlastCommand
        sublime.active_window().show_quick_panel(
            cmd.blast_apps, setApplication)


class SelectBlastFormat(sublime_plugin.WindowCommand):

    def run(self):
        cmd = RemoteBlastCommand
        sublime.active_window().show_quick_panel(cmd.blast_formats, setFormat)


def setFormat(index):
    global blast_format
    cmd = RemoteBlastCommand

    if index > -1:
        blast_format = cmd.blast_formats[index]


def setDatabase(index):
    global blast_db
    cmd = RemoteBlastCommand

    if index > -1:
        blast_db = cmd.blast_dbs[index]


def setApplication(index):
    global blast_app
    cmd = RemoteBlastCommand

    if index > -1:
        blast_app = cmd.blast_apps[index]


def validate_nt(seq):
    # {'G', 'T', 'U', 'C', 'A'}
    valid_bases = set(IUPAC.unambiguous_dna.letters +
                      IUPAC.unambiguous_rna.letters)

    seq_arr = list(seq.upper())
    invalid = list(set(seq_arr) - valid_bases)
    return invalid


def validate_aa(seq):

    seq_arr = list(seq.upper())
    invalid = list(set(seq_arr) - set(IUPAC.protein.letters))
    return invalid
