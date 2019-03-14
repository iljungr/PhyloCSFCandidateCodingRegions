#!/usr/bin/env python
# Copyright 2019 Irwin Jungreis
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
CommonUtils.py
Common utilities.
"""
from __future__ import division, print_function
import sys, os, shutil, hashlib
if sys.version_info.major == 2 and sys.version_info.minor == 7 and sys.version_info.micro < 10 :
    # Overriding update and __get_item__ of OrderedDict is broken in early versions of
    # python 2.7, so use 2.7.10 version.
    from collections2_7_10 import OrderedDict
else :
    from collections import OrderedDict # Only available in python 2.7 and later

assert sys.version_info[:2] == (2, 7), 'This is intended to run in Python 2.7.'

Tab = '\t'
NL  = '\n'

class MyException(Exception) :
    "Use for any error that shouldn't be handled, but can print OS-level string."
    pass
    
def stop(status = 1):
    sys.exit(status)
    
def raiser(ex): raise ex # Use this to raise an exception inside an expression

def err_msg(message) :
    print(message, file = sys.stderr)
    sys.stderr.flush()

def get_absolute_path(path) :
    "Convert ., .., and ~"
    return os.path.abspath(os.path.expanduser(path))

pjoin = os.path.join
dirname = os.path.dirname

def ls(dirName) :
    return os.listdir(get_absolute_path(dirName))

def pwd():
    return os.getcwd()

def cd(newDir, dump = False) :
    newDir = get_absolute_path(newDir)
    if dump :
        err_msg('cd ' + newDir)
    os.chdir(newDir)

chdir = cd

DirStack = []

def pushd(newDirName = None, dump = False) :
    DirStack.append(pwd())
    if newDirName != None :
        cd(newDirName, dump)

def popd(dump = False) :
    if len(DirStack) == 0 :
        raise MyException('popd: directory stack empty')
    cd(DirStack.pop(), dump)
    return pwd()

def dirs() :
    return [pwd()] + DirStack

def myopen(path, *args) :
    path = get_absolute_path(path)
    writing = len(args) > 0 and 'w' in args[0]
    if not writing and not file_exists(path) and file_exists(path + '.gz') :
        path = path + '.gz'
    if path[-3:] == '.gz' :
        import gzip
        return gzip.open(path, *args)
    elif path[-4:] == '.bgz' :
        import bgzf # Make sure sys.path has appropriate dir
        return bgzf.open(path, *args)
    else :
        return open(path, *args)

def file_exists(fileName) :
    return os.path.exists(get_absolute_path(fileName))

def is_directory(fileName) :
    return os.path.isdir(get_absolute_path(fileName))

isdir = is_directory

def mv(source, destination, dump = False, onlyDump = False) :
    if dump or onlyDump :
        err_msg('mv %s %s' % (source, destination))
    if onlyDump :
        return
    # Note that os.rename is similar to shutil.move but doesn't work
    #    if destination is an existing directory.
    shutil.move(get_absolute_path(source), get_absolute_path(destination))

def cp(source, destination, dump = False, onlyDump = False) :
    if dump or onlyDump :
        err_msg('cp %s %s' % (source, destination))
    if onlyDump :
        return
    shutil.copy(get_absolute_path(source), get_absolute_path(destination))

def ln_minus_s(fileName, linkName, dump = False, onlyDump = False) :
    linkName = get_absolute_path(linkName)
    if dump or onlyDump :
        err_msg('ln -s %s %s' % (fileName, linkName))
    if onlyDump :
        return
    os.symlink(fileName, get_absolute_path(linkName))

symlink = ln_minus_s

def rmdir(path, dump = False, onlyDump = False) :
    if dump or onlyDump :
        err_msg('rmdir %s' % path)
    if onlyDump :
        return
    os.rmdir(get_absolute_path(path))
    
def rm_minus_r(path, dump = False, onlyDump = False) :
    "Like rm -rf path"
    if dump or onlyDump :
        err_msg('rm -rf %s' % path)
    if onlyDump :
        return
    shutil.rmtree(path)
    
def rm(path, dump = False, onlyDump = False) :
    if dump or onlyDump :
        err_msg('rm %s' % path)
    if onlyDump :
        return
    os.remove(get_absolute_path(path))

def mkdir(path, dump = False, onlyDump = False) :
    if dump or onlyDump :
        err_msg('mkdir ' + path)
    if onlyDump :
        return
    return os.mkdir(get_absolute_path(path))

def assure_dir(path, dump = False) :
    if not is_directory(path) :
        mkdirs(path, dump)

def mkdirs(path, dump = False):
    "Create a directory and all parent directories."
    path = get_absolute_path(path)
    if os.path.isdir(path):
        pass
    elif os.path.isfile(path):
        raise OSError("cannot create directory, file already exists: '%s'" % path)
    else:
        head, tail = os.path.split(path)
        if head and not os.path.isdir(head):
            mkdirs(head, dump)
        if tail:
            mkdir(path, dump)

def same_file(path1, path2) :
    "Return True if the two paths point to the same file."
    return os.path.samefile(get_absolute_path(path1), get_absolute_path(path2))

def equal_files(fileName1, fileName2) :
    "Return True if the two files have equal contents."
    return checksum(fileName1) == checksum(fileName2)

def checksum(fileName) :
    """
    If file is too big to fit in memory, see this thread for an implementation that
    reads 4096-byte chunks and feeds them to md5 one at a time:
    https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    """
    return hashlib.md5(myopen(fileName).read()).hexdigest()

def temp_fileName(suffix = '', prefix = 'tmp', dir = None) :
    """
    Return a file name guaranteed to be unique, starting with prefix, ending with suffix,
        in dir or the default temporary directory.
    This just calls mkstemp and then closes and deletes the file it creates.
    """
    import tempfile
    fhandle, fname = tempfile.mkstemp(prefix = prefix, suffix = suffix, dir = dir)
    os.close(fhandle)
    rm(fname)
    return fname

def prepend_to_file(stringToPrepend, fileName) :
    # Prepend the string to the specified file
    tempFileName = temp_fileName()
    with myopen(tempFileName, 'w') as tempFile :
        tempFile.write(stringToPrepend)
        with myopen(fileName) as inFile :
            tempFile.write(inFile.read())
    mv(tempFileName, fileName)

def parse_arg(argName, shortArgName = None) :
    # If --argName is in sys.argv, return True and remove it.
    # If shortArgName (which should be a single character) is in an argument of the 
    #   form -string, where string doesn't start with -, then return True and remove it.
    # Otherwise return False
    if '--' + argName in sys.argv[1:] :
        del sys.argv[sys.argv.index('--' + argName, 1)]
        return True
    elif shortArgName != None :
        for pos, arg in enumerate(sys.argv) :
            if pos == 0 or arg[0] != '-' or (len(arg) > 1 and arg[1] == '-') :
                continue
            if shortArgName in arg[1:] :
                sys.argv[pos] = sys.argv[pos].replace(shortArgName, '')
                if sys.argv[pos] == '-' :
                    del sys.argv[pos]
                return True
    return False
        
def parse_arg_value(argName, argType = str, default = None) :
    # If sys.argv includes a string of the form --argName=value, 
    #   return the value, casted to type argType, and remove it from sys.argv.
    # Otherwise return default.
    for pos, arg in enumerate(sys.argv) :
        if pos == 0 or not arg.startswith('--' + argName + '=') :
            continue
        del sys.argv[pos]
        return argType(arg[len(argName) + 3 :])
    return default
                    
def without_trailing_digits(word) :
    while word != '' and word[-1].isdigit() :
        word = word[: -1]
    return word
remove_trailing_digits = without_trailing_digits
    
def without_trailing_char(string, char) :
    return string[:-1] if string.endswith(char) else string

def without_trailing_nl(line) :
    "Remove trailing \n, \r, or \r\n, if present"
    line = without_trailing_char(line, '\n')
    line = without_trailing_char(line, '\r')
    return line
    
remove_trailing_nl = without_trailing_nl
strip_nl = without_trailing_nl

def neighbors(iter) :
    # Return a list of the consecutive pairs of elements in the iterator.
    return list(neighbors_iter(iter))

def neighbors_iter(iter) :
    # Generate all consecutive pairs of elements in the iterator.
    firstTime = True
    for item in iter :
        if firstTime :
            firstTime = False
        else :
            yield (prev, item)
        prev = item

def antistrand(strand) :
    if strand == '+' :
        return '-'
    elif strand == '-' :
        return '+'
    else :
        return strand

def anti_interval(interval, strandOfInterval) :
    """
    Return interval on other strand that shares 3rd codon positions, and other strand.
    Explanation: if interval is coding, then the anti-interval often gets a "ghost"
        evolutionary signature of coding. Comparison might distinguish which is real.
    We could instead return an interval that shares the first two codon positions, i.e.,
        shift back 1 instead of forward 2. At first glance, that would make sense,
        because the coding evolutionary signature comes from suppression of substitutions
        in the first two positions rather than excess substitutions in the third position.
        However, when carefully considering the purpose of looking at the anti-interval,
        it seems like the right shift is slightly better. The goal is to distinguish
        the real signal from the ghost by comparing evolutionary signatures of the two. We
        want an antisense interval that will get a bad score if it is a ghost and a good
        one if it is real. In most cases it's a wash, but for the common case that the
        original interval ends just before a stop, there's a slight advantage to the 
        right shift. That's because in that case the anti-interval will be enriched
        for stops if it is a ghost (because a TAA or TAG codon on the original strand
        becomes TAx on that anti-strand), whereas these will have been suppressed if the
        anti-interval is coding.
    """
    antiOffset = 2 if strandOfInterval == '+' else -2
    return (type(interval)([interval[0] + antiOffset, interval[1] + antiOffset]),
            antistrand(strandOfInterval))

def anti_intervals(intervals, strandOfIntervals) :
    """
    Return intervals on other strand that share 3rd codon positions, and other strand.
    See explanation in anti_interval.
    Results at splice boundaries are silly.
    """
    return (type(intervals)([anti_interval(interval, strandOfIntervals)[0]
                             for interval in intervals]),
            antistrand(strandOfIntervals))

dnaComplementDict = {'G':'C', 'C':'G', 'A':'T', 'T':'A',
                     'g':'c', 'c':'g', 'a':'t', 't':'a'}
rnaComplementDict = {'G':'C', 'C':'G', 'A':'U', 'U':'A',
                     'g':'c', 'c':'g', 'a':'u', 'u':'a'}

def reverse_complement(string, rna = False) :
    """
    Reverse complement the string, complementing upper or lower case a, c, g, and t (or u 
        if rna is True), preserving case, and leaving all other characters the same, e.g., 
        'N' or '-'.
    Fail if string contains t or T if rna or if it contains u or U and not rna.
    """
    if rna :
        assert('t' not in string and 'T' not in string), string
        complementDict = rnaComplementDict
    else :
        assert('u' not in string and 'U' not in string), string
        complementDict = dnaComplementDict
    return ''.join(complementDict.get(char, char) for char in string[::-1])

def map_to_strand(string, strand, rna = False) :
    return reverse_complement(string, rna) if strand == '-' else string

def regionString_to_triples(regionStr) :
    # Parse a string of form chrom:START-END[+chrom:START-END]*, 
    #    where the numbers can include commas.
    # Return [(chrom, START, END)]
    triples = []
    for intervalStr in regionStr.split('+') :
        try :
            chrom, intervalPart = intervalStr.split(':')
            intervalPart = intervalPart.replace(',', '')
            uminus = intervalPart.startswith('-')
            if uminus :
                intervalPart = intervalPart[1:]
            start, end = map(int, intervalPart.split('-', 1)) # Keeps unary minus on end
            if uminus :
                start = -start
            triples.append((chrom, start, end))
        except ValueError, IndexError :
            raise MyException('Invalid interval string: %s' % intervalStr)
    return triples

def get_intervals_length(intervals, includesChrom = True) :
    # Return number of bases in a list of intervals: [(CHROM, START, END)...] 
    # (or [(START, END)...] if not includesChrom)
    if includesChrom :
        startIndex = 1
    else :
        startIndex = 0
    return sum(interval[startIndex + 1] - interval[startIndex] + 1
               for interval in intervals)
    
def intervals_prefix(intervals, strand, numBases, includesChrom = True) :
    # Return subset of intervals numBases long starting at beginning relative to strand.
    # Intervals is a list: [(CHROM, START, END)]
    # (or [(START, END)...] if not includesChrom)
    if numBases < 0 or numBases > get_intervals_length(intervals, includesChrom) :
        raise ValueError, ('intervals_prefix: invalid numBases (%s) for length %s' %
                            (numBases, get_intervals_length(intervals, includesChrom)))
    result = []
    reverse = strand == '-'
    startInd = [0, 1][includesChrom]
    for interval in intervals[::[1, -1][reverse]] :
        intervalLen = interval[startInd + 1] - interval[startInd] + 1
        if numBases >= intervalLen :
            result.insert([len(result), 0][reverse], interval)
            numBases -= intervalLen
            continue
        if numBases > 0 :
            if reverse :
                result.insert(0, interval[:startInd] +
                                 (interval[startInd + 1] - numBases + 1,
                                  interval[startInd + 1]))
            else :
                result.insert(len(result), interval[:startInd] +
                                           (interval[startInd],
                                            interval[startInd] + numBases - 1))
        return result
    return result

def intervals_suffix(intervals, strand, numBases, includesChrom = True) :
    # Return subset of intervals numBases long ending at end relative to strand.
    # Intervals is a list: [(CHROM, START, END)]
    # (or [(START, END)...] if not includesChrom)
    return intervals_prefix(intervals, ['-', '+'][strand == '-'], numBases, includesChrom)

def bed_line_to_intervals(line) :
    """
    Input a line in a BED format file.
    Return (name, chrom, intervals, strand),
       where intervals = [[start1, end1], [start2, end2], ...] in the input order, which
       should satisfy start1 <= end1 < start2 <= end2 < ...
    Note that the resulting intervals treat chromosomes as starting at position 1,
       and include both endpoints in the interval (as is done in gtf/gff files,
       the genome browser, and PhyloCSF), rather than treating chromosomes as starting at
       position 0 and ending beyond the end of the interval, as is done in input BED file.
       In other words, the starts will differ by 1 from the numbers in the BED file. 
       For example, an interval containing just the first base of a chromosome would be 
       0, 1 in the BED file, and 1, 1 here. (End is ignored for BED12, but not for BED6.)
    For Bed6 files, return one interval for the whole thing.
    """
    words = map(str.strip, line.split('\t')) # Remove trailing \n, as well as spaces
    chrom = words[0]
    firstStart = int(words[1]) + 1
    name = words[3]
    strand = words[5]
    if len(words) > 11 : # BED12
        numExons = int(words[9])
        exonSizes = map(int, without_trailing_char(words[10], ',').split(','))
        exonRelStarts = map(int, without_trailing_char(words[11], ',').split(','))
        intervals = [[firstStart + relStart, firstStart + relStart + size - 1]
                     for size, relStart in zip(exonSizes, exonRelStarts)]
        # Checks:
        if len(exonSizes) != numExons or len(exonRelStarts) != numExons :
            raise AssertionError, 'Mismatched counts: ' + line
        if int(words[2]) != intervals[-1][1] :
            raise AssertionError, 'Mismatched ends: ' + line
    elif len(words) <= 9 : # BED6
        firstEnd = int(words[2])
        intervals = [[firstStart, firstEnd]]
    else :
        raise AssertionError, 'Line with invalid number of fields: ' + line   
    return name, chrom, intervals, strand

def is_bed_comment(line) :
    # Return true if line from bed file is a comment or browser metadata 
    return line == '' or line[0] == '#' or line.split()[0] in ['browser', 'track']

def csv_to_columns(file) :
    # Read a comma-separated-value file.
    # Return a dictionary and a list: {columnName:[values]}, [columnNames in order].
    # columnNames and values will be stripped of leading and trailing spaces.
    valuesDict = None
    for line in file :
        if line == '' or line == '\n' or line[0] == '#' :
            continue # Skip comments and blank lines
        words = [word.strip() for word in line.split(',')] # strip removes trailing \n
        if valuesDict == None : # Header row
            valuesDict = dict((word, []) for word in words)
            columnNames = words
            numCols = len(words)
        else :
            numWords = len(words)
            if numWords < numCols :
                words += '' * (numCols - numWords)
            elif numWords > numCols :
                err_msg('read_csv: ignoring extra words.')
                del words[numCols :]
            for word, columnName in zip(words, columnNames) :
                valuesDict[columnName].append(word)
    return valuesDict, columnNames

def csv_to_table(file, indexColumnName) :
    # Read a comma-separated-value file.
    # Return a Table, where row names are taken from the column with name indexColumnName.
    # columnNames and values will be stripped of leading and trailing spaces.
    table = None
    for line in file :
        if line == '' or line == '\n' or line[0] == '#' :
            continue # Skip comments and blank lines
        words = [word.strip() for word in line.split(',')] # strip removes trailing \n
        if table == None : # Header row
            indexColumn = words.index(indexColumnName)
            table = Table(words[: indexColumn] + words[indexColumn + 1 :])
            numCols = len(words)
        else :
            numWords = len(words)
            if numWords < numCols :
                words += '' * (numCols - numWords)
            elif numWords > numCols :
                raise ValueError, 'csv_to_table: too many words (%d vs %d).' % (numWords,
                                                                                numCols)
            table.add_row(words[indexColumn],
                          words[: indexColumn] + words[indexColumn + 1 :])
    return table
    
def table_to_csv(file, table, indexColumnName = 'name') :
    # Print the Table as a comma-separated file (with index column first)
    columnNames = table.column_names()
    print(','.join([indexColumnName] + list(columnNames)), file = file)
    for rowName in table.row_names() :
        print(','.join([rowName] + [table[rowName][columnName]
                                    for columnName in columnNames]), file = file)

class Table(object) :
    """
    Class for storing a table with (ordered) row and column names so that a cell can be
        accessed as atable[rowName][columnName].
    Only access by methods -- implementation is subject to change.
    """
    def __init__(self, columnNames) :
        self.colNames = tuple(columnNames)
        self.colDict = dict((colName, colInd)
                            for colInd, colName in enumerate(self.colNames))
        self.rowNames = []
        self.rowDict = {} # {rowName : orderedValues, ...}
    def add_row(self, rowName, orderedValues) :
        # Values must be in order of column names
        self.rowNames.append(rowName)
        self.rowDict[rowName] = list(orderedValues)
    def delete_row(self, rowName) :
        self.rowNames.remove(rowName)
        del self.rowDict[rowName]
    def row_names(self) :
        return self.rowNames
    def column_names(self) :
        return self.colNames
    def __getitem__(self, rowName) :
        return _TableRow(self, self.rowDict[rowName])
        
class _TableRow(object) :
    # For internal use by Table. Implementation subject to change.
    def __init__(self, table, row) :
        self.table = table
        self.row = row
    def __getitem__(self, colName) :
        return self.row[self.table.colDict[colName]]
    def __setitem__(self, colName, value) :
        self.row[self.table.colDict[colName]] = value

def split(strToSplit, sep = None) :
    """
    Equivalent to strToSplit.split(sep) except that split('') -> [] instead of [''].
    Ideally split would be inverse of join: split(sep.join(listOfStrs), sep) -> listOfStrs
    However, sep.join([]) == sep.join(['']) == '', so it is impossible for split to always
        be inverse to join. Python makes the inverse work for [''] but not for []. This
        split function instead makes the inverse work for [] but not for [''].
    """
    return strToSplit.split(sep) if strToSplit != '' else []

def split_str(strToSplit, maxLength = 70) :
    """
    Split string into lines of length at most maxLength. Useful for writing fasta files.
    If maxLength is None, just return the string without splitting it.
    """
    if maxLength == None :
        yield strToSplit
        return
    start = 0
    while start < len(strToSplit) - maxLength :
        yield strToSplit[start : start + maxLength]
        start += maxLength
    yield strToSplit[start : ]
    
def write_to_fasta(outFile, header, sequence, maxLength = 70) :
    print('>' + header, file = outFile)
    for line in split_str(sequence, maxLength) :
        print(line, file = outFile)

def iter_fasta(faFile) :
    # Iterate through a fasta file returning pairs: sequence name, sequence
    seqName = None
    seq = ''
    for line in faFile :
        line = strip_nl(line)
        if len(line) == 0 :
            continue
        if line[0] == '>' :
            if seqName != None :
                yield seqName, seq
            seqName = line[1:]
            seq = ''
        else :
            seq += line
    if seqName != None :
        yield seqName, seq

def is_string_int(string) :
    try :
        int(string)
    except ValueError :
        return False
    return True

def is_string_number(s):
    try:
        float(s)
    except ValueError:
        return False
    return True

def get_backup_names(filePathAndName) :
    """
    Return the list of files that are backups of filePathAndName in increasing order,
        and ending with the (non-existent) file that would be the next backup.
    Backups are filePathAndName with .NNN inserted before .EXT or .EXT.gz or .EXT.bgz for
        whatever file extension is there, if there is one.
    """
    filePathAndName = get_absolute_path(filePathAndName)
    filePath, fileName = os.path.split(filePathAndName)
    fileNameNoZip = (fileName[:-3] if fileName.endswith('.gz')  else
                     fileName[:-4] if fileName.endswith('.bgz') else
                     fileName)
    dotPos = fileNameNoZip.rfind('.') if '.' in fileNameNoZip else len(fileNameNoZip)
    prefix = fileName[:dotPos]
    suffix = fileName[dotPos:]
    dirList = os.listdir(filePath)
    usedVers = [int(name[dotPos + 1 : dotPos + 4])
                    for name in dirList
                    if dotPos + 4 <= len(name) and
                       name[:dotPos] == prefix and
                       name[dotPos + 4:] == suffix and
                       is_string_int(name[dotPos + 1 : dotPos + 4])]
    usedVers.sort()
    nextVer = usedVers[-1] + 1 if usedVers else 1
    assert nextVer < 1000, 'version %d too high' % nextVer
    backupNames = [pjoin(filePath, prefix + '.%03d'% ver + suffix)
                   for ver in usedVers + [nextVer]]
    assert not file_exists(backupNames[-1]), backupNames[-1]
    return backupNames

def backup_prev_file(filePathAndName) :
    # If file exists, rename it to something of the form Prefix.NNN.Suffix
    filePathAndName = get_absolute_path(filePathAndName)
    if file_exists(filePathAndName) and filePathAndName != '/dev/null' :
        backupNames = get_backup_names(filePathAndName)
        mv(filePathAndName, backupNames[-1])

def open_with_backup(filePathAndName, additionalMode = '') :
    # Opens the file for writing. Takes any file with the same name and renames it
    # to something of the form Prefix.NNN.Suffix
    backup_prev_file(filePathAndName)
    return myopen(filePathAndName, 'w' + additionalMode)

###### DictClass ############

class DictClass(OrderedDict) :
    # Exactly like an OrderedDict, except that you can also access using . syntax.
    # E.g. DictClass({'a':7}).a -> 7
    def __init__(self, *pArgs, **kArgs) :
        # Record that we are initializing when calling OrderedDict.__init__ so that when
        # it creates members, __setattr__ treats them as attributes, not as items.
        self.__dict__['__initializing_the_OrderedDict'] = True
        super(DictClass, self).__init__(*pArgs, **kArgs)
        self.__dict__['__initializing_the_OrderedDict'] = False
    def __getattr__(self, memberName) :
        try :
            # First see if it is a real attribute (e.g., of OrderedDict); if not then
            #     treat items of the OrderedDict as if they were attributes.
            return self.__dict__.get(memberName, self[memberName])
        except KeyError :
            raise AttributeError('%r object has no attribute %r' % (type(self).__name__,
                                                                    memberName))
    def __setattr__(self, memberName, value) :
        # Warning: If OrderedDict tries to create new members after __init__
        #     this won't work because we'll make them items instead of members.
        if self.__dict__['__initializing_the_OrderedDict'] :
            self.__dict__[memberName] = value
        else :
            self[memberName] = value
