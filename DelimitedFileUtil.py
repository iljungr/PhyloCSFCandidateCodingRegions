#!/usr/bin/env python
# Copyright 2019 Irwin Jungreis
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,x
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""
DelimitedFile.py
Utilities for working with delimited files.
"""
from __future__ import division
from __future__ import print_function
import sys, itertools
from CommonUtils import DictClass, myopen, open_with_backup, strip_nl

assert sys.version_info[:2] == (2, 7), 'This is intended to run in Python 2.7.'

class DelimitedFileReader(object) :
    """
    Iterator for reading a delimited file, where fields are separated by the delimiter
        (e.g. Tab) and (optionally) the first (header) line contains the field names.
        Lines starting with '#' may be skipped.
    After opening the file, and, if desired, specifying mapping functions for particular
        fields, call next (or iterate through the class instance), to read the lines after
        the header.
    Tip: to get a list of all the records, call list(reader).
    """
    def __init__(self, fileName, delimiter = '\t', fieldNames = None,
                 intFields = None, floatFields = None, skipComments = True,
                 allowDupFields = False, skipQuotes = False, skipEmptyLines = False) :
        """ If fieldNames is None, get from the first line of the file;
                o.w., assume there is no header line.
            Integer and float fields may be added later using add_mapper instead of here.
            If skipComments, ignore lines that start with '#'.
            If allowDupFields, allow several fields to have the same name; in that case
                the result of accessing that field is ambiguous.
            If skipQuotes, strip any double quote characters from start or end of string.
                Excel sometimes adds these when saving Tab delimited files.
            If skipEmptyLines, allow blank lines, and ignore them.
        """
        self.inFile = myopen(fileName)
        self.fileIter = (_skip_pounds(self.inFile, skipEmptyLines) if skipComments else
                         self.inFile)
        self.delimiter = delimiter
        if fieldNames == None :
            self.prevLine = self.fileIter.next()
            self.fieldNames = strip_nl(self.prevLine).split(self.delimiter)
        else :
            self.prevLine = None
            self.fieldNames = list(fieldNames)
        if not allowDupFields :
            for name in self.fieldNames :
                assert sum(n == name for n in self.fieldNames) == 1, 'Dup field %s' % name
        self.fieldMappers = {}
        if intFields != None :
            self.add_mapper(intFields, int)
        if floatFields != None :
            self.add_mapper(floatFields, float)
        self.skipQuotes = skipQuotes
        self.skipEmptyLines = skipEmptyLines
    def get_fieldNames(self) :
        return self.fieldNames[:] # Return a copy so changing it won't affect this DFR
    def add_mapper(self, fieldNameOrNames, mapper) :
        # fieldNameOrNames can be a string or a collection of strings
        # mapper = lambda valueString : value, e.g., int
        # If mapper can't map it should raise ValueError so field will stay str.
        if isinstance(fieldNameOrNames, str) :
            fieldNameOrNames = [fieldNameOrNames]
        for fieldName in fieldNameOrNames :
            self.fieldMappers[fieldName] = mapper

    class DelimitedFileRecord(DictClass) :
        """
        This class holds the information from one line of a delimited file. The fields can
            be accessed using [fieldName] as if it were a dictionary,
            or using the .fieldName syntax, as if it were a class.
        The fields are mapped at access time, not when the record is created (which is why
            we don't just use a DictClass). Reasons for this delayed mapping:
            - So when we write, we can write the unmapped fields. That way 1.10 doesn't
              become 1.1.
            - If a value is changed, the mapping will be applied to the new value.
        """
        def __init__(self, itemsBeforeMap, fieldMappers) :
            DictClass.__init__(self, itemsBeforeMap)
            # self.fieldMappers = fieldMappers would be wrong because DictClass
            #     overrided __setattr__.
            object.__setattr__(self, 'fieldMappers', fieldMappers)
        def __getitem__(self, key) :
            return self._map_value(key, DictClass.__getitem__(self, key))
        def items(self) :
            "Warning: won't work in 3.0 because it returns a list intead of an iteration."
            return list(self._yielditems())
        def values(self) :
            "Warning: won't work in 3.0 because it returns a list intead of an iteration."
            return [value for key, value in self._yielditems()]
        def _yielditems(self) :
            for key, value in DictClass.items(self) :
                yield key, self._map_value(key, value)
        def _map_value(self, key, value) :
            if key in self.fieldMappers :
                try :
                    return self.fieldMappers[key](value) # OK; . only calls __getattr__
                                                         # for nonmembers.
                except ValueError : # int('') -> ValueError
                    pass
                except TypeError : # int(None) -> TypeError. None can't be read but could
                                   # be added later
                    pass
            return value # Keep as str things that can't be mapped.
        def get(self, key, default = None) :
            # Override "get" so that it does mapping
            return self[key] if self.has_key(key) else default

    def next(self) :
        "Read a line and return a DictClass like {fieldName : mappedFieldValue, ...}"
        self.prevLine = self.fileIter.next()
        fieldValues = strip_nl(self.prevLine).split(self.delimiter)
        if self.skipQuotes :
            fieldValues = [value.strip('"') for value in fieldValues]
        items = itertools.izip_longest(# Handle missing fields at the end
                                       self.fieldNames, fieldValues, fillvalue = '')
        return self.DelimitedFileRecord(items, self.fieldMappers)
    def get_prev_line(self) : # Return the most recently read line
        return self.prevLine
    def __iter__(self) :
        return self
    def close(self) :
        self.inFile.close()

DFR = DelimitedFileReader

class DelimitedFileWriter(object) :
    """
    Class for writing a delimited file.
    """
    def __init__(self, fileName, fieldNames, delimiter = '\t', writeHeader = True,
                 backup = False) :
        self.outFile = open_with_backup(fileName) if backup else myopen(fileName, 'w')
        self.delimiter = delimiter
        self.fieldNames = fieldNames[:] # Make a copy
        for name in self.fieldNames :
            assert sum(n == name for n in self.fieldNames) == 1, 'Duplicate field %s.' % name
        if writeHeader :
            print(delimiter.join(fieldNames), file = self.outFile)
    def write_line(self, rec) :
        "rec is a DelimitedFileRecord, DictClass, or dict with the values for each field."
        "We use dict.__getitem__ to bypass mapping, so, e.g., 1.10 stays 1.10."
        def get_it(fieldName) :
            try :
                value = dict.__getitem__(rec, fieldName)
                if type(value) is unicode :
                    # Convert to str, since no encoding was specified when opening outFile
                    value = value.encode('utf8')
                    # Maybe we should instead convert unicode to str with ? for non-ascii
                    # chars via value = value.encode('ascii', 'replace') ...
                else :
                    value = str(value)
                return value
            except KeyError :
                return ''
        print(self.delimiter.join(map(get_it, self.fieldNames)), file = self.outFile)
    def write_comment(self, comment) :
        "Write a line containing comment preceded by '#'"
        print('#' + comment, file = self.outFile)
    def close(self) :
        self.outFile.close()
    def flush(self) :
        self.outFile.flush()

DFW = DelimitedFileWriter

def _skip_pounds(strIter, skipEmptyLines = False) :
    """
    Return an iterator like strIter except it skips strings that start with "#".
    If skipEmptyLines, also skip lines consisting of nothing but a newline character.
    """
    return itertools.ifilterfalse(lambda s : s.startswith('#') or
                                             (len(strip_nl(s)) == 0 if skipEmptyLines else
                                              False),
                                  strIter)

