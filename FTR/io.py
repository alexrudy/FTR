# -*- coding: utf-8 -*-
"""
Base classes for system IO
"""

import os

class IOBase(object):
    """A base class for Input/Output of filters"""
    
    def write(self, filename, **kwargs):
        """Write to a file with a given name.
        
        :param filename: The name of the file to write.
        :keyword format: The file format to use.
        
        The file format will be deduced from the extension if it isn't explicitly provided.
        
        """
        format = kwargs.pop('format', os.path.splitext(filename)[1][1:])
        
        function = "__to_{0:s}__".format(format)
        if not hasattr(self, function):
            raise ValueError("Unknown write format '{0}' for {1}".format(
                format, self.__class__.__name__))
        return getattr(self, function)(filename, **kwargs)
        
    @classmethod
    def read(cls, filename, **kwargs):
        """Read from a file with a given name.
        
        :param filename: The name of the file to read.
        :keyword format: The file format to use.
        
        The file format will be deduced from the extension if it isn't explicitly provided.
        
        """
        format = kwargs.pop('format', os.path.splitext(filename)[1][1:])
        
        function = "__from_{0:s}__".format(format)
        if not hasattr(cls, function):
            raise ValueError("Unknown read format '{0}' for {1}".format(
                format, cls.__name__))
        return getattr(cls, function)(filename, **kwargs)
        
    @classmethod
    def formats(cls, readonly=False, writeonly=False):
        """Return a dictionary of formats, with tuples for read/write."""
        formats = {}
        for method in dir(cls):
            if not (method.endswith("__") and len(method) > 7):
                continue
            if method.startswith("__to_"):
                format = method[5:-2]
                read = formats.get(format, (False, False))[0]
                formats[format] = (read, True)
            elif method.startswith("__from_") and len(method) > 9:
                format = method[7:-2]
                write = formats.get(format, (False, False))[1]
                formats[format] = (True, write)
            
        
        # Hanlde the return value.
        if readonly:
            readers = []
            for key, (read, write) in formats.values():
                if read:
                    readers.append(key)
            formats = readers
        elif writeonly:
            writers = []
            for key, (read, write) in formats.values():
                if write:
                    writers.append(key)
            formats = writers
        
        return formats
                