#!/usr/bin/env python
import sys
import os
import re

CPP_SUFFIXES = ['.cpp','.cc','.hpp','.h']

'''
Example usage:

    ./scoped_lock_checker.py [foldername1 [foldername 2 ... ]]


A common error with scoped locks is to construct the scoped lock, but
fail to assign it to a variable.  Eg.,

   Mutex m;
   scoped_lock(m);

instead of

   Mutex m;
   scoped_lock sl (m);

In the first case, the block of code instantly constructs and destroys
a scoped lock object.  Subsequent code that is executed is unguarded
by the mutex.  In contrast, in the second example, code executed after
constructing the scoped lock is guarded until sl goes out of scope.

This script tries to check for this error.  It accepts a series of
folders and examines all c++ files in that folder and its subfolders
looking for the error described above.


Notes:

   1) Script probably does not handle symbolic links correctly, and
      may go into an infinite loop if it encounters links.

   2) False positives: the script may incorrectly identify an error in
      a few particular cases:

        * A programmer names a variable/method strangely.  Eg.,
        
            void do_something_while_in_scoped_lock();
            int holder_of_scoped_lock(3);

        * Strange, but valid line-breaking.  Eg.,
        
            scoped_lock l =
                scoped_lock(m);
            
        * Constructors and destructors for the scoped_lock class as
          errors.  In general though, you probably should only be
          running this on your application code, rather than its
          associated libraries.

'''



class Error(object):
    '''
    Just holds information necessary for reporting info back to client.
    '''
    def __init__(self,filename,incorrect_line):
        self.filename = filename
        self.incorrect_line = incorrect_line

    def __str__(self):
        return 'Error in %s: %s' % (self.filename,self.incorrect_line)

def list_cpp_files(file_list):
    '''Finds all c++ files in file_list and returns them.
    
    Args:
        file_list: {list} each element is a filename string.
    
    Returns:

        A list of strings from file_list whose suffixes indicate
        the file is a c++ file (cc, cpp, h, hpp)
    '''
    to_return = []
    for filename in file_list:
        for cpp_suffix in CPP_SUFFIXES:
            if filename.endswith(cpp_suffix):
                to_return.append(filename)
                break
    return to_return

def run_file(filename):
    '''Searches for incorrect uses of scoped_lock in program text.

    Args:
        filename: {string} Fully-qualified name of file to search for
        potential scoped lock errors in program.
    
    Returns:
        {list}: Each element of list is an Error object of potential
        scoped lock errors found in this file.
    '''
    regex = '.*?scoped_lock\s*\('
    errors = []
    
    with open(filename,'r') as f:
        contents = f.read()
        questionable_instances = re.findall(regex,contents)
        error_text = filter(
            lambda questionable:
                None if ('=' in questionable) else questionable,
            questionable_instances)

        errors += map(
            lambda text: Error(filename,text),
            error_text)
        
    return errors
        
def run_folder(foldername):
    '''Searches for incorrect uses of scoped_lock in all cpp files in
    folder.

    Args:
        foldername: {string} For all c++ files in foldername and its
        subdirectories, check if they might have scoped lock errors.
    
    Returns:
        {list}: Each element of list is an Error object of potential
        scoped lock errors found in this folder or its subfolders.
    '''
    folder_errors = []
    for root,dirs,files in os.walk(foldername):
        # check files in this directory
        cpp_files = list_cpp_files(files)
        for filename in cpp_files:
            fq_filename = os.path.join(root,filename)
            folder_errors += run_file(fq_filename)
            
    return folder_errors

def run_list_of_folders(folder_list):
    '''Searches for incorrect uses of scoped_lock in all cpp files in
    each folder in folder_list and its subdirectories.  Prints
    potential errors.

    Args:
        foldername: {list} For all c++ files in every folder specified
        in folder_list and their subdirectories, check if they might
        have scoped lock errors.
    
    Returns:
        None
    '''
    all_errors = []
    for folder in folder_list:
        all_errors += run_folder(folder)

    if len(all_errors) == 0:
        print '\nNo errors\n'
    else:
        print '\nErrors found!\n'
        for error in all_errors:
            print (str(error))



if __name__ == '__main__':
    folder_list = sys.argv[1:]
    run_list_of_folders(folder_list)
