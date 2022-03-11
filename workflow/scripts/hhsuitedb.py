#!/usr/bin/env python3

"""
    hhbsuite.py
    Creates HH-suite database files from A3M and HHM files 
    Usage: Usage: python hhsuite_db.py -o <db_name> [-ia3m <a3m_dir>] [-ihhm <hhm_dir>] [-ics <cs_dir>] [more_options]

    HHsuite version 3.0.0 (15-03-2015)

    Reference: 
    Remmert M., Biegert A., Hauser A., and Soding J.
    HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
    Nat. Methods, epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

    (C) Johannes Soeding, 2012

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    We are very grateful for bug reports! Please contact us at soeding@mpibpc.mpg.de
"""


from glob import glob
import subprocess
import shlex

import tempfile
import os
import sys
import shutil

# Get the script folder of hhsuite
HHLIB = ""

for path in os.environ['PATH'].split(os.pathsep):
    if os.path.isfile(os.path.join(path.replace('bin', 'scripts'), 'hhsuitedb.py')):
        HHSCRIPT = path.replace('bin', 'scripts')
        HHLIB = os.path.dirname(path)

sys.path.insert(0, HHSCRIPT)

import ffindex
import a3m

HHM_FOCUS = snakemake.params.comparison

##########################################################################
sys.stderr = sys.stdout = open(snakemake.log[0], "w")
##########################################################################

def write_glob_to_file(glob_expr, filename):
    """Write a glob expression to a file .

    Args:
        glob_expr (string): Global expression to write
        filename (string)): path to file
    """
    fh = open(filename, "w")
    for f in glob(glob_expr):
        fh.write(f+"\n")
    fh.close()

##########################################################################

def write_set_to_file(set, filename):
    """Write a set of strings to a file .

    Args:
        set (set): set of file
        filename (string): path to file
    """
    fh = open(filename, "w")
    for f in set:
        fh.write(f+"\n")
    fh.close()

##########################################################################

def execute(command, env={}):
    """Execute a shell command and return its output .

    Args:
        command (string): Sheel command to run
        env (dict): optional, dict of environment variable (default: {})
    Returns:
        bytes: Output and error of the command
    """

    print(f'Executing {command}')
    cmd = shlex.split(command)
    if env:
        process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
    else:
        process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return stdout, stderr

##########################################################################

def remove_files_from_index(files_file, database_index_path):
    """Remove files from the database

    Args:
        files_file (list of string): Files to remove from the database
        database_index_path (string): Path to the database
    """

    if (os.path.exists(database_index_path) and os.path.getsize(database_index_path) > 0):
        stdout, stderr = execute(f"ffindex_modify -us -f {files_file}, {database_index_path}")
        print(f"----ffindex_modify - stdout----\n{stdout.decode('utf8')}\n----ffindex_modify - stderr----\n{stderr.decode('utf8')}\n")

##########################################################################

def build_ffindex_database(files_file, data_path, index_path):
    """Build the database for the ffindex database .

    Args:
        files_file (list of string): Files to remove from the database
        data_path (string): path to the ffdata file
        index_path (string): path to the ffindex file
    """
    stdout, stderr = execute(f"ffindex_build -s {data_path} {index_path} -f {files_file}")
    print(f"----ffindex_build - stdout----\n{stdout.decode('utf8')}\n----ffindex_build - stderr----\n{stderr.decode('utf8')}\n")

##########################################################################

def optimize_database(database_data_path, database_index_path):
    """Optimize the database using the size of the file and the size of the database .

    Args:
        database_data_path (string): path to the ffdata file
        database_index_path (string): path to the ffindex file
    """

    test_data = os.path.exists(database_data_path)
    size_data = os.path.getsize(database_data_path)
    test_index = os.path.exists(database_index_path)
    size_index = os.path.getsize(database_index_path)

    if test_data and test_index and size_data > 0 and size_index > 0: 
        stdout, stderr = execute(f"ffindex_build -as {database_data_path}.optimized {database_index_path}.optimized -d {database_data_path} -i {database_index_path}")
        print(f"----ffindex_build - stdout----\n{stdout.decode('utf8')}\n----ffindex_build - stderr----\n{stderr.decode('utf8')}\n")
        shutil.move(database_data_path+".optimized", database_data_path)
        shutil.move(database_index_path+".optimized", database_index_path)

##########################################################################

def get_large_a3ms(a3m_base_path):
    """Reads the a3m database and returns a list of names that are large enough for large alignment files .

    Args:
        a3m_base_path (string): Path to the a3m database

    Returns:
        [type]: [description]
    """

    entries = ffindex.read_index(f"{a3m_base_path}.ffindex")
    data = ffindex.read_data(f"{a3m_base_path}.ffdata")
  
    large_alignments = set()
    for entry in entries:
        lines = ffindex.read_lines(entry, data)
        alignment = a3m.A3M_Container()
        try:
            alignment.read_a3m_from_lines(lines)
      
            if alignment.number_sequences > 50:
                large_alignments.add(entry.name)
        except:
            print(f"Warning: A3M {entry.name} is corrupted!\n")

    if len(entries) > 0 and len(large_alignments) == 0:
        large_alignments.add(entries[0].name)

    return large_alignments

##########################################################################

def calculate_hhm(threads, a3m_base_path, hhm_base_path):
    """Calculate hhm using a3m

    Args:
        threads (int): Number of cpus to use
        a3m_base_path (string): Path to a3m database
        hhm_base_path (string): Path to hhm database
    """

    tmp_dir = tempfile.mkdtemp()
    
    large_a3ms = get_large_a3ms(a3m_base_path)
    large_a3m_index = os.path.join(tmp_dir, "large.ffindex")

    a3m_index = read_ffindex(f"{a3m_base_path}.ffindex")
    write_subset_index(a3m_index, large_a3ms, large_a3m_index)
    
    stdout, stderr = execute(f"mpirun -np {threads} ffindex_apply {a3m_base_path}.ffdata {large_a3m_index} -d {hhm_base_path}.ffdata -i {hhm_base_path}.ffindex -- hhmake -v 0 -i stdin -o stdout")
    print(f"----ffindex_apply + hhmake - stdout----\n{stdout.decode('utf8')}\n----ffindex_build + hhmake - stderr----\n{stderr.decode('utf8')}\n")

##########################################################################

def calculate_cs219(a3m_db_base, cs_db_base):
    """Calculate CSU219 using the HHLIB framework .

    Args:
        threads (int): number of cpu to use
        a3m_db_base (string): Path to the a3m database
        cs_db_base (string): Path to the hhr database
    """
    hhlib_environment = HHLIB if HHLIB else os.environ['HHLIB']
    print(hhlib_environment)
    if(hhlib_environment):
        #TODO: check
        path_cs219_lib = os.path.join(hhlib_environment, "data", "cs219.lib")
        path_context_data = os.path.join(hhlib_environment, "data", "context_data.lib")
        stdout, stderr = execute(f"cstranslate -A {path_cs219_lib} -D {path_context_data} -x 0.3 -c 4 --ffindex -i {a3m_db_base} -o {cs_db_base} -I a3m -b")
        print(f"----cstranslate - stdout----\n{stdout.decode('utf8')}\n----cstranslate - stderr----\n{stderr.decode('utf8')}\n")
    else:
        print("ERROR: HHLIB environment variable not set! See manual!\n")
        sys.exit(1)
        #TODO throw error

##########################################################################

def merge_databases(source_data_path, source_index_path, dest_data_path, dest_index_path):
    """Merge two databases into one .

    Args:
        source_data_path (string): Path to the source database ffdata
        source_index_path (string): Path to the source database ffindex
        dest_data_path (string): Path to the dest database ffdata
        dest_index_path (string): Path to the dest database ffindex
    """
    stdout, stderr = execute(f"ffindex_build -as -d {source_data_path} -i {source_index_path} {dest_data_path} {dest_index_path}")
    print(f"----ffindex_build - stdout----\n{stdout.decode('utf8')}\n----ffindex_build - stderr----\n{stderr.decode('utf8')}\n")

##########################################################################

def sort_database(data_path, index_path):
    """Sort database using ffindex .

    Args:
        data_path (string): Path to the database ffdata
        index_path (string): Path to the database ffindex
    """
    stdout, stderr = execute(f"ffindex_build -as {data_path} {index_path}")
    print(f"----ffindex_build - stdout----\n{stdout.decode('utf8')}\n----ffindex_build - stderr----\n{stderr.decode('utf8')}\n")

##########################################################################

def read_ffindex(file):
    """Read a ffindex and return a list of all the lines in the file .

    Args:
        file (string): path to the ffindex

    Returns:
        list of string: The file read line by line
    """
    fh = open(file, "r")
    
    index = []
    for line in fh:
        index.append(line.rstrip().split())
    
    fh.close()
    
    return index

##########################################################################

def write_subset_index(index, subset, output_file):
    """Write index in a subset of a subset .

    Args:
        index (list of ffindex): the ffindex value to write
        subset (list of string): list of value to write
        output_file (string): new ffindex file
    """

    fh = open(output_file, "w")
    for entry in index:
        if entry[0] in subset:
            name = entry[0]
            offset = entry[1]
            length = entry[2]
            fh.write(f"{name:.64}\t{offset}\t{length}\n")
    fh.close()

##########################################################################

def is_sorted(index):
    """Check if a list is sorted by the index .

    Args:
        index (list): All lines in the ffindex

    Returns:
        [type]: [description]
    """
    for i in range(len(index)-1):
        if(index[i][0] > index[i+1][0]):
            return False
    return True

##########################################################################

def get_duplicates(index):
    """Returns a list of duplicates in the index .

    Args:
        index (list): All lines in the ffindex

    Returns:
        list: All the duplicated line
    """
    duplicates = set()
    entries = set()
    
    for entry in index:
        if(entry[0] in entries):
            duplicates.add(entry[0])
        
        entries.add(entry[0])

    return duplicates

##########################################################################

def get_missing(index, cmp_index):
    """Get the set of entries that are not in the index and the index that are not in the index .

    Args:
        index (list): Lines on the ffindex
        cmp_index (list): Lines in the other ffindex

    Returns:
        list: List of missing
    """
    missing = set()
    
    index_set = set()
    for entry in index:
        index_set.add(entry[0])
        
    cmp_index_set = set()
    for entry in cmp_index:
        cmp_index_set.add(entry[0])
    
    for key in cmp_index_set:
        if key not in index_set:
            missing.add(key)
    
    return missing

##########################################################################

def get_overhead(index, cmp_index):
    """Get the overhead between two ffindex files.

    Args:
        index (list): List of line in the ffindex
        cmp_index (list): List of line in the ffindex

    Returns:
        list: List of overheads
    """
    overhead = set()
    
    index_set = set()
    for entry in index:
        index_set.add(entry[0])
        
    cmp_index_set = set()
    for entry in cmp_index:
        cmp_index_set.add(entry[0])
    
    for key in index_set:
        if key not in cmp_index_set:
            overhead.add(key)
    
    return overhead

##########################################################################

def handle_duplicates(suffix, calculate, threads, db_basename, force_mode):
    """This function handles duplicates and writes a temporary file to a temporary file .

    Args:
        suffix (string): type of data
        calculate (function): calculate function depending on the data
        threads (int): number of cpus
        db_basename (string): Path to the database folder with the prefix
        force_mode (bool): if value need to be forced or not
    """
    index = read_ffindex(f"{db_basename}_{suffix}.ffindex")
    duplicates = get_duplicates(index)
    
    if(suffix == "a3m" and len(duplicates) > 0):
        sys.stderr.write(f"ERROR: {db_basename}_a3m.ffindex contains duplicates!\n")
        sys.stderr.write("ERROR: Your database is broken!\n")
        exit(1);
        
    if(len(duplicates) == 0):
        return
    
    for duplicate in duplicates:
        sys.stderr.write(f"WARNING: {db_basename}_{suffix}.ffindex contains duplicate {duplicate}!\n")

    
    if force_mode:
        tmp_dir = tempfile.mkdtemp()
        
        try:
            sys.stderr.write("WARNING: remove duplicates and recalculate them from their corresponding a3m's!\n")
            a3m_index = read_ffindex(f"{db_basename}_a3m.ffindex")
            
            duplicates_a3m_base = os.path.join(tmp_dir, "duplicates_a3m")
            duplicates_a3m_index_file = os.path.join(tmp_dir, "duplicates_a3m.ffindex")
            duplicates_a3m_data_file = os.path.join(tmp_dir, "duplicates_a3m.ffdata")
            write_subset_index(a3m_index, duplicates, duplicates_a3m_index_file)
            os.symlink(f"{db_basename}_a3m.ffdata", duplicates_a3m_data_file)
            
            duplicates_new_base = os.path.join(tmp_dir, "duplicates_new")
            duplicates_new_index_file = os.path.join(tmp_dir, "duplicates_new.ffindex")
            duplicates_new_data_file = os.path.join(tmp_dir, "duplicates_new.ffdata")
            
            duplicates_index_file = os.path.join(tmp_dir, "duplicates.dat")
            write_set_to_file(duplicates, duplicates_index_file)
            remove_files_from_index(duplicates_index_file, f"{db_basename}_{suffix}.ffindex")
            
            calculate(duplicates_a3m_base, duplicates_new_base)
            sort_database(duplicates_new_data_file, duplicates_new_index_file)
            
            merge_databases(duplicates_new_data_file, duplicates_new_index_file, 
                            f"{db_basename}_cs219.ffdata", f"{db_basename}_cs219.ffindex")
            sort_database(f"{db_basename}_cs219.ffdata", f"{db_basename}_cs219.ffindex")
            optimize_database(f"{db_basename}_cs219.ffdata", f"{db_basename}_cs219.ffindex")
        finally:
            shutil.rmtree(tmp_dir)
    else:
        sys.stderr.write("You may try to use the option --force to fix the database!\n")

##########################################################################
    
def handle_unsorted(suffix, db_basename, force_mode):
    """Check if an index is sorted and sort it.

    Args:
        suffix (string): name of the database suffix
        db_basename (string): name of the database
        force_mode (bool): if you want to force mode on
    """
    index = read_ffindex(f"{db_basename}_{suffix}.ffindex")

    if(not is_sorted(index)):
        sys.stderr.write(f"Index {db_basename}_{suffix}.ffindex is unsorted!\n")
        if force_mode:
            sys.stderr.write(f"Try to sort unsorted index {db_basename}_{suffix}.ffindex!\n")
            sort_database(f"{db_basename}_{suffix}.ffdata", f"{db_basename}_{suffix}.ffindex")
        else:
            sys.stderr.write("You may try to use the option --force to fix the database!\n")

##########################################################################
    
def handle_missing(suffix, calculate, db_basename, threads, force_mode):
    """Calculate missing entries and write them to a temporary file .

    Args:
        suffix (string): type of data
        calculate (function): calculate function depending on the data
        db_basename (string): Path to the database folder with the prefix
        threads (int): number of cpus
        force_mode (bool): if value need to be forced or not
    """
    index = read_ffindex(f"{db_basename}_{suffix}.ffindex")
    a3m_index = read_ffindex(f"{db_basename}_a3m.ffindex")
    
    missing = get_missing(index, a3m_index)
    
    if(len(missing) == 0):
        return
    
    for f in missing:
        sys.stderr.write(f"WARNING: Missing entry {f} in {db_basename}_{suffix}.ff{{data,index}}!\n")
    
    if force_mode:
        sys.stderr.write("WARNING: Try to calculate missing entries!\n")
        tmp_dir = tempfile.mkdtemp()
        
        try:
            missing_base = os.path.join(tmp_dir, f"missing_{suffix}")
            missing_index_file = os.path.join(tmp_dir, f"missing_{suffix}.ffindex")
            missing_data_file = os.path.join(tmp_dir, f"missing_{suffix}.ffdata")
            
            missing_a3m_base = os.path.join(tmp_dir, "missing_a3m")
            missing_a3m_index_file = os.path.join(tmp_dir, "missing_a3m.ffindex")
            missing_a3m_data_file = os.path.join(tmp_dir, "missing_a3m.ffdata")
            write_subset_index(a3m_index, missing, missing_a3m_index_file)
            os.symlink(os.path.abspath(f"{db_basename}_a3m.ffdata"), missing_a3m_data_file)
            
            calculate(missing_a3m_base, missing_base)
            
            merge_databases(missing_data_file, missing_index_file, f"{db_basename}_{suffix}.ffdata", f"{db_basename}_{suffix}.ffindex")
            optimize_database(f"{db_basename}_{suffix}.ffdata", f"{db_basename}_{suffix}.ffindex")
        finally:
            shutil.rmtree(tmp_dir)
    else:
        sys.stderr.write("You may try to use the option --force to fix the database!\n")

##########################################################################

def handle_overhead(suffix, db_basename, force_mode):
    """Fix the overhead entries in a3m database .

    Args:
        suffix (string): name of the database suffix
        db_basename (string): name of the database
        force_mode (bool): if you want to force mode on
    """
    index = read_ffindex(f"{db_basename}_{suffix}.ffindex")
    a3m_index = read_ffindex(f"{db_basename}_a3m.ffindex")

    overhead = get_overhead(index, a3m_index)

    #delete overhead cs219 files
    if(len(overhead) == 0):
        return

    for f in overhead:
        sys.stderr.write(f"WARNING: Entry {f} from {db_basename}_{suffix}.ff{{data,index}} has no corresponding entry in the a3m database!\n")

    if force_mode:
        sys.stderr.write("WARNING: Try to fix overhead entries!\n")
        tmp_dir = tempfile.mkdtemp()
        
        try:
            index_file = os.path.join(tmp_dir, "to_delete.dat")
            write_set_to_file(overhead, index_file)
            remove_files_from_index(index_file, f"{db_basename}_{suffix}.ffindex")
            optimize_database(f"{db_basename}_{suffix}.ffdata", f"{db_basename}_{suffix}.ffindex")
        finally:
            shutil.rmtree(tmp_dir)
    else:
        sys.stderr.write("You may try to use the option --force to fix the database!\n")

##########################################################################

def check_a3m_format(db_basename, force_mode):
    """Check the format of the A3M file for corrupted data .

    Args:
        db_basename (string): name of the database
        force_mode (bool): if you want to force mode on
    """
    entries = ffindex.read_index(f"{db_basename}_a3m.ffindex")
    data = ffindex.read_data(f"{db_basename}_a3m.ffdata")
    
    corrupted_alignments = set()
    for entry in entries:
        lines = ffindex.read_lines(entry, data)
        alignment = a3m.A3M_Container()
        try:
            alignment.read_a3m_from_lines(lines)
        except:
            corrupted_alignments.add(entry.name)
            sys.stderr.write(f"Warning: A3M {entry.name} is corrupted!\n")
    
    if len(corrupted_alignments) == 0:
        return
    
    if force_mode:
        tmp_dir = tempfile.mkdtemp()
        
        try:
            sys.stderr.write("WARNING: remove corrupted a3m's!\n")
            
            corrupted_index_file = os.path.join(tmp_dir, "corrupted.dat")
            write_set_to_file(corrupted_alignments, corrupted_index_file)
            
            for suffix in ["a3m", "cs219", "hhm"]:
                remove_files_from_index(corrupted_index_file, f"{db_basename}_{suffix}.ffindex")
                sort_database(f"{db_basename}_{suffix}.ffdata", f"{db_basename}_{suffix}.ffindex")
                optimize_database(f"{db_basename}_{suffix}.ffdata", f"{db_basename}_{suffix}.ffindex")
        finally:
            shutil.rmtree(tmp_dir)
    else:
        sys.stderr.write("You may try to use the option --force to fix the database!\n")
    
##########################################################################

def check_database(output_basename, threads, force_mode):
    """Check for the fasta database and the metadata of the fasta database .

    Args:
        output_basename (string): Path to the database folder with the prefix
        threads (int): number of cpus
        force_mode (bool): if value need to be forced or not
    """
    if not (os.path.exists(f"{output_basename}_a3m.ffindex") and os.path.exists(f"{output_basename}_a3m.ffdata")):
        sys.stderr.write("Error: No a3m database found!")
        exit(1)
        
    if not (os.path.exists(f"{output_basename}_cs219.ffindex") and os.path.exists(f"{output_basename}_cs219.ffdata")):
        open(f"{output_basename}_cs219.ffindex", 'a').close()
        open(f"{output_basename}_cs219.ffdata", 'a').close()
        
    if not (os.path.exists(f"{output_basename}_hhm.ffindex") and os.path.exists(f"{output_basename}_hhm.ffdata")):
        calculate_hhm(threads, f"{output_basename}_a3m", f"{output_basename}_hhm")

    if HHM_FOCUS == 'hhm':
        new_content = ''
        with open(f"{output_basename}_a3m.ffindex", 'rt') as r_file:
            for line in r_file:
                new_content += line.replace('a3m', 'hhm')
        with open(f"{output_basename}_a3m.ffindex", 'wt') as w_file:
            w_file.write(new_content)

    check_a3m_format(output_basename, force_mode)
    #TODO check hhm's...
    #TODO check cs219...
    
    handle_unsorted("a3m", output_basename, force_mode)
    handle_unsorted("hhm", output_basename, force_mode)
    handle_unsorted("cs219", output_basename, force_mode)
    
    handle_duplicates("a3m", calculate_hhm, threads, output_basename, force_mode)
    handle_duplicates("hhm", calculate_hhm, threads, output_basename, force_mode)
    handle_duplicates("cs219", calculate_cs219, threads, output_basename, force_mode)

    handle_missing("cs219", calculate_cs219, output_basename, threads, force_mode)
    
    handle_overhead("hhm", output_basename, force_mode)
    handle_overhead("cs219", output_basename, force_mode)

##########################################################################

def add_new_files(globular_expression, suffix, output_basename):
    """Add new files to the database .

    Args:
        globular_expression (string): globular expression of the files
        suffix (string): trype of data
        output_basename (string): path to database file
    """

    tmp_dir = tempfile.mkdtemp()

    files = set(glob(globular_expression))
    
    if(len(files) == 0):
        return

    try:
        files_index = os.path.join(tmp_dir, "files.dat")
        write_set_to_file(files, files_index)
        
        new_base = os.path.join(tmp_dir, "new")
        new_index_file =  f"{new_base}.ffindex"
        new_data_file = f"{new_base}.ffdata"
        build_ffindex_database(files_index, new_data_file, new_index_file)
    
        output_index_file = f"{output_basename}_{suffix}.ffindex"
        output_data_file = f"{output_basename}_{suffix}.ffdata"
        remove_files_from_index(files_index, output_index_file)
        merge_databases(new_data_file, new_index_file, output_data_file, output_index_file)
        optimize_database(output_data_file, output_index_file)
        
        if(suffix == "a3m"):
            for other_suffix in ["hhm", "cs219"]:
                output_index_file = f"{output_basename}_{other_suffix}.ffindex"
                output_data_file = f"{output_basename}_{other_suffix}.ffdata"
                if os.path.exists(output_index_file) and os.path.exists(output_data_file) and os.path.getsize(output_data_file) > 0:
                    remove_files_from_index(files_index, output_index_file)
                    optimize_database(output_data_file, output_index_file)
        
    finally:
        shutil.rmtree(tmp_dir)

##########################################################################

def main():
  a3m_files = snakemake.params.glob_a3ms
  hhm_files = snakemake.params.glob_hhms

  #Important to do a3m's first... deleting out of date hhm's and cs219's
  add_new_files(a3m_files, "a3m", snakemake.params.database_name)
    
  add_new_files(hhm_files, "hhm", snakemake.params.database_name)

  check_database(snakemake.params.database_name, snakemake.threads, True)
  
##########################################################################

if __name__ == "__main__":
  main()
  
