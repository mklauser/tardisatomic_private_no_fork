#!/usr/bin/env python
__author__ = 'michi'



import tardisatomic

import argparse
import h5py
import os
import sys
import sqlite3
import re
from tardisatomic import import_chianti,  import_regemorter


try:
    import sqlparse
    sqlparse_available = True
except ImportError:
    sqlparse_available = False


basic_atom_data = h5py.File(os.path.join(os.path.dirname(tardisatomic.__file__), 'data', 'atom_data_basic.h5'))['basic_atom_data']
symbol2z = dict(zip(basic_atom_data['symbol'], basic_atom_data['atomic_number']))


parser = argparse.ArgumentParser()

parser.add_argument('input_db', help='HDF5 file that is created')
parser.add_argument('atomic_number', nargs='+')
parser.add_argument('-o', '--output_db', default=None, help='Output Database - if not set will overwrite input db')


args = parser.parse_args()

if args.output_db is None:
    print """WARNING OVERWRITING CURRENT DATABASE - please specify an output db if this is not wished"""
    continue_prompt = raw_input('Continue [y/N]')
    if not continue_prompt.lower().startswith('y'):
        sys.exit(0)
    conn = sqlite3.connect(args.input_db)
else:
    os.system('cp %s %s' % (args.input_db, args.output_db))
    conn = sqlite3.connect(args.output_db)

curs = conn.cursor()
if curs.execute('pragma table_info(collision_data)').fetchall() == []:
    print "collision_data table does not exist - creating it"
    import_chianti.create_collision_data_table(conn)

for atomic_number in args.atomic_number:
    atomic_number = int(atomic_number)
    print('computing the rates for %d'%atomic_number)
    import_regemorter.add_regemorter(conn, atomic_number)

conn.close()

