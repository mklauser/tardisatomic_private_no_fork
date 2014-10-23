__author__ = 'michi'

import numpy as np
import scipy.special as scisp
from astropy import table, constants, units

def add_regemorter(conn, atomic_number):
    curs = conn.cursor()
    curs.execute('SELECT T FROM T_grid')
    T = np.array(curs.fetchall())
    curs.execute('SELECT id, atom, ion, wl, f_lu, level_id_upper, level_id_lower, g_upper, g_lower  FROM lines WHERE atom=?', (str(atomic_number)))
    lines = curs.fetchall()
    for line in lines:
        atom = int(line[1])
        ion = int(line[2])
        nu = units.Unit('angstrom').to('Hz',line[3], units.spectral())
        f_lu = line[4]
        level_number_upper = int(line[5])
        level_number_lower = int(line[6])
        g_upper = line[7]
        g_lower = line[8]

        c_lu = compute_van_regemorter(T, f_lu, nu)
        C_ul_conversion = g_upper / float(g_lower)
        #Delete old
        delete_stmt = 'DELETE FROM collision_data WHERE atom=? AND ion=? AND level_number_lower = ? AND level_number_upper=?'
        curs.execute(delete_stmt,(str(atom),str(ion),str(level_number_lower), str(level_number_upper)))

        #Build T insert
        insert_stmt = 'insert into collision_data(atom, ion, level_number_lower, level_number_upper, %s, c_ul_conversion)' \
                      ' SELECT %s ' \
                      + 'WHERE NOT EXISTS(SELECT 1 FROM collision_data WHERE atom=%d AND ion=%d AND ' \
                        'level_number_lower=%d AND level_number_upper=%d)'%(atom, ion, level_number_lower, level_number_upper)


        c_lu = [float(i) for i in c_lu]
        collision_data = [atom, ion, level_number_lower, level_number_upper] + c_lu + [C_ul_conversion]
        values_str = ','.join([str(i) for i in collision_data])
        T_str = ', '.join(['t%06d' % item for item in T])
        insert_stmt = insert_stmt%(T_str,values_str)
        curs.execute(insert_stmt)

    curs.close()
    conn.commit()







def compute_van_regemorter(T, f_lu, nu_lu):
    g = 0.2 # This value is set to 2. We should select the value based on the main quantum number
    u = constants.h.cgs.value * nu_lu / constants.k_B.cgs.value / T
    I = 13.6 # eV
    c0 = 5.46510e-11
    integ = 0.276 * np.exp(u) * scisp.exp1(u)
    gamma = (g,integ )
    c= c0 * T**(0.5) *14.5 * (I / constants.h.cgs.value / nu_lu )* f_lu* constants.h.cgs.value * nu_lu \
       / constants.k_B.cgs.value / T * np.exp(- u)
    return c




