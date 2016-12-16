from __future__ import print_function,division
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import ipdb
fn_table = '47tuc_orig.csv'

table = pd.read_csv(fn_table,sep=',',comment='#')
table.rename(columns= lambda x: x.strip(),inplace=True)

names =  [table['DATASET'][ii].strip() for ii in range(len(table)) if ii%6==0]
coords = [SkyCoord(np.float(table['WFC-MEUR'][ii])*360./24.,
                   np.float(table['HRC-MEUR'][ii]),
                   frame='icrs',unit='deg')
          for ii in range(len(table)) if ii%6==0]
#dec=     [table['HRC-MEUR'][ii] for ii in range(len(table)) if ii%6==0]
pma =    [np.float(table['MEYLANe1'][ii])*1000. for ii in range(len(table)) if ii%6==5]
pmd =    [np.float(table['GILLILU2'][ii])*1000. for ii in range(len(table)) if ii%6==5]
magv =   [np.float(table['GILLILU1'][ii]) for ii in range(len(table)) if ii%6==0]
magu =   [np.float(table['MEYLANe3'][ii]) for ii in range(len(table)) if ii%6==0]
magf475w=[np.float(table['GILLILU2'][ii]) for ii in range(len(table)) if ii%6==0]

ra = []
dec= []
for ii in range(len(coords)):
    ra.append(str(int(coords[ii].ra.hms[0]))+' '+\
              str(int(coords[ii].ra.hms[1]))+' '+\
              str(coords[ii].ra.hms[2]))
    dec.append(str(int(coords[ii].dec.dms[0]))+' '+\
               str(int(abs(coords[ii].dec.dms[1])))+' '+\
              str(abs(coords[ii].dec.dms[2])))

table_new = pd.DataFrame({'MAIN_ID':names,
                          'RA':ra,
                          'DEC':dec,
                          'PMRA':pma,
                          'PMDEC':pmd,
                          'MAG_V':magv,
                          'MAG_U':magu,
                          'MAG_USE':magf475w
                          })
table_bright  = table_new.loc[table_new['MAG_USE']<=165 ]

table_bright.to_csv('47tuc.csv',sep=',',columns=['MAIN_ID','RA','DEC','PMRA','PMDEC','MAG_V','MAG_U','MAG_USE'],index=False)

print('DONE.')
ipdb.set_trace()


