from optparse import OptionParser
import time
import numpy as np
from Observations import *
from Merge_Files import *

parser = OptionParser()
parser.add_option("-L", "--combi_list", type="string", default="None", help="total list of combis [%default]")

opts, args = parser.parse_args()


combi_file= open(opts.combi_list,'rb')
r=[]
for line in combi_file.readlines():
    res=[]
    for val in line.strip().split(' '):
        res.append(int(val))
    r.append(tuple(res))

pos_tab=np.rec.fromrecords(r,names=('fielda','fieldb','fieldc'))
thedir='../Ana_Cadence/OpSimLogs/WFD'
for val in pos_tab:
    fielda=val['fielda']
    fieldb=val['fieldb']
    fieldc=val['fieldc']
    obsa=Observations(fielda,filename=thedir+'/Observations_WFD_'+str(fielda)+'.txt')
    obsb=Observations(fieldb,filename=thedir+'/Observations_WFD_'+str(fieldb)+'.txt')
    obsc=Observations(fieldc,filename=thedir+'/Observations_WFD_'+str(fieldc)+'.txt')
    timeref=time.time()
    print 'Processing',fielda,fieldb,fieldc
    Merge_Fields([obsa,obsb,obsc])
    print 'after',time.time()-timeref
