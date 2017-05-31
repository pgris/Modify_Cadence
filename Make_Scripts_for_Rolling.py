#!/usr/bin/python

import os
import time
import sys
from optparse import OptionParser
import glob
import numpy as np


def submit(combi_list):

    scriptdir="scripts"
    if not os.path.isdir(scriptdir):
        os.makedirs(scriptdir)
    
    procname="fake_rolling"

    procnum=len(glob.glob(scriptdir+'/*.sh'))
    procnum+=1
    procid=procname+"_"+str(procnum)
    filename = scriptdir+ "/" + procid+ ".list"
    theFile = open(filename,"w+")
    for val in combi_list:
        print>>theFile, val['fielda'],val['fieldb'],val['fieldc']
    

    cwd = os.getcwd()
    dirLog = cwd + "/logs"
    if not os.path.isdir(dirLog) :
        os.makedirs(dirLog)
    logname=procid+'.log'
    log = dirLog + "/" + logname
    
    cmd = 'python run_modify.py --combi_list '+filename +' >& '+log
    scriptName = "scripts/" + procid+ ".sh"
    script = open(scriptName,"w")
    script.write("#!/usr/local/bin/bash\n")
    script.write("cd " + cwd +"\n")
    #script.write("bash" + "\n")
    script.write("source setups.sh" + "\n")
    script.write(" " + cmd + "\n")
    script.close()
    os.system("chmod +x " + scriptName)
    #os.system("./"+scriptName)
    time.sleep(1)


parser = OptionParser()
parser.add_option("-N", "--mods", type="int", default=10, help="number of combi[%default]")
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

for spl in np.array_split(pos_tab,opts.mods):
    submit(spl)





