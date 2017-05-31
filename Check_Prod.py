import numpy as np
import os
from Observations import *
import pylab as plt

bands=['u','g','r','i','z','y']

def cadence(obs):
    res={}

    for key, val in obs.seasons.items():
        cad={}
        for band in bands:
            sel=val[np.where(val['band']=='LSSTPG::'+band)]
            #print len(sel)
            if len(sel) > 2:
                diff=sel['mjd'][1:]-sel['mjd'][:-1]
                cad[band]=np.mean(diff)
            else:
                cad[band]=999.
        res[key]=cad

    return res

def Plot_Super(obs_old,obs_new,valx='mjd',valy='airmass',legx='MJD',legy='airmass'):

    figb, axb = plt.subplots(ncols=1, nrows=1, figsize=(15,9))
    
    for key,val in obs_old.seasons.items():
        axb.plot(val[valx],val[valy],'bo')
        
    for key,val in obs_roll.seasons.items():
        axb.plot(val[valx],val[valy],'r*')

    axb.set_xlabel(legx)
    axb.set_ylabel(legy)

def Plot_Indiv(ax,cad_old,cad_roll,band):

    per=cad_old.keys()
    #print per
    val_old=[]
    val_roll=[]
    
    for key,val in cad_old.items():
        val_old.append(val[band])
    ax[0].plot(per,val_old,'b+')
        
    for key,val in cad_roll.items():
        val_roll.append(val[band])
    print 'hello',per,val_roll
    ax[0].plot(per,val_roll[:len(per)],'r+')

    cadence_ratio=[x/y for x, y in zip(val_roll,val_old)]
    
    ax[1].plot(per,cadence_ratio,'k+')
    

def Plot_Cadence(cad_list_old,cad_list_new,band,field_type):

    figb, axb = plt.subplots(ncols=1, nrows=2, figsize=(15,9))
    
    for i in range(len(cad_list_old)):
        Plot_Indiv(axb,cad_list_old[i],cad_list_new[i],band)
    
    axb[0].set_xlabel('Season')
    axb[0].set_ylabel('Cadence$^{-1}$')
    axb[0].set_ylim(0, 20.)

    axb[1].set_xlabel('Season')
    axb[1].set_ylabel('Cadence$^{-1}$ ratio')

    axb[1].set_ylim(0,1)
    figb.suptitle('Field type:'+field_type+' - band '+band)

thelist='List_Combi.txt'

thedir='OpSimLogs/WFD_Rolling'
prefix='Observations_WFD_'
suffix='.txt'

thedir_orig='../Ana_Cadence/OpSimLogs/WFD'

sfile=open(thelist,'r')

obs={}

print 'Loading files'
r=[]
for line in sfile.readlines():
    rb=[]
    for field in line.strip().split(' '):
        fname=thedir+'/'+prefix+str(field)+suffix
        
        if not os.path.isfile(fname):
            print "Missing :",field
        else:
            rb.append(field)
    if len(rb) > 0:
        r.append(tuple(rb))

tab_fields=np.rec.fromrecords(r,names=('fielda','fieldb','fieldc'))


sel=tab_fields['fielda']

field=1811


cad_old=[]
cad_roll=[]

for val in range(100):
    field=sel[val]

    obs_old=Observations(field,filename=thedir_orig+'/'+prefix+str(field)+suffix)
    obs_roll=Observations(field,filename=thedir+'/'+prefix+str(field)+suffix)
    
    cado=cadence(obs_old)
    cadr=cadence(obs_roll)

    if len(cado) != len(cadr):
        print 'This is strange! ',field,len(cado),len(cadr),len(obs_old.seasons),len(obs_roll.seasons)

    cad_old.append(cado)
    cad_roll.append(cadr)

    """
    for key, val in cad_old.items():
        print key, val
        print key, cad_roll[key]
    """

#Plot_Super(obs_old,obs_roll)
Plot_Cadence(cad_old,cad_roll,'g','a')


plt.show()
