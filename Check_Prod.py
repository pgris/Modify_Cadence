import numpy as np
import os
from Observations import *
import pylab as plt
import astropy.coordinates as coord
import astropy.units as u

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

def Plot_Super(obs_old,obs_new,valx='mjd',valy='airmass',legx='MJD',legy='airmass',field=0):

    figb, axb = plt.subplots(ncols=1, nrows=1, figsize=(15,9))
    
    """
    y=[1.1,1.55]
    """
    for key,val in obs_old.seasons.items():
        axb.plot(val[valx],val[valy],'bo')
        print key,np.min(val[valx]),np.max(val[valx])
       

    for key,val in obs_roll.seasons.items():
        axb.plot(val[valx],val[valy],'r*')
        print key,np.min(val[valx]),np.max(val[valx])
        """
        minx=np.min(val[valx])
        maxx=np.max(val[valx])
        axb.plot([minx,minx],y,linestyle='-',color='k')
        axb.plot([maxx,maxx],y,linestyle='--',color='k')
        """
    axb.set_xlabel(legx)
    axb.set_ylabel(legy)
    figb.suptitle('Field '+str(field))
 
    plt.gcf().savefig('Plots/Airmass_vs_mjd_%s.png' % str(field), 
                         bbox_inches='tight')
    plt.gcf().savefig('Plots/Airmass_vs_mjd_%s.pdf' % str(field),
                         bbox_inches='tight')

def Plot_Indiv(ax,cad_old,cad_roll,band):

    per=cad_old.keys()
    per=[pp+1 for pp in per]
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

    figb, axb = plt.subplots(ncols=1, nrows=2, figsize=(12,12))
    
    for key,val in cad_list_old.items():
        Plot_Indiv(axb,val,cad_list_new[key],band)
    
    axb[0].set_xlabel('Season')
    axb[0].set_ylabel('Cadence$^{-1}$')
    axb[0].set_ylim(0, 20.)

    axb[1].set_xlabel('Season')
    axb[1].set_ylabel('Cadence$^{-1}$ ratio')

    axb[1].set_ylim(0,1)
    figb.suptitle('Field type:'+field_type+' - band '+band)

def Plot_Cadence_Hist(cad_list_old,cad_list_new,band,tab_fields,therange):

    values=Get_Values_for_Hist(cad_list_old,band,tab_fields)
    valuesb=Get_Values_for_Hist(cad_list_new,band,tab_fields)

    sela=[val for val in values if val < max(therange)]
    selb=[val for val in valuesb if val < max(therange)]
    print 'Band',band,'Wide : ',np.mean(sela),np.std(sela),'Rolling : ',np.mean(selb),np.std(selb)
   
    figb, axb = plt.subplots(ncols=1, nrows=1, figsize=(15,9))

    axb.hist(values,bins=2*max(therange),range=therange,histtype='step',color='k',linewidth=1.5, label='Wide Survey')
    axb.hist(valuesb,bins=2*max(therange),range=therange,histtype='step',color='r',linewidth=1.5,label='Rolling Cadence')
    axb.set_ylabel('Number of Entries', fontsize=18)
    axb.set_xlabel('Cadence$^{-1}$ [days]', fontsize=18)
    axb.set_xlim(therange[0],therange[1])
    axb.legend(loc='upper right',prop={'size':18})
    axb.tick_params(axis='x', labelsize=15)
    axb.tick_params(axis='y', labelsize=15)
    #figb.suptitle('Field type:'+field_type+' - band '+band)
    figb.suptitle(band)


    plt.gcf().savefig('Cadence_Hist_%s.png' % band, 
                         bbox_inches='tight')
    plt.gcf().savefig('Cadence_Hist_%s.pdf' % band,
                         bbox_inches='tight')

def Plot_msky_Hist(cad_list_old,cad_list_new,band,tab_fields,therange):

    values=Get_Values_median_for_Hist(cad_list_old,band,tab_fields,'sky')
    valuesb=Get_Values_median_for_Hist(cad_list_new,band,tab_fields,'sky')
   
    #print values,valuesb
    print 'Band',band,'Wide : ',np.median(values),'Rolling : ',np.median(valuesb),'Wide : ',np.mean(values),np.std(values),'Rolling : ',np.mean(valuesb),np.std(valuesb)

    figb, axb = plt.subplots(ncols=1, nrows=1, figsize=(15,9))

    #axb.hist(values,bins=2*max(therange),range=therange,histtype='step',color='k', label='Wide Survey')
    #axb.hist(valuesb,bins=2*max(therange),range=therange,histtype='step',color='r',label='Rolling Cadence')
    axb.hist(values,bins=20,histtype='step',color='k', label='Wide Survey')
    axb.hist(valuesb,bins=20,histtype='step',color='r',label='Rolling Cadence')
    axb.set_ylabel('Number of Entries')
    axb.set_xlabel('msky [mag]')
    #axb.set_xlim(therange[0],therange[1])
    axb.legend(loc='upper right')
    #figb.suptitle('Field type:'+field_type+' - band '+band)
    figb.suptitle(band)

    """
    plt.gcf().savefig('Cadence_Hist_%s.png' % band, 
                         bbox_inches='tight')
    plt.gcf().savefig('Cadence_Hist_%s.pdf' % band,
                         bbox_inches='tight')
    """

def Get_Values_for_Hist(cad_list,band,tab_fields):

    seasons={'fielda':[1,4,7],'fieldb':[2,5,8],'fieldc':[3,6,9]}
    values=[]

    for field_type in ['fielda','fieldb','fieldc']:
        sela={k:v for (k,v) in cad_list.items() if k in tab_fields[field_type]}
        
        #print sela.keys()

       
        for key,val in sela.items():
            for seas in seasons[field_type]:
            #print key,seas,val[seas]
                if val.has_key(seas):
                    if val[seas][band] < 990.:
                        values.append(val[seas][band])

    return values
    

def Get_Values_median_for_Hist(cad_list,band,tab_fields,what):

    seasons={'fielda':[1,4,7],'fieldb':[2,5,8],'fieldc':[3,6,9]}
    values=[]

    for field_type in ['fielda','fieldb','fieldc']:
        sela={k:v for (k,v) in cad_list.items() if k in tab_fields[field_type]}
        
        #print sela.keys()

        reca=None
        for key,val in sela.items():
            for seas in seasons[field_type]:
            #print key,seas,val[seas]
                #print 'hello',val.seasons[seas].dtype
                if val.seasons.has_key(seas):
                    index=val.seasons[seas]['band']=='LSSTPG::'+band
                    if reca is None:
                        reca=val.seasons[seas][index]
                    #print 'ici',reca
                    else:
                        reca=np.append(reca,val.seasons[seas][index])
            #print reca.dtype,val.seasons[seas][index]
            #print reca[what]
            values.append(np.median(reca[what]))
                

    return values


def Plot_Mollweid(tab_fields,tab_radec):
    
    #sel=tab_fields['fielda']
    #index = tab_radec['fieldid'] == [fi for fi in tab_fields['fielda']]
    #sel=np.select(tab_radec['fieldid'],tab_fields['fielda'])

    colors=['b','r','k']

    fig = plt.figure(figsize=(8,6))

    ax = fig.add_subplot(111, projection="mollweide")

    for i, tt in enumerate(['a','b','c']):
        index=np.in1d(tab_radec['fieldid'],tab_fields['field'+tt])
        sel=tab_radec[index]
 
        ra = coord.Angle(np.rad2deg(sel['Ra']),unit=u.degree)
        ra = ra.wrap_at(180*u.degree)
        dec = coord.Angle(np.rad2deg(sel['Dec']),unit=u.degree)
        ax.scatter(ra.radian,dec.radian,color=colors[i],marker='+')
    
    ax.grid(True) # afficher la grille

    plt.gcf().savefig('Rolling_fields.png', 
                         bbox_inches='tight')
    plt.gcf().savefig('Rolling_fields.pdf', 
                         bbox_inches='tight')

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

Plot_msky=False

if Plot_Cadence:

    obs_old={}
    obs_roll={}
    
#for val in range(100):
    for fields in tab_fields[:]:
        for field in fields:
        #field=sel[val]

            obs_old[field]=Observations(field,filename=thedir_orig+'/'+prefix+str(field)+suffix)
            obs_roll[field]=Observations(field,filename=thedir+'/'+prefix+str(field)+suffix)

    ranges={'u':[0,40],'g':[0,30],'r':[0,20],'i':[0,20],'z':[0,20],'y':[0,20]}

    for band in bands:
        Plot_msky_Hist(obs_old,obs_roll,band,tab_fields,ranges[band])




Plot_Cadence=True

if Plot_Cadence:

    cad_old={}
    cad_roll={}
    
#for val in range(100):
    for fields in tab_fields[:]:
        for field in fields:
        #field=sel[val]

            obs_old=Observations(field,filename=thedir_orig+'/'+prefix+str(field)+suffix)
            obs_roll=Observations(field,filename=thedir+'/'+prefix+str(field)+suffix)
            
            cado=cadence(obs_old)
            cadr=cadence(obs_roll)
        
            if len(cado) != len(cadr):
                print 'This is strange! ',field,len(cado),len(cadr),len(obs_old.seasons),len(obs_roll.seasons)

            cad_old[field]=cado
            cad_roll[field]=cadr

    """
    for key, val in cad_old.items():
        print key, val
        print key, cad_roll[key]
    """


    ranges={'u':[0,40],'g':[0,30],'r':[0,20],'i':[0,20],'z':[0,20],'y':[0,20]}

    for band in bands:
        Plot_Cadence_Hist(cad_old,cad_roll,band,tab_fields,ranges[band])

#for band in bands:
#    Plot_Cadence(cad_old,cad_roll,band,'a')


Plot_Sup=False

if Plot_Sup:
    field=511
    obs_old=Observations(field,filename=thedir_orig+'/'+prefix+str(field)+suffix)
    obs_roll=Observations(field,filename=thedir+'/'+prefix+str(field)+suffix)
    
    Plot_Super(obs_old,obs_roll,field=field)
#Plot_Cadence(cad_old,cad_roll,'g','a')

Plot_Mollweid=False

if Plot_Mollweid:

    r=[]
    for val in tab_fields[:]:
        for tt in ['a','b','c']:
            field=val['field'+tt]
            obs=Observations(field,filename=thedir_orig+'/'+prefix+str(field)+suffix)
            r.append((obs.fieldid,obs.Ra,obs.Dec))

    tab_radec=np.rec.fromrecords(r,names=('fieldid','Ra','Dec'))

    Plot_Mollweid(tab_fields,tab_radec)

plt.show()
