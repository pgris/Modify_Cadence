from Observations import *
from Merge_Files import *
import time
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table
import healpy as hp

def Get_Combis(pos_tab,ra_grid=2.,n_merger=3):
        
    ra_step=ra_grid # degrees
    n_dec_zone=n_merger # number of region in dec = number of fields to be merged

    pos_tab.sort(order='fieldRA')
    ra_min=np.min(pos_tab['fieldRA'])
    
    ra_dec_strip={}
    istrip=-1
    
    while ra_min < 360.-ra_step:
        istrip+=1
        ra_dec_strip[istrip]={}
        ra_max=ra_min+ra_step
        sel=pos_tab[np.where(np.logical_and(np.rad2deg(pos_tab['fieldRA'])>=ra_min,np.rad2deg(pos_tab['fieldRA'])<ra_max))]
        sel.sort(order='fieldDec')
        num_per_part=len(sel)/n_dec_zone
        ntag=0
        for count in range(n_dec_zone):
            if count == n_dec_zone-1:
                ra_dec_strip[istrip][count]=sel[ntag:]
            else:
                ra_dec_strip[istrip][count]=sel[ntag:ntag+num_per_part] 
            #print count, len(sel),num_per_part,len(ra_dec_strip[istrip][count])
            ntag+=num_per_part

        ra_min+=ra_step
             #break


    ra_dec_final_combi={}
    
    for iv in range(len(ra_dec_strip)):
        #print 'o yes',iv
        ra_dec_final_combi[iv]={}
        strip_copy={}
        icombi=-1
        for i in range(1,n_dec_zone):
            strip_copy[i]=ra_dec_strip[iv][i].copy()
        
        for val in ra_dec_strip[iv][0]:
            r=[]
            icombi+=1
            #restable=Table(names=('fieldID','fieldRA','fieldDec'), dtype=('i4', 'f8','f8'))
            r.append(tuple(val))
            for  iu in range(1,n_dec_zone):
                strip_copy[iu],resu=Get_Nearest(strip_copy[iu],val)
                r.append(tuple(resu))
            ra_dec_final_combi[iv][icombi]=np.rec.fromrecords(r,names=('fieldID','fieldRA','fieldDec'))

    all_combo=[]
    for key,vval in ra_dec_final_combi.items():
        for key,val in vval.items():
            local=[]
            for ik in range(len(val)):
                local.append(val['fieldID'][ik])
            all_combo.append(local)
         
    return all_combo
 
def Get_Nearest(orig_table, val):
   
    table=Table(orig_table)
    c = SkyCoord(ra=val['fieldRA']*u.radian, dec=val['fieldDec']*u.radian)  
    catalog = SkyCoord(ra=table['fieldRA']*u.radian, dec=table['fieldDec']*u.radian)  
    idx, d2d, d3d = c.match_to_catalog_sky(catalog)

    #print 'astropy matching',idx,d2d,d3d,len(table),type(idx),table
    theres=[table['fieldID'][int(idx)],table['fieldRA'][int(idx)],table['fieldDec'][int(idx)]]
    table.remove_row(int(idx))
    return table,theres

thedir='../Ana_Cadence/OpSimLogs/WFD'
thelist='List_WFD_test.txt'

sfile=open(thelist,'r')

obs={}

print 'Loading files'
r=[]
for line in sfile.readlines():
    spli=line.split(' ')
    iv=int(spli[0])
    obs[iv]=Observations(int(spli[0]),filename=thedir+'/'+spli[1].strip())
    r.append((iv,obs[iv].Ra,obs[iv].Dec))

#print r

#This is to get the combination of fields
"""
pos_tab=np.rec.fromrecords(r,names=('fieldID','fieldRA','fieldDec'))

combi_tot=Get_Combis(pos_tab)

for val in combi_tot:
    print val[0],val[1],val[2]
"""
sfile=open('List_Combi.txt','r')

obs={}

print 'Loading files'
r={}
for i in range(3):
    r[i]=[]

for line in sfile.readlines():
    spli=line.split(' ')
    for i,field in enumerate(spli): 
        obs=Observations(int(field),filename=thedir+'/Observations_WFD_'+field.strip()+'.txt')
        r[i].append((iv,obs.Ra,obs.Dec))

tabs={}
for i in range(3):
    tabs[i]=np.rec.fromrecords(r[i],names=('fieldID','fieldRA','fieldDec'))


fig = plt.figure(figsize=(8,6))

axa = fig.add_subplot(111, projection="mollweide")
colors=['k','r','b']
for i in range(3):
    axa.scatter(tabs[i]['fieldRA'],tabs[i]['fieldDec'],color=colors[i])

plt.show()


"""
combi_tot=[[1074,1873,2572]]

for combi in combi_tot:
    print 'before',time.time()
    timeref=time.time()
    Merge_Fields([obs[combi[0]],obs[combi[1]],obs[combi[2]]])
    print 'after',time.time()-timeref
"""            
