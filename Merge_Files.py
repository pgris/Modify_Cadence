import numpy as np
from Observations import *
import pylab as plt
from Sky_Brightness_With_Moonlight import SkyBrightness
import math
import palpy as pal
from Parameters import parameters
from lsst.sims.photUtils import SignalToNoise

DEG2RAD = math.pi / 180.    # radians = degrees * DEG2RAD
RAD2DEG = 180. / math.pi    # degrees = radians * RAD2DEG
TWOPI = 2 * math.pi
DAY = 86400.
simEpoch=59580

class Merge_Fields:
    def __init__(self,fields=[],merge_factor=0.8):

        self.merge_factor=merge_factor
        self.output_dir='OpSimLogs/WFD_Rolling'
        self.prefix='Observations_WFD_'

        self.mjdCol='mjd'
        self.filterCol='band'
        bands_orig=['u','g','r','i','z','y']
        self.bands=[]
        for band in bands_orig:
            self.bands.append('LSSTPG::'+band)


        fields_new=[]
        season_length=[]
        for field in fields:
            fields_new.append(Observations(field.fieldid,field.Ra,field.Dec,nseasons=len(field.seasons)))
            #print 'new field',field.fieldid,'Nseasons=',len(field.seasons)
            season_length.append(len(field.seasons))

        min_season=np.min(season_length)

        #season 0: same fields

        for i in range(len(fields_new)):
            fields_new[i].seasons[0]=fields[i].seasons[0]

        #Loop on seasons and make the merging

        combi_orig=[i for i in range(len(fields))]

        combi=list(combi_orig)

        #print seasons_to_tag
        for j in range(min_season-1):
            #print combi
            self.Concat(combi,j+1,fields,fields_new)
            combi=combi[1:]+combi[:1]

                #break
            #break
        #remaining periods
        
        for idx in combi_orig:
            for io in range(len(fields[idx].seasons)):
                #print 'completing here',idx,io,fields_new[idx].fieldid
                if fields_new[idx].seasons[io] is None and fields[idx].seasons[io] is not None:
                    #print 'passed'
                    fields_new[idx].seasons[io]=fields[idx].seasons[io]
        
        #print 'hello',fielda
        #self.Split_field(fielda,0,0.2)
        


        Plot=False
        if Plot:
            for idx,field in enumerate(fields):
                print field.fieldid,len(field.seasons),len(fields_new[idx].seasons)
                self.Plot(fields,fields_new,field.fieldid)
            """
            self.Plot(fields,fields_new,1716)
            self.Plot(fields,fields_new,2658)
            """
            #self.Plot_Filters(fields,fields_new,1074,valx='sky',legx='msky')
            plt.show()
        
        Write=True

        if Write:
            vars=['band','mjd','exptime','rawSeeing','seeing','moon_frac','sky','kAtm','airmass','m5sigmadepth','Nexp','Ra','Dec'] 
            bandeau='# band : \n# mjd : \n# exptime : \n# rawSeeing : \n# seeing : [was FWHMeff] \n# moon_frac : \n# sky : \n# kAtm : \n# airmass : \n# m5sigmadepth : \n# Nexp : \n# Ra : \n# Dec : \n# end\n'

            """
            for field in fields:
                print 'field',field.fieldid
                for key, val in field.seasons.items():
                    print key,len(val)
            print 'new fields'
            """
            for field in fields_new:
                #print 'field',field.fieldid
                outfile = open(self.output_dir+'/'+self.prefix+str(field.fieldid)+'.txt','wb')
                outfile.write(bandeau)
                recap=None
                for key, val in field.seasons.items():
                    #print key,len(val)
                    if recap is None:
                        recap=val
                    else:
                        recap=np.append(recap,val)
                
                for val in recap:
                    tot=''
                    for var in vars:
                        tot+=str(val[var])+' '
                    outfile.write(tot+' \n')
                #print recap
                #recap.dump(outfile)
                outfile.close()



    def Print_Nobs(self,fields):
        print "Number of observations per season and per band"
        for field in fields:
            print 'Field',field.fieldid
            for key, val in field.seasons.items():
                res=''
                for band in self.bands:
                    if val is not None:
                        sel=val[np.where(val[self.filterCol]==band)]
                        res+=band[-1]+': '+str(len(sel))+' '
                if res!='':
                    print key, res

    def Split_Field(self,season,frac=1.):
        
        if frac > 0.99:
            return season, None
        
        season.sort(order='mjd')

        n_tobe_excluded=(1.-frac)*len(season)
        if n_tobe_excluded < 1 and  n_tobe_excluded > 0.:
            n_tobe_excluded=1
        
        
        excluded=[]
        while len(excluded)< np.rint(n_tobe_excluded):
            aleat=np.random.randint(0,len(season))
            if aleat not in excluded:
                excluded.append(aleat)

        remain = list(set([x for x in range(len(season))]) - set(excluded))
        """
        print excluded,remain
        print len(excluded),len(remain)
        print 'Remaining'
        print season[remain]
        print 'excluded'
        print season[excluded]
        """
        return season[remain],season[excluded]
        
    def Concat(self,combi,iseason,fields_input,fields_rolling):

        percent=[self.merge_factor]*len(combi)
        percent[0]=1
        
        shift=[]

        #print fields_input[combi[0]].seasons
        ref_season=fields_input[combi[0]].seasons[iseason]

        for val in combi:
            shift.append(self.Get_Shift(ref_season,fields_input[val].seasons[iseason]))

        #print 'hello',shift 
        
        season_keep_all=[]
        season_rem_all=[]

        
        for j,val in enumerate(combi):
            season_keep=None
            season_remain=None
            for band in self.bands:
                season=fields_input[val].seasons[iseason]
                season_band=season[np.where(season[self.filterCol]==band)]
                
                season_k,season_rem=self.Split_Field(season_band,percent[j])
                
                if season_keep is None:
                    season_keep=season_k
                else:
                    season_keep=np.append(season_keep,season_k)
                
                if season_remain is None:
                    if season_rem is not None:
                        season_remain=season_rem
                else:
                    season_remain=np.append(season_remain,season_rem)
            
            #print 'resultat',len(season_keep)
            season_keep_all.append(season_keep)
            season_rem_all.append(season_remain)

            #print 'alors',len(season_keep_all[0]),
        fields_rolling[combi[0]].seasons[iseason]=season_keep_all[0]
        ra=fields_rolling[combi[0]].Ra
        dec=fields_rolling[combi[0]].Dec
        for j in range(1,len(season_keep_all)):
            fields_rolling[combi[0]].seasons[iseason]=np.append(fields_rolling[combi[0]].seasons[iseason],self.Shift(season_keep_all[j],shift[j],ra,dec))
            fields_rolling[combi[j]].seasons[iseason]=season_rem_all[j]

        #and do a Reshuffling of the master field

        #self.Print_Nobs(fields_rolling)
        fields_rolling[combi[0]].seasons[iseason]=self.Reshuffle(fields_rolling[combi[0]].seasons[iseason],fields_input[combi[0]].seasons[iseason])

    def Get_Shift(self,array_ref,array_add):

        mean_ref=np.mean(array_ref[self.mjdCol])
        mean_add=np.mean(array_add[self.mjdCol])

        resu = int(mean_ref - mean_add)
            
        return resu

    def Shift(self,season,shift,ra,dec):

        season[self.mjdCol][:]+=shift
        season['Ra'][:]=ra
        season['Dec'][:]=dec
        return season

    def Plot(self, fields, fields_new,fieldref,valx='mjd',valy='airmass',legx='MJD',legy='airmass',fontsize=12):

        #fieldref=1074
        figa, axa = plt.subplots(ncols=1, nrows=1, figsize=(6,5))
        for field in fields:
            if field.fieldid==fieldref:
               for key, val in field.seasons.items():
                   axa.plot(val[valx],val[valy],'bo')
                   print key,len(val)
               for field in fields_new:
                   if field.fieldid==fieldref:
                       for key, val in field.seasons.items():
                           if val is not None:
                               print key,len(val)
                               axa.plot(val[valx],val[valy],'r*')
        axa.set_xlabel(legx,{'fontsize': fontsize})
        axa.set_ylabel(legy,{'fontsize': fontsize})
        #plt.show()

    def Plot_Filters(self, fields, fields_new,fieldref,valx='mjd',legx='MJD',legy='Number of Entries',fontsize=12):

        field=None
        field_new=None
        
        for fielda in fields:
            if fielda.fieldid==fieldref:
                field=fielda
                break
        for fielda in fields_new:
            if fielda.fieldid==fieldref:
                field_new=fielda
                break


        figb, axb = plt.subplots(ncols=2, nrows=3, figsize=(10,9))

        for j,band in enumerate(self.bands):
            if j<2:
                k=0
            if j>= 2 and j < 4:
                k=1
            if j>=4:
                k=2
                
            arr_a=None
            arr_b=None
            for key, val in field.seasons.items():
                if arr_a is None:
                    arr_a=val[np.where(val[self.filterCol]==band)]
                else:
                    arr_a=np.append(arr_a,val[np.where(val[self.filterCol]==band)])
                 
            for key, val in field_new.seasons.items():
                if val is not None:
                    print key,len(val) 
                    if arr_b is None:
                        arr_b=val[np.where(val[self.filterCol]==band)]
                    else:
                        arr_b=np.append(arr_b,val[np.where(val[self.filterCol]==band)])

            axb[k][j%2].hist(arr_a[valx],bins=15,histtype='step',fill=False,color='r')
            axb[k][j%2].hist(arr_b[valx],histtype='step',fill=False,color='k')
            axb[k][j%2].set_xlabel(legx,{'fontsize': fontsize})
            axb[k][j%2].set_ylabel(legy,{'fontsize': fontsize})
        #plt.show()    



    def Reshuffle(self,season,season_ref):
        
        #for band in ['LSSTPG::r']:
        season_reshuffled=None
        for band in self.bands:
            seas_band=season[np.where(season[self.filterCol]==band)]
        
            #print 'size',len(seas_band)
            periods=None
            if len(seas_band) > 1:
                mjd_min=np.min(seas_band[self.mjdCol])
                mjd_max=np.max(seas_band[self.mjdCol])
                ra=np.unique(seas_band['Ra'])[0]
                dec=np.unique(seas_band['Dec'])[0]
                if mjd_max-mjd_min > 1.:
                    periods=self.Get_Periods(mjd_min,mjd_max,ra,dec,band[-1],len(seas_band))
            if len(seas_band) > 0:    
                if season_reshuffled is None:
                    season_reshuffled=self.Reshuffle_band(seas_band,periods)
                else:
                    season_reshuffled=np.append(season_reshuffled,self.Reshuffle_band(seas_band,periods))

        return season_reshuffled

    def Reshuffle_band(self,season,periods):

        if periods is None:
            return season

        mean_per_band=len(season)
        if len(periods) > 1:
            mean_per_band=len(season)/(len(periods)-1)

        iobs=-1
        airmass=1.2
        ra=np.unique(season['Ra'])[0]
        dec=np.unique(season['Dec'])[0]
        filtre=np.unique(season[self.filterCol])[0][-1]
        
        season_reshuffled=season.copy()
        toremove=[]
        for key, vals in periods.items():
            #print key,len(vals)
            for mjd in vals:
                iobs+=1
                visitexptime=season_reshuffled['exptime'][iobs]
                seeing=season_reshuffled['rawSeeing'][iobs]
                for delta in [0.,1.,2.,-1.,-2.]:
                    res, mjd_res, airmass,filtskybrightness,moonfrac =self.Get_mjd(mjd+delta-0.5,mjd+delta+0.5,ra,dec,filtre)
                    if res: 
                        m5=self.Recompute_fiveSigmaDepth(airmass=airmass,visitexptime=visitexptime,filtskybrightness=filtskybrightness,filtre=filtre,seeing=seeing)
                        season_reshuffled['mjd'][iobs]=mjd_res
                        season_reshuffled['airmass'][iobs]=airmass
                        season_reshuffled['sky'][iobs]=filtskybrightness
                        season_reshuffled['m5sigmadepth'][iobs]=m5
                        season_reshuffled['moon_frac'][iobs]=moonfrac
                        break
                if not res:
                    print 'could not find any',airmass
                    toremove.append(iobs)

        for val in toremove:            
            np.delete(season_reshuffled,val)
           

        return season_reshuffled

    def Get_mjd(self, mjdmin, mjdmax, ra, dec,filtre):
        
        res=False
        mjd_rand=-1

        for tirage in range(50):
            mjd_rand=np.random.uniform(mjdmin,mjdmax)
            airmass_rand= self.Get_airmass(mjd_rand,ra,dec)
            dateCol_rand=0
            #print 'hello',airmass_rand
            mysky=SkyBrightness(ra,dec,mjd_rand,dateCol_rand,filtre,airmass_rand)
            if mysky.Twilight == False and airmass_rand<1.5 and airmass_rand >=1. and mysky.moonZD_DEG > 85. and mysky.moonPhase < 60:
                return True, mjd_rand,airmass_rand,mysky.new_skybrightness(),mysky.moonPhase
        return res, mjd_rand,airmass_rand,mysky.new_skybrightness(),mysky.moonPhase


    def Get_Periods(self,mjd_min,mjd_max,ra,dec,filtre, nobs):

        periods={}
        r=[]
        airmass=1.3
        for mjd in np.arange(mjd_min,mjd_max,1./48.):
            mysky=SkyBrightness(ra,dec,mjd,0.,filtre,airmass)
             #if mysky.moonPhase < 60.:
            if mysky.Twilight == False and mysky.target_airmass<1.5 and mysky.moonZD_DEG > 85. and mysky.moonPhase < 60:
                r.append((mjd,mysky.moonPhase,mysky.distance2moon_DEG,mysky.new_skybrightness(),mysky.moonZD_DEG,mysky.target_airmass))
        #print mjd_min,mjd_max,mjd_max-mjd_min,len(r)

        if len(r) == 0:
            return None
        tabc=np.rec.fromrecords(r,names=['mjd','moonPhase','dist2Moon','sky','moonZD_DEG','airmass'])

        thediff=tabc['mjd'][1:]-tabc['mjd'][:-1]
        idx,=np.where(thediff > 10.)
        lidx=[val+1 for val in list(idx)]
        lidx.insert(0,0)
        lidx.append(len(tabc['mjd'])-1)
        #print lidx
        iper=-1
        for val in range(len(lidx)-1):
            iper+=1
            #print iper,iper+1,len(lidx)
            periods[iper]=[tabc['mjd'][lidx[iper]],tabc['mjd'][lidx[iper+1]-1]]
        

        n_per_period=nobs/len(periods)

        n_last=nobs-n_per_period*(len(periods)-1)

        npper=[n_per_period]*(len(periods)-1)
        npper.append(n_last)

        #print npper

        periods_compl={}
        for key, val in periods.items():
            periods_compl[key]=[va for va in np.linspace(val[0],val[1],npper[key])]


        return periods_compl
             

    def Reshuffle_old(self,season,season_ref):

        #for band in self.bands:

        vala=[]
        valb=[]
        for band in ['LSSTPG::r']:
            seas_band=season[np.where(season[self.filterCol]==band)]
            
            seas_band_ref=season_ref[np.where(season_ref[self.filterCol]==band)]
            print 'size',len(seas_band),len(seas_band_ref)
            mjd_min=np.min(seas_band[self.mjdCol])
            mjd_max=np.max(seas_band[self.mjdCol])
            for obs in seas_band:
                ra=obs['Ra']
                dec=obs['Dec']
                mjd=obs[self.mjdCol]
                filtre=obs[self.filterCol][-1]
                airmass=obs['airmass']
                mysky=SkyBrightness(ra,dec,mjd,0.,filtre,airmass)
                print filtre,mysky.distance2moon_DEG,mysky.moonPhase,mysky.new_skybrightness(),ra,dec
                vala.append((mjd,mysky.moonPhase,mysky.distance2moon_DEG,mysky.new_skybrightness()))

            for obs in seas_band_ref:
                ra=obs['Ra']
                dec=obs['Dec']
                mjd=obs[self.mjdCol]
                filtre=obs[self.filterCol][-1]
                airmass=obs['airmass']
                mysky=SkyBrightness(ra,dec,mjd,0.,filtre,airmass)
                #print filtre,mysky.distance2moon_DEG,mysky.moonPhase,mysky.new_skybrightness(),mysky.new_skybrightness()-obs['sky'],airmass-mysky.target_airmass
                valb.append((mjd,mysky.moonPhase,mysky.distance2moon_DEG,mysky.new_skybrightness()))

        taba=np.rec.fromrecords(vala,names=['mjd','moonPhase','dist2Moon','sky'])
        tabb=np.rec.fromrecords(valb,names=['mjd','moonPhase','dist2Moon','sky'])


        figa, axa = plt.subplots(ncols=2, nrows=2, figsize=(15,9))
        
        self.Plot_Indiv(axa[0][0],taba,valx='mjd',valy='moonPhase',mcol='bo',legx='MJD',legy='MoonPhase')
        self.Plot_Indiv(axa[0][0],tabb,valx='mjd',valy='moonPhase',mcol='r*')

        self.Plot_Indiv(axa[0][1],taba,valx='mjd',valy='dist2Moon',mcol='bo',legx='MJD',legy='dist2Moon [DEG]')
        self.Plot_Indiv(axa[0][1],tabb,valx='mjd',valy='dist2Moon',mcol='r*')

        self.Plot_Indiv(axa[1][0],taba,valx='sky',valy='moonPhase',mcol='bo',legx='msky',legy='MoonPhase')
        self.Plot_Indiv(axa[1][0],tabb,valx='sky',valy='moonPhase',mcol='r*')

        self.Plot_Indiv(axa[1][1],taba,valx='sky',valy='dist2Moon',mcol='bo',legx='msky',legy='dist2Moon [DEG]')
        self.Plot_Indiv(axa[1][1],tabb,valx='sky',valy='dist2Moon',mcol='r*')
        

        print 'hello',mjd_min,mjd_max

        

        ra=np.unique(season['Ra'])[0]
        dec=np.unique(season['Dec'])[0]
        filtre='r'
        airmass=1.3
        r=[]
        for mjd in np.arange(mjd_min,mjd_max,1./24.):
             mysky=SkyBrightness(ra,dec,mjd,0.,filtre,airmass)
             #if mysky.moonPhase < 60.:
             if mysky.Twilight == False and mysky.target_airmass<1.5 and mysky.moonZD_DEG > 85. and mysky.moonPhase < 60:
                 r.append((mjd,mysky.moonPhase,mysky.distance2moon_DEG,mysky.new_skybrightness(),mysky.moonZD_DEG,mysky.target_airmass))
             

        tabc=np.rec.fromrecords(r,names=['mjd','moonPhase','dist2Moon','sky','moonZD_DEG','airmass'])

        figb, axb = plt.subplots(ncols=2, nrows=2, figsize=(15,9))
        
        self.Plot_Indiv(axb[0][0],tabc,valx='mjd',valy='moonPhase',mcol='bo',legx='MJD',legy='MoonPhase')
        self.Plot_Indiv(axb[0][1],tabc,valx='moonPhase',valy='moonZD_DEG',mcol='bo',legx='MoonPhase',legy='moonZD_DEG')
        self.Plot_Indiv(axb[1][0],tabc,valx='sky',valy='moonPhase',mcol='bo',legx='msky',legy='MoonPhase')
        self.Plot_Indiv(axb[1][1],tabc,valx='sky',valy='moonZD_DEG',mcol='bo',legx='msky',legy='moonZD_DEG')
        
        thediff=tabc['mjd'][1:]-tabc['mjd'][:-1]
        idx,=np.where(thediff > 10.)
        lidx=[val+1 for val in list(idx)]
        lidx.insert(0,0)
        lidx.append(len(tabc['mjd'])-1)
        print lidx,tabc['mjd'][lidx]
        
        plt.show()

    def Plot_Indiv(self,ax,tab,valx,valy,mcol='bo',legx='',legy='',fontsize=12.):
        ax.plot(tab[valx],tab[valy],mcol)

        if legx!='':
            ax.set_xlabel(legx,{'fontsize': fontsize})
        if legy!='':
            ax.set_ylabel(legy,{'fontsize': fontsize})

    def Get_Sidereal_Time_at_Site(self,mjd):
        
        lon_RAD=-70.7494 *DEG2RAD #obs site long (from sims_operations/conf/system/SiteCP.conf)
        
        lst_RAD = pal.gmst(mjd) + lon_RAD
        if lst_RAD < 0:
            lst_RAD += TWOPI

        return lst_RAD

    def Get_airmass(self,mjd,ra_RAD,dec_RAD):
        
        lha_RAD = self.Get_Sidereal_Time_at_Site(mjd)-ra_RAD

        lat_RAD= -30.2444* DEG2RAD #obs site lat (from sims_operations/conf/system/SiteCP.conf)

        (az_RAD, d1, d2, alt_RAD, d4, d5, pa_RAD, d7, d8) = pal.altaz(lha_RAD,dec_RAD,lat_RAD)

        # Altitude -> Zenith distance (degrees)
        targetZD_DEG = 90. - (alt_RAD* RAD2DEG)
        # Altitude -> Zenith distance (radian)
        zd_RAD = 1.5707963 - alt_RAD
        # Airmass
        #am = slalib.sla_airmas (zd_RAD)
        am = pal.airmas(zd_RAD)
        return am

    def Recompute_fiveSigmaDepth(self,airmass=1.1,visitexptime=30.,filtskybrightness=21,filtre='r',seeing=0.05,FWHMeff=-1):
        
        param=parameters()
        if FWHMeff < 0:
            Filter_Wavelength_Correction = np.power(500.0 / param.filterWave[filtre], 0.3)
            Airmass_Correction = math.pow(airmass,0.6)
            FWHM_Sys = param.FWHM_Sys_Zenith * Airmass_Correction
            FWHM_Atm = seeing * Filter_Wavelength_Correction * Airmass_Correction
            finSeeing = param.scaleToNeff * math.sqrt(np.power(FWHM_Sys,2) + param.atmNeffFactor * np.power(FWHM_Atm,2))
            FWHMeff = SignalToNoise.FWHMgeom2FWHMeff(finSeeing)

        Tscale = visitexptime/ 30.0 * np.power(10.0, -0.4*(filtskybrightness - param.msky[filtre]))
        dCm = param.dCm_infinity[filtre] - 1.25*np.log10(1 + np.power(10.,0.8*param.dCm_infinity[filtre]- 1.)/Tscale)

        m5_recalc=dCm+param.Cm[filtre]+0.5*(filtskybrightness-21.)+2.5*np.log10(0.7/finSeeing)-param.kAtm[filtre]*(airmass-1.)+1.25*np.log10(visitexptime/30.)
        
        return m5_recalc
