from astropy.wcs import WCS
import mastcasjobs
import pandas as pd

class CASjobs_sources(object):
    
    def __init__(self,info,maglim=20,band='i',context='ps1', name = None, path=None):
        if type(info) == str:
            self.wcs = WCS(info)
        elif type(info) == astropy.wcs.wcs.WCS:
            self.wcs = info
        elif type(info) == astropy.io.fits.header.Header:
            self.wcs = WCS(info)
        self.ra = None
        self.dec = None
        self.rad = None
        self.context = context
        self.name = name
        self.maglim = maglim
        self.band = band
        self.table = None
        self.query = None
        if path is None:
            self.path = './'
        else:
            self.path = path
        
    def get_coords(self):
        image  = np.zeros(self.wcs.array_shape)
        centre = [image.shape[1]//2,image.shape[0]//2]
        self.ra, self.dec = wcs.all_pix2world(centre[0],centre[1],0)
        
        size = np.sqrt((image.shape[1]-centre[0])**2 + (image.shape[0]-centre[1])**2) + 5
        pix_size = np.max(abs(self.wcs.pixel_scale_matrix))
        self.rad = size * pix_size * 60 # size in arc minutes
        return 
    
    def _check_params(self):
        m = ''
        message = "\n {} not defined."
        if self.name is None:
            print(message.format(name) + 'assigning default name.')
            self.name = 'default'
        if self.ra is None:
            m += message.format(ra)
        if self.dec is None:
            m += message.format(dec)
        if self.rad is None:
            m += message.format(rad)
        if len(m) > 2:
            raise ValueError(m)
        return
    
    def get_query(self):
        
        self._check_params()
        
        if self.context.lower() == 'ps1':
            ps1_query = """
                        select 
                        o.raMean, o.decMean,o.gMeanPSFMag,o.rMeanPSFMag,
                        o.iMeanPSFMag,o.zMeanPSFMag,o.yMeanPSFMag,o.iMeanKronMag
                        into mydb.[{dbname}]
                        from fGetNearbyObjEq({ra},{dec},{rad}) x
                        JOIN MeanObjectView o on o.ObjID=x.ObjId
                        LEFT JOIN StackObjectAttributes AS soa ON soa.objID = x.objID
                        WHERE o.nDetections > 5
                        AND soa.primaryDetection > 0
                        AND o.{band}MeanPSFMag < {maglim}
                        """
            self.name = 'ps1_'+self.name
            self.query = ps1_query.format(dbname=self.name,ra=self.ra,
                                          dec=self.dec,rad=self.rad,
                                          band=self.band,maglim=self.maglim)
        if self.context.lower() == 'gaia':
            gaia_query = """
                         select T.ra,T.dec,T.phot_g_mean_mag as gaia
                         into mydb.[{dbname}]
                         from fGetNearbyObjEq({ra},{dec},{rad}) as n
                         join GAIApublicVOview T on T.source_id = n.objID
                         WHERE T.phot_g_mean_mag < {maglim}
                         """
            self.name = 'gaia_'+self.name
            self.query = ps1_query.format(dbname='gaia_'+self.name,ra=self.ra,
                              dec=self.dec,rad=self.rad,
                              band=self.band,maglim=self.maglim)
        return 
    
    def submit_query(self,reset= True):
        if self.context == 'ps1':
            c = 'PanSTARRS_DR2'

        elif self.context == 'gaia':
            c = 'GAIA_DR2'
            
        else:
            raise ValueError('Only gaia and ps1 available now.')


        jobs = mastcasjobs.MastCasJobs(context=c)
        if reset:
            jobs.drop_table_if_exists(self.name)

        job_id = jobs.submit(self.query)
        status = jobs.monitor(job_id)
        print(status)
        self.table = jobs.get_table('tester_footprint',format='CSV').to_pandas()
        return 

    def save_space(self):
        """
        Creates a pathm if it doesn't already exist.
        """
        try:
            if not os.path.exists(self.path):
                os.makedirs(self.path)
        except FileExistsError:
            pass
        return

    def save_table(self):
        self.save_space()
        self.table.to_csv(self.path+self.name + '.csv', index = False)

    def get_table(self,save=None):
        self.get_coords()
        self.get_query()
        self.submit_query()
        if save is not None:
            self.save_table()
        return







