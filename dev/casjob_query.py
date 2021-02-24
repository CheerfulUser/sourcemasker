from astropy.wcs import WCS
from astropy.wcs.wcs import WCS as wcs_class
from astropy.io.fits.header import Header as header_class
import mastcasjobs
import pandas as pd
import numpy as np

class CASjobs_sources(object):
    """
    Query an entire image region for sources, on casjobs. Takes the image WCS to calculate
    image bounds and search area.

    Inputs
    ------
        info : wcs/str/header
            the wcs, file name, or header for the image to query the region of.
    Options
    -------
        maglim : float
            magnitude limit of the query 
        band : str
            ps1 band to do the magnitude cut 
        context : str
            casjobs query context, currently only ps1 and gaia are avaialble 
        name : str
            name of the database that will be used to save on casjobs and locally 
        path : str
            path to the save directory 
    """
    
    def __init__(self,info,maglim=20,band='i',context='ps1', name = None):
        if type(info) == str:
            self.wcs = WCS(info)
        elif type(info) == wcs_class:
            self.wcs = info
        elif type(info) == header_class:
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
        
    def get_coords(self):
        """
        Get the centre coordinates and query radius from the wcs 
        """
        dim1, dim2  = self.wcs.array_shape
        centre = [dim1//2,dim2//2]
        self.ra, self.dec = self.wcs.all_pix2world(centre[0],centre[1],0)
        
        size = np.sqrt((dim1-centre[0])**2 + (dim2-centre[1])**2) + 5
        pix_size = np.max(abs(self.wcs.pixel_scale_matrix))
        self.rad = size * pix_size * 60 # size in arc minutes
        return 
    
    def _check_params(self):
        """
        Check all relevant variables are defined
        """
        m = ''
        message = "\n {} not defined"
        if self.name is None:
            print(message.format('name') + ' assigning default name.')
            self.name = 'default'
        if self.ra is None:
            m += message.format('ra')
        if self.dec is None:
            m += message.format('dec')
        if self.rad is None:
            m += message.format('rad')
        if len(m) > 2:
            raise ValueError(m)
        return
    
    def get_query(self):
        """
        Get the query string to submit to casjobs
        """
        
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
            self.query = gaia_query.format(dbname=self.name,ra=self.ra,
                              dec=self.dec,rad=self.rad,
                              band=self.band,maglim=self.maglim)
        return 
    
    def submit_query(self,reset= True):
        """
        Submit the query and download the resulting table
        """
        if self.context == 'ps1':
            c = 'PanSTARRS_DR2'

        elif self.context == 'gaia':
            c = 'GAIA_DR2'
            
        else:
            raise ValueError('Only gaia and ps1 available now.')

        jobs = mastcasjobs.MastCasJobs(context=c)
        if reset:
            jobs.drop_table_if_exists(self.name)
        else:
            try:
                self.table = jobs.get_table(self.name,format='CSV').to_pandas()
                if self.context == 'ps1':
                    self.table = self.table.replace(-999,np.nan)
                print('loading existing table')
                return 
            except:
                pass

        job_id = jobs.submit(self.query)
        status = jobs.monitor(job_id)
        print(status)
        if status[0] != 5:
            raise ValueError('No table created')
        self.table = jobs.get_table(self.name,format='CSV').to_pandas()

        if self.context == 'ps1':
            self.table = self.table.replace(-999,np.nan)
        
        return 

    def save_space(self):
        """
        Creates a path if it doesn't already exist.
        """
        try:
            if not os.path.exists(self.path):
                os.makedirs(self.path)
        except FileExistsError:
            pass
        return

    def save_table(self,save):
        """
        Save the query output 
        """
        self.save_space()
        self.table.to_csv(save, index = False)

    def get_table(self,reset=True,save=None):
        """
        Runs all functions to get the table.
        """
        self.get_coords()
        self.get_query()
        self.submit_query(reset=reset)
        if save is not None:
            self.save_table(save)
        return







