import numpy as np, healpy as hp
import os, sys
from tqdm import tqdm
import glob, h5py
import yaml
from sklearn.neighbors import BallTree
import fitsio
import joblib
from numpy.lib.recfunctions import stack_arrays

from constants import MEDSCONF
from simulating import End2EndSimulation

TMP_DIR = os.environ['TMPDIR']

class SuppressPrint:
    def __enter__(self):
        self._original_stdout = sys.stdout  # Save a reference to the original standard output
        sys.stdout = open(os.devnull, 'w')  # Redirect standard output to a null device

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()  # Close the null device
        sys.stdout = self._original_stdout  # Restore the original standard output

def match_catalogs(path, tilename, bands, config):

    #with SuppressPrint(): #To remove some annoying log10() and divide invalid errors
    #    SIM = End2EndSimulation.__new__(End2EndSimulation)
    #SIM.gal_kws  = config['gal_kws']
    #SIM.star_kws = config['star_kws']
    #SIM.source_rng = np.random.default_rng(seed = 42)
    
    #Input_catalog  = SIM._make_sim_catalog()

    '/project/chihway/dhayaa/DECADE/Balrog/v08_ProductionRun3/Input_DES1210+0043-cat.fits',
    '/project/chihway/dhayaa/DECADE/Balrog/v08_ProductionRun3/SrcExtractor_DES1003-3206_i-cat.fits',
    '/project/chihway/dhayaa/DECADE/Balrog/v08_ProductionRun3/OldSrcExtractor_DES1301-3540_i-cat.fits',
    '/project/chihway/dhayaa/DECADE/Balrog/v08_ProductionRun3/metacal_DES1451-0124.fits',

    
    Input_catalog = fitsio.read(r'%s/SplicedSim_Input_%s-cat.fits' % (path, tilename), ext = 1)
    
    #Get all paths
    sof_path   = r'%s/fitvd_%s.fits' % (path, tilename)
    Truth_path = r'%s/Input_%s-cat.fits' % (path, tilename)
    OCat_path  = [r'%s/OldSrcExtractor_%s_%s-cat.fits' % (path, tilename, band) for band in bands] #Path to original SrcExtractor
    BCat_path  = [r'%s/SrcExtractor_%s_%s-cat.fits' % (path, tilename, band) for band in bands] #Path to new SrcExtractor
    
    #Read mcal, truth, and srcext catalogs
    sof   = fitsio.read(sof_path,   ext = 1)
    Truth = fitsio.read(Truth_path, ext = 1)
    Ocat  = [fitsio.read(i, ext = 1) for i in OCat_path]
    Bcat  = [fitsio.read(i, ext = 1) for i in BCat_path]
    
    if len(sof) != len(Bcat[0]):
        print("SOF NOT SAME LENGTH AS BALROG SE CAT. SKIPPING TILE", tilename)
        return None

    #STEP 1: match SrcExtractor objects with injected objects. Bcat[0] is r-band
    tree = BallTree(np.vstack([Ocat[0]['DELTAWIN_J2000'],   Ocat[0]['ALPHAWIN_J2000']]).T * np.pi/180, leaf_size=40, metric="haversine")
    d, j = tree.query(np.vstack([Bcat[0]['DELTAWIN_J2000'], Bcat[0]['ALPHAWIN_J2000']]).T * np.pi/180)

    d, j = d[:, 0], j[:, 0] #convert to 1d array
    d    = d * 180/np.pi * 60*60 #convert to arcsec
    
    #Keep only objects that have no matches below 1 arcsec
    Mask = d > 1
    j    = j[Mask]
    Nobj = sum(Mask)
    
    #STEP 2: Find the nearest injected object to each one
    tree   = BallTree(np.vstack([Truth['dec'], Truth['ra']]).T * np.pi/180, leaf_size=40, metric="haversine")
    d2, j2 = tree.query(np.vstack([Bcat[0]['DELTAWIN_J2000'], Bcat[0]['ALPHAWIN_J2000']]).T * np.pi/180)

    d2, j2 = d2[:, 0], j2[:, 0] #convert to 1d array
    d2     = d2 * 180/np.pi * 60*60 #convert to arcsec    
    
    #STEP 3: Construct the catalog
        
    #declare type of the output array
    dtype  = [('meas_detected', int), ('meas_d_truth_arcsec', float), ('truth_ra', float), ('truth_dec', float)]
    dtype += [(f'truth_{name}', dtype) for name, dtype in Input_catalog.dtype.descr]
    dtype += [(f'meas_{X[0]}', X[1]) if len(X) == 2 else (f'meas_{X[0]}', X[1], X[2]) for X in sof.dtype.descr]
    
    for b in 'GRIZ':
        dtype += [(f'meas_{X[0]}_{b}', X[1]) if len(X) == 2 else (f'meas_{X[0]}_{b}', X[1], X[2]) for X in Bcat[0].dtype.descr]
    
    dtype  = np.dtype(dtype)
    output = np.zeros(Nobj, dtype = dtype)
    
    for n in sof.dtype.names:
        output[f'meas_{n}'] = sof[n][Mask]
    
    for n in Bcat[0].dtype.names:
        for b, B in zip('GRIZ', Bcat):
            output[f'meas_{n}_{b}'] = B[n][Mask]

    for n in Input_catalog.dtype.names:
        output[f'truth_{n}'] = Input_catalog[n][Truth['ind']][j2][Mask]

    output['truth_ra']             = Truth['ra'][j2][Mask]
    output['truth_dec']            = Truth['dec'][j2][Mask]
    output['meas_d_truth_arcsec']  = d2[Mask]
    
    return output


if __name__ == "__main__":
    

    name     = os.path.basename(os.path.dirname(__file__))
    BROG_DIR = os.environ['BALROG_DIR']
    PATH     = BROG_DIR + '/' + name
    config   = yaml.load(open(os.path.dirname(__file__) + '/config.yaml', 'r'), Loader=yaml.Loader)
    print('GETTING BALROG FILES FROM:', PATH)
    
    files = sorted(glob.glob(PATH + '/fitvd_*'))
    tiles = [f[-17:-5] for f in files] #Tilenames

    print(f"I HAVE FOUND {len(files)} FILES")
    
    FINAL_CAT = [None] * len(files)
    tilenames = [None] * len(files)
    
    def my_func(i):
        f = files[i]
        tile = os.path.basename(f).split('_')[1].split('.')[0]
        cat  = match_catalogs(os.path.dirname(f), tile, 'griz', config)
        
        return i, cat, [tile] * len(cat)
        
    jobs = [joblib.delayed(my_func)(i) for i in range(len(files))]

    with joblib.parallel_backend("loky"):
        outputs = joblib.Parallel(n_jobs = 1, verbose=10)(jobs)
        
        for o in outputs:
            if o is None: continue
            FINAL_CAT[o[0]] = o[1]
            tilenames[o[0]] = o[2]
                    
    FINAL_CAT = np.concatenate([f for f in FINAL_CAT if f is not None], axis = 0)
    tilenames = np.concatenate([t for t in tilenames if t is not None], axis = 0)
    
    BITMASK = hp.read_map('/project/chihway/dhayaa/DECADE/Foreground_Masks/GOLD_Ext0.2_Star5_MCs2.fits')
    bmask   = BITMASK[hp.ang2pix(hp.npix2nside(BITMASK.size), FINAL_CAT['truth_ra'], FINAL_CAT['truth_dec'], lonlat = True)]

    with h5py.File(PATH + '/BalrogOfTheDwarves_DiffCatalog_20250406.hdf5', 'w') as f:
    
        for i in tqdm(FINAL_CAT.dtype.names, desc = 'Making HDF5'):

            if i in ['meas_flagstr']: continue #Ignore some fields
            f.create_dataset(i, data = FINAL_CAT[i])
        
        
        f.create_dataset('FLAGS_FOREGROUND', data = bmask); print("FINISHED ADDING SUPPLEMENTARY DATA: FLAG_FOREGROUND")            
        f.create_dataset('tilename', data = tilenames.astype('S'), dtype = h5py.special_dtype(vlen=str)); print("FINISHED ADDING SUPPLEMENTARY DATA: TILENAME")
        
        
        #Deredden quantities
        for name in ['SFD98', 'Planck13']:

            if name == 'SFD98':
                EXTINCTION = hp.read_map('/project/chihway/dhayaa/DECADE/Extinction_Maps/ebv_sfd98_nside_4096_ring_equatorial.fits')
                R_SFD98    = EXTINCTION[hp.ang2pix(4096, f['meas_ra'][:], f['meas_dec'][:], lonlat = True)]
                Ag, Ar, Ai, Az = R_SFD98*3.186, R_SFD98*2.140, R_SFD98*1.569, R_SFD98*1.196

            elif name == 'Planck13':
                EXTINCTION = hp.read_map('/project/chihway/dhayaa/DECADE/Extinction_Maps/ebv_planck13_nside_4096_ring_equatorial.fits')
                R_PLK13    = EXTINCTION[hp.ang2pix(4096, f['meas_ra'][:], f['meas_dec'][:], lonlat = True)]
                Ag, Ar, Ai, Az = R_PLK13*4.085, R_PLK13*2.744, R_PLK13*2.012, R_PLK13*1.533

            f.create_dataset('Ag_' + name.lower(), data = Ag)
            f.create_dataset('Ar_' + name.lower(), data = Ar)
            f.create_dataset('Ai_' + name.lower(), data = Ai)
            f.create_dataset('Az_' + name.lower(), data = Az); print(f"FINISHED ADDING SUPPLEMENTARY DATA: EXTINCTION f{name}")

    print(f"FINISHED MAKING CATALOG")
        
    
