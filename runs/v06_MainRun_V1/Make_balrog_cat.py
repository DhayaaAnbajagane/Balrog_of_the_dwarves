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

def match_catalogs(path, tilename, seed, bands, config):

    Input_catalog = fitsio.read(r'%s/SplicedSim_Input_%s_seed%d-cat.fits' % (path, tilename, seed), ext = 1)
     
    #Get all paths
    sof_path   = r'%s/fitvd_%s_seed%d.fits' % (path, tilename, seed)
    Truth_path = r'%s/Input_%s_seed%d-cat.fits' % (path, tilename, seed)
    OCat_path  = [r'%s/OldSrcExtractor_%s_%s-cat.fits' % (path, tilename, band) for band in bands] #Path to original SrcExtractor
    BCat_path  = [r'%s/SrcExtractor_%s_seed%d_%s-cat.fits' % (path, tilename, seed, band) for band in bands] #Path to new SrcExtractor
    
    print(path, tilename, seed)
    #Read mcal, truth, and srcext catalogs
    sof   = fitsio.read(sof_path,   ext = 1)
    Truth = fitsio.read(Truth_path, ext = 1)
    Ocat  = [fitsio.read(i, ext = 1) for i in OCat_path]
    Bcat  = [fitsio.read(i, ext = 1) for i in BCat_path]
    cnum  = int(path[-1])
    
    if len(Truth) == 0:
        print("NO INJECTIONS IN TILE", tilename)
        return None
         
    if len(sof) != len(Bcat[0]):
        print("SOF NOT SAME LENGTH AS BALROG SE CAT. SKIPPING TILE", tilename)
        return None
    
    if len(sof) == 0:
        print("NO DETECTIONS????? IN TILE", tilename)
        return None

    #STEP 1: match SrcExtractor objects with injected objects. Bcat[0] is r-band
    tree = BallTree(np.vstack([Bcat[0]['DELTAWIN_J2000'], Bcat[0]['ALPHAWIN_J2000']]).T * np.pi/180, leaf_size=40, metric="haversine")
    d, j = tree.query(np.vstack([Truth['dec'], Truth['ra']]).T * np.pi/180)

    d, j = d[:, 0], j[:, 0] #convert to 1d array
    d    = d * 180/np.pi * 60*60 #convert to arcsec
    
    #Keep only ids below 0.5 arcsec
    Mask = d < 0.5
    j    = j[Mask]
    Nobj = len(Mask)
    
    #STEP 2: Take old, original SrcExtractor, for each truth object ask how close a nearby object is.
    tree   = BallTree(np.vstack([Ocat[0]['DELTAWIN_J2000'], Ocat[0]['ALPHAWIN_J2000']]).T * np.pi/180, leaf_size=40, metric="haversine")
    d2, j2 = tree.query(np.vstack([Truth['dec'], Truth['ra']]).T * np.pi/180)

    d2, j2 = d2[:, 0], j2[:, 0] #convert to 1d array
    d2     = d2 * 180/np.pi * 60*60 #convert to arcsec
    
    
    #STEP 3: Construct the catalog
        
    #declare type of the output array
    dtype  = [('meas_detected', int), ('meas_d_contam_arcsec', float), ('truth_ra', float), ('truth_dec', float), ('truth_uniq', int),
              ('seed', '>i8'), ('cnum', '>i8')]
    dtype += [(f'truth_{name}', dtype) for name, dtype in Input_catalog.dtype.descr]
    dtype += [(f'meas_{X[0]}', X[1]) if len(X) == 2 else (f'meas_{X[0]}', X[1], X[2]) for X in sof.dtype.descr]
    
    for b in 'GRIZ':
        dtype += [(f'meas_{X[0]}_{b}', X[1]) if len(X) == 2 else (f'meas_{X[0]}_{b}', X[1], X[2]) for X in Bcat[0].dtype.descr]
    dtype  = np.dtype(dtype)
    output = np.zeros(Nobj, dtype = dtype)
    output['seed'] = seed
    output['cnum'] = cnum
    
    assert np.all(Truth['ID'] == Input_catalog['ID'][Truth['ind']]), "Something went wrong in matching"
    
    for n in Input_catalog.dtype.names:
        output[f'truth_{n}'] = Input_catalog[n][Truth['ind']]
    
    for n in sof.dtype.names:
        output[f'meas_{n}'][Mask] = sof[n][j]
    
    for n in Bcat[0].dtype.names:
        for b, B in zip('GRIZ', Bcat):
            output[f'meas_{n}_{b}'][Mask] = B[n][j]

    output['truth_ra']             = Truth['ra']
    output['truth_dec']            = Truth['dec']
    output['truth_uniq']           = Truth['uniq']
    output['meas_detected']        = Mask.astype(int)
    output['meas_d_contam_arcsec'] = d2
    
    dtype = [('seed', '>i8'), ('cnum', '>i8')]
    for b in 'GRIZ':
        dtype += [(f'{X[0]}_{b}', X[1]) if len(X) == 2 else (f'meas_{X[0]}_{b}', X[1], X[2]) for X in Bcat[0].dtype.descr]
    
    Bout  = np.zeros(len(Bcat[0]), dtype = dtype)
    Bout['seed'] = seed
    Bout['cnum'] = cnum
    for n in Bcat[0].dtype.names:
        for b, B in zip('GRIZ', Bcat):
            Bout[f'{n}_{b}'] = B[n]
            
    dtype = Truth.dtype.descr + Input_catalog.dtype.descr[1:] + [('seed', '>i8'), ('cnum', '>i8')]
    Tout  = np.zeros(len(Truth), dtype = dtype)
    Tout['seed'] = seed
    Tout['cnum'] = cnum
    for n in Truth.dtype.names:         Tout[n] = Truth[n]
    for n in Input_catalog.dtype.names: Tout[n] = Input_catalog[n][Truth['ind']]
            
    return output, Tout, Bout


if __name__ == "__main__":
    

    name     = os.path.basename(os.path.dirname(__file__))
    BROG_DIR = os.environ['BALROG_DIR']
    PATH     = BROG_DIR + '/' + name
    config   = None
    print('GETTING BALROG FILES FROM:', PATH)
    
    files = sorted(glob.glob(PATH + '_config*' + '/fitvd_*'))
    
    print(f"I HAVE FOUND {len(files)} FILES")
    
    FINAL_CAT = [None] * len(files)
    TRUTH_CAT = [None] * len(files)
    BOUT_CAT  = [None] * len(files)
    tilenames = [None] * len(files)
    
    def my_func(i):
        f = files[i]
        tile = os.path.basename(f).split('_')[1].split('.')[0]
        seed = int(os.path.basename(f).split('_')[2].split('.')[0][4:])
        cat  = match_catalogs(os.path.dirname(f), tile, seed, 'griz', config)
        N    = len(cat) if cat is not None else 0 
        return i, cat, [tile] * N
        
    jobs = [joblib.delayed(my_func)(i) for i in range(len(files))]

    with joblib.parallel_backend("loky"):
        outputs = joblib.Parallel(n_jobs = 2, verbose=10)(jobs)
        
        for o in outputs:
            if o is None: continue
            FINAL_CAT[o[0]] = o[1][0]
            TRUTH_CAT[o[0]] = o[1][1]
            BOUT_CAT[o[0]]  = o[1][2]
            tilenames[o[0]] = o[2]
                    
    FINAL_CAT = np.concatenate([f for f in FINAL_CAT if f is not None], axis = 0)
    TRUTH_CAT = np.concatenate([f for f in TRUTH_CAT if f is not None], axis = 0)
    BOUT_CAT  = np.concatenate([f for f in BOUT_CAT  if f is not None], axis = 0)
    tilenames = np.concatenate([t for t in tilenames if t is not None], axis = 0)
    
    BITMASK = hp.read_map('/project/chihway/dhayaa/DECADE/Foreground_Masks/GOLD_Ext0.2_Star5_MCs2.fits')
    bmask   = BITMASK[hp.ang2pix(hp.npix2nside(BITMASK.size), FINAL_CAT['truth_ra'], FINAL_CAT['truth_dec'], lonlat = True)]

    with h5py.File(PATH + '/BalrogOfTheDwarves_Catalog_20250512.hdf5', 'w') as b:
        
        f = b.create_group('Truth')
        for i in tqdm(TRUTH_CAT.dtype.names, desc = 'Making TRUTH HDF5'):
            f.create_dataset(i, data = TRUTH_CAT[i])
            
        f = b.create_group('Output')
        for i in tqdm(BOUT_CAT.dtype.names, desc = 'Making OUTPUT HDF5'):
            f.create_dataset(i, data = BOUT_CAT[i])
        
        
        f = b.create_group('Matched')
        for i in tqdm(FINAL_CAT.dtype.names, desc = 'Making MATCHED HDF5'):

            if i in ['meas_flagstr']: continue #Ignore some fields
            f.create_dataset(i, data = FINAL_CAT[i])
        
        
        f.create_dataset('FLAGS_FOREGROUND', data = bmask); print("FINISHED ADDING SUPPLEMENTARY DATA: FLAG_FOREGROUND")            
        f.create_dataset('tilename', data = tilenames.astype('S'), dtype = h5py.special_dtype(vlen=str)); print("FINISHED ADDING SUPPLEMENTARY DATA: TILENAME")
        
        
        #Deredden quantities
        for name in ['SFD98', 'Planck13']:

            if name == 'SFD98':
                EXTINCTION = hp.read_map('/project/chihway/dhayaa/DECADE/Extinction_Maps/ebv_sfd98_nside_4096_ring_equatorial.fits')
                R_SFD98    = EXTINCTION[hp.ang2pix(4096, f['truth_ra'][:], f['truth_dec'][:], lonlat = True)]
                Ag, Ar, Ai, Az = R_SFD98*3.186, R_SFD98*2.140, R_SFD98*1.569, R_SFD98*1.196

            elif name == 'Planck13':
                EXTINCTION = hp.read_map('/project/chihway/dhayaa/DECADE/Extinction_Maps/ebv_planck13_nside_4096_ring_equatorial.fits')
                R_PLK13    = EXTINCTION[hp.ang2pix(4096, f['truth_ra'][:], f['truth_dec'][:], lonlat = True)]
                Ag, Ar, Ai, Az = R_PLK13*4.085, R_PLK13*2.744, R_PLK13*2.012, R_PLK13*1.533

            f.create_dataset('Ag_' + name.lower(), data = Ag)
            f.create_dataset('Ar_' + name.lower(), data = Ar)
            f.create_dataset('Ai_' + name.lower(), data = Ai)
            f.create_dataset('Az_' + name.lower(), data = Az); print(f"FINISHED ADDING SUPPLEMENTARY DATA: EXTINCTION f{name}")

    print(f"FINISHED MAKING CATALOG")
        
    
