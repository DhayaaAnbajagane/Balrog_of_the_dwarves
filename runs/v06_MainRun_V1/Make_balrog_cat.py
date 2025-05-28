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
    
    #v06_MainRun_V1_config3/SrcExtractor_DES2111-0207_seed2761987536_z-cat.fits
    #v06_MainRun_V1_config3/Input_DES2111-0207_seed2091518858-cat.fits
    #v06_MainRun_V1_config3/fitvd_DES1108+1751_seed805984730.fits
    #v06_MainRun_V1_config3/SplicedSim_Input_DES1916-1624_seed1953607505-cat.fits
    
    name     = os.path.basename(os.path.dirname(__file__))
    BROG_DIR = os.environ['BALROG_DIR']
    PATH     = BROG_DIR + '/' + name
    config   = None
    print('GETTING BALROG FILES FROM:', PATH)
    
    files = sorted(glob.glob(PATH + '_config*' + '/fitvd_*'))
    
    
    fitvd = fitsio.read(files[0])
    SE    = fitsio.read(files[0].replace('fitvd_', 'SrcExtractor_').replace('.fits', '_r-cat.fits'))
    Input = fitsio.read(files[0].replace('fitvd_', 'Input_').replace('.fits', '-cat.fits'))
    Splic = fitsio.read(files[0].replace('fitvd_', 'SplicedSim_Input_').replace('.fits', '-cat.fits'))
    
    
    def fitvd_it(key): 
        
        if key == 'tilename':
            f = lambda i : [os.path.basename(files[i]).split('_')[1].split('.')[0]] * fitsio.read(files[i], columns = 'id').size
        elif key == 'seed':
            f = lambda i : [int(os.path.basename(files[i]).split('_')[2].split('.')[0][4:])] * fitsio.read(files[i], columns = 'id').size
        elif key == 'cnum':
            f = lambda i : [int(os.path.dirname(files[i])[-1])] * fitsio.read(files[i], columns = 'id').size
        else:
            f = lambda i : fitsio.read(files[i], columns = key)
        
        return np.concatenate(
                    joblib.Parallel(n_jobs = -1, verbose=10)(
                                    [joblib.delayed(f)(i) for i in range(len(files))]
                                )
                        )
    
    
    def SE_it(key, band): 
        
        return np.concatenate(
                    joblib.Parallel(n_jobs = -1, verbose=10)(
                                    [joblib.delayed(lambda i : fitsio.read((files[i]
                                                                            .replace('fitvd_', 'SrcExtractor_')
                                                                            .replace('.fits', f'_{band}-cat.fits')), columns = key))(i) 
                                     for i in range(len(files))]
                                )
                        )
    
    
    def Truth_it(key): 
        
        if key == 'tilename':
            f = lambda i : ([os.path.basename(files[i]).split('_')[1].split('.')[0]] * 
                            fitsio.read(files[i].replace('fitvd_', 'Input_').replace('.fits', f'-cat.fits'), columns = 'id').size)
        elif key == 'seed':
            f = lambda i : ([int(os.path.basename(files[i]).split('_')[2].split('.')[0][4:])] * 
                            fitsio.read(files[i].replace('fitvd_', 'Input_').replace('.fits', f'-cat.fits'), columns = 'id').size)
        elif key == 'cnum':
            f = lambda i : ([int(os.path.dirname(files[i])[-1])] * 
                            fitsio.read(files[i].replace('fitvd_', 'Input_').replace('.fits', f'-cat.fits'), columns = 'id').size)
        else:
            f = lambda i : fitsio.read(files[i].replace('fitvd_', 'Input_').replace('.fits', f'-cat.fits'), columns = key)
            
        return np.concatenate(
                    joblib.Parallel(n_jobs = -1, verbose=10)(
                                    [joblib.delayed(f)(i) for i in range(len(files))]
                                )
                        )
    
    
    def Splic_it(key): 
        
        def my_matched_function(i):
            
            Input = fitsio.read(files[i].replace('fitvd_', 'Input_').replace('.fits', '-cat.fits'), columns = 'ind')
            Splic = fitsio.read(files[i].replace('fitvd_', 'SplicedSim_Input_').replace('.fits', '-cat.fits'), columns = key)
            Splic = Splic[Input]
            
            return Splic
            
        return np.concatenate(
                    joblib.Parallel(n_jobs = -1, verbose=10)(
                                    [joblib.delayed(my_matched_function)(i) for i in range(len(files))]
                                )
                        )
    
    
    with h5py.File(PATH + '/BalrogOfTheDwarves_Catalog_20250527.hdf5', 'w') as b:
        
        f = b.create_group('Truth')
        for i in tqdm(Input.dtype.names, desc = 'Adding Truth'):
            if i in ['ID']: continue
            f.create_dataset(i, data = Truth_it(i))
            
        f.create_dataset('tilename', data = Truth_it('tilename').astype('S'), dtype = h5py.special_dtype(vlen=str))
        f.create_dataset('cnum',     data = Truth_it('cnum'))
        f.create_dataset('seed',     data = Truth_it('seed'))
            
        for i in tqdm(Splic.dtype.names, desc = 'Adding Input'):
            f.create_dataset(i, data = Splic_it(i))
            
        f = b.create_group('Output')
        for i in tqdm(fitvd.dtype.names, desc = 'Adding fitvd'):
            if i in ['flagstr']: continue
            f.create_dataset(i, data = fitvd_it(i))
            
        f.create_dataset('tilename', data = fitvd_it('tilename').astype('S'), dtype = h5py.special_dtype(vlen=str))
        f.create_dataset('cnum',     data = fitvd_it('cnum'))
        f.create_dataset('seed',     data = fitvd_it('seed'))
            
        for band in 'GRIZ':
            for i in tqdm(SE.dtype.names, desc = f'Adding SourceExtractor BAND {band}'):
                f.create_dataset(i + f'_{band}', data = SE_it(i, band = band.lower()))
        
        
    print(f"FINISHED MAKING CATALOG")
        
    
