import os, io, subprocess as sp
import jinja2
import yaml
import pandas as pd
import argparse
import datetime as dt
import numpy as np
import glob
import time


my_parser = argparse.ArgumentParser()

my_parser.add_argument('--Initialize', action='store_true', default = False)
my_parser.add_argument('--Maintain',   action='store_true', default = False)
my_parser.add_argument('--TileStart',  action='store', type = int, default = 0)
my_parser.add_argument('--TileEnd',    action='store', type = int, default = 10)
my_parser.add_argument('--Nseeds',     action='store', type = int, default = 1)

my_parser.add_argument('--MaxConcurrentJobs', action='store', type = int, default = 10)
my_parser.add_argument('--MaxCutoffTime',     action='store', type = int, default = 3600*48) #Maxcutoff time in seconds

args = vars(my_parser.parse_args())

#Print args for debugging state
print('-------INPUT PARAMS----------')
for p in args.keys():
    print('%s : %s'%(p.upper(), args[p]))
print('-----------------------------')
print('-----------------------------')

TILENAME_SEED = 42
tiles = pd.read_csv(os.environ['BALROG_RUN_DIR'] + '/data/Tilelist_DR3_All.csv')
tilenames = list(np.random.default_rng(TILENAME_SEED).choice(tiles['TILENAME'].values, size = len(tiles), replace = False))[args['TileStart']:args['TileEnd']]

print(tilenames)
if __name__ == '__main__':
    
    #Automatically get folder name (assuming certain folder structure)
    name = os.path.basename(os.path.dirname(__file__))

    #Create output directory for metacal
    BALROG_DIR = os.environ['BALROG_DIR']
    os.makedirs(BALROG_DIR +'/'+ name, exist_ok=True)

        
    ########################################################################3
    
    def is_finished(tilename, confignum, seed):
        
        args  = {'name' : name, 'tile' : tilename, 'confignum' : confignum, 'seed' : seed}
        plus  = os.path.join(os.environ['BALROG_DIR'], "/%(name)s_v%(confignum)d/fitvd_%(tile)s_seed%(seed)d.fits" % args)
        
        return os.path.isfile(plus)
    
    def is_job(tilename, confignum, seed):
        
        return os.path.isfile('job_%s_v%d_s%d.sh' % (tilename, confignum, gal_seed))
    
            
    def create_job(tilename, confignum, gal_seed):
        
        #Now create all job.sh files for running sims
        with open('job.sh.temp', 'r') as fp:
            tmp = jinja2.Template(fp.read())
        
        if not is_finished(tilename, confignum, gal_seed):
            with open('job_%s_v%d_s%d.sh' % (tilename, confignum, gal_seed), 'w') as fp:
                fp.write(tmp.render(tilename=tilename, model_name = name, seed_galsim=gal_seed, seed_mcal=42, config_num = confignum))
            os.system('chmod u+x job_%s_v%d_s%d.sh' % (tilename, confignum, gal_seed))
        
    def prep_job(tilename, confignum, gal_seed):
        
        print("-------------------------------")

        commands = []
        for b in 'griz':
            commands.append('python $BALROG_RUN_DIR/run_sims.py prep --tilename="%s" --bands="%s" --output-desdata="$PREP_DIR/%s_config%s/outputs_%s_seed%s" --config-file="config0.yaml"'%(tilename, b, name, confignum, tilename, gal_seed))
        
        # Run all commands in parallel using subprocess.Popen
        processes = [sp.Popen(command, shell = True) for command in commands]

        # Wait for all processes to complete
        for process in processes: process.wait()
            
        print("-------------------------------")
        print("DONE PREPPING")
        print("-------------------------------")
        
        
    def current_job_count():
        
        x = sp.check_output("squeue --format='%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R' --sort=+i -u dhayaa", shell = True)
        j = pd.read_csv(io.StringIO(x.decode("utf-8")), delim_whitespace=True)
        
        count = 0
        tiles = []
        
        #metacal_DES0849+0252_seed0_gminus
        for i in range(len(j)):
            
            condition1 = 'balrog' in j['NAME'].values[i]
            condition2 = j['STATE'].values[i] == 'RUNNING'
            
            if condition1:
                count += 1
                tiles.append(j['NAME'].values[i][7:])
                
        for t in tiles:
            print("CURRENTLY RUNNING: %s" % t)
                
        return count
        
    #job_DES1916-1624_v0_s1763835248.sh
        
    def tilename_from_jobfilename(job): return job[4:4+12]
        
    def confignum_from_jobfilename(job): return job[18:18+1]
    
    def seednum_from_jobfilename(job): return job[21:-3]    
                
    ################################################################
    
    #In initial step we just create every job that we need
    if args['Initialize']:
        
        Nconfigs = len(glob.glob('config*.yaml'))
        
        for i, tilename in enumerate(tilenames):
            
            tile_seed = np.random.default_rng(seed = i).integers(2**32, size = Nconfigs)
            
            for confignum in range(Nconfigs):    
                
                seeds = np.random.default_rng(seed = tile_seed[confignum]).integers(2**32, size = args['Nseeds'])
                
                for seed in seeds:

                    create_job(tilename, confignum, seed)
    
    #In next step we keep track of all jobs and add to queue when we can
    elif args['Maintain']:
        
        job_list = sorted(glob.glob('job_*.sh'))
        start = dt.datetime.now()

        while (dt.datetime.now() - start).seconds < 3600*96: #4 days max
            
            if len(job_list) == 0:
                print("ALL JOBS HAVE BEEN STARTED. NOTHING TO MAINTAIN")
                break
                
            if current_job_count() >= args['MaxConcurrentJobs']:
                
                print("---------------------------")
                print(dt.datetime.now())
                print("---------------------------") 
                time.sleep(60*5)

            else:
                
                j = job_list[0]
                t = tilename_from_jobfilename(j)
                c = confignum_from_jobfilename(j)
                s = seednum_from_jobfilename(j)
                
                prep_job(t, c, s)
                
                os.system('sbatch %s' % j)
                os.system('rm %s' % j)
                
                print("SUBMITTED JOB %s" % j)
                
                #This just removes the top item from the list
                #which we do as we just ran that job
                job_list = job_list[1:]
                
        else:
            
            print("GONE OVERTIME. SHUTTING DOWN SCRIPT")            
                
    
