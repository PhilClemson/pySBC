import argparse
import os
import subprocess
import numpy as np
import time
from math import ceil
import copy
import glob
from numpy.random import randint
from threading import Timer
import signal

def parse_args():
    parser = argparse.ArgumentParser(description="Perform Simulation Based Calibration.")
    parser.add_argument("--exe", dest="exe", action="store",
                        help="Location of compiled Stan executable.")
    parser.add_argument("--nt", dest="num_threads", action="store",
                        help="Number of threads.",type=int)
    parser.add_argument("--warmup", dest="num_warmup", action="store",
                        help="Number of warmup samples.",type=int)
    parser.add_argument("--samples", dest="num_samples", action="store",
                        help="Number of samples after warmup.",type=int)
    parser.add_argument("--data", dest="data_file", action="store",
                        help="Location of data file.")
    parser.add_argument("--J", dest="num_bins", action="store",
                        help="Number of bins in the rank statistic histogram.",type=int)
    parser.add_argument("--ej", dest="expected_bin_count", action="store",
                        help="Expected number in each rank statistic bin (if uniformly distributed). Should be no less than 5.",type=int)
    parser.add_argument("--output", dest="output_file", action="store",
                        help="Name of output file.")
    parser.add_argument("--t", dest="time_out", action="store",
                        help="Maximum time to wait for a chain to complete in seconds.",type=int)
    return parser.parse_args()

class chain:
    # store information about the various chains here
    def __init__(self, exe, num_samples, num_warmup, init, thinning_factor, data_file, time_out, n):
        self.base_command = "{} method=sample num_samples={} num_warmup={} thin={} data file={} output file=sbc_output/sbc_output{}.out refresh=0 init={}".format(exe, num_samples, num_warmup, thinning_factor, data_file, n, init)
        self.time_out = time_out
        self.ind = n
        
    def generate_seed_and_run(self):
        # generate random seed
        #seed = randint(0,65535,1)
        seed = randint(0,2147483647,1)
        self.seed = seed[0]
        self.shexec(self.base_command + " random seed={}".format(self.seed))
        self.t = Timer(self.time_out, self.shkill, []) # kill if time specified by time_out passes
        self.t.start()
        self.timer_flag = 1
        
    def run_with_seed(self,seed):
        self.seed = seed
        self.shexec(self.base_command + " random seed={}".format(seed))
        self.timer_flag = 0
        
    def shexec(self,command):
        self.process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=os.setsid)
        self.block_poll = 0
    
    def shkill(self):
    	self.block_poll = 1
        try:
            os.killpg(os.getpgid(self.process.pid), signal.SIGTERM)
            print('Chain {} timed out. Re-running with a different seed.'.format(self.ind))
        except:
            print('Chain {} did not exit correctly. Re-running with a different seed.'.format(self.ind))
        self.generate_seed_and_run()
    
    def get_seed(self):
    	return self.seed
    	
    def get_ind(self):
    	return self.ind
    	
    def get_poll(self):
    	# Check we're not in the middle of restarting the chain
    	if self.block_poll == 0:
            return self.process.poll()
        else:
            return None
    
    def get_code(self):
    	# Check we're not in the middle of restarting the chain
    	if self.block_poll == 0:
    	    if self.timer_flag == 1:
    	    	#print("Cancelling chain {}".format(self.ind))
    	        self.t.cancel()
    	        self.timer_flag = 0
    	    return self.process.returncode
    	else:
    	    return None
    	
    def print_output(self):
        print(self.process.stderr.read())
        print(self.process.stdout.read())
    
    def wait(self):
        return self.process.wait()

if __name__ == "__main__":
    args = parse_args()
    exe = args.exe
    num_threads = args.num_threads
    num_warmup = args.num_warmup
    num_samples = args.num_samples
    data_file = args.data_file
    num_bins = args.num_bins
    expected_bin_count = args.expected_bin_count
    output_file = args.output_file
    time_out = args.time_out
    
    # set maximum thinning factor after ESS calculation
    max_thin = 100
    
    # store number of samples prior to changes due to ESS calculation
    num_actual_samples = copy.deepcopy(num_samples)
    
    wd = os.getcwd()
    
    # remove output files from a previous run
    for filename in glob.glob("sbc_output/sbc_output*"):
        os.remove(filename) 
    
    num_sims = int(ceil(float(expected_bin_count * num_bins) / float(num_threads)))
    
    total_sims = int(num_sims*num_threads)
    
    chains = [None]*num_threads
    seeds = [None]*total_sims
    thinning_factor = [None]*total_sims

    used_threads = 0
    
    print("Starting SBC chains...")
    i = 0
    n = 1
    while (used_threads < num_threads):
        # start sampling processes
        chains[i] = chain(exe, num_samples, num_warmup, 2, 1, data_file, time_out, n)
        chains[i].generate_seed_and_run()
        i += 1
        n += 1
        used_threads += 1
    
    ind = 0
    
    flag_finished = 0
    finished_chains = []
    
    # wait for processes to finish and start new ones
    while (flag_finished == 0):
    	flag_finished = 1
        for i in range(0,num_threads):
            poll = chains[i].get_poll()
	    if poll is not None:
		if chains[i].get_code() != 0:
		    flag_finished = 0
		    #chains[i].print_output()
		    print('Chain {} encountered errors. Re-running with a different seed.'.format(chains[i].get_ind()))
		    # bad random seed? try with another one
		    chains[i].generate_seed_and_run()
		else:
		    if (ind < total_sims):
		        chain_ind = chains[i].get_ind()
		        seeds[chain_ind - 1] = copy.deepcopy(chains[i].get_seed())
		        if (n < total_sims + 1):
		            flag_finished = 0
		            ind += 1
		            print("Finished chain {} ({} of {})".format(chain_ind,ind,total_sims))
			    # start new process
			    chains[i] = chain(exe, num_samples, num_warmup, 2, 1, data_file, time_out, n)
			    chains[i].generate_seed_and_run()
			    n += 1
			else:
			    if (chain_ind not in finished_chains):
			        flag_finished = 0
		        	ind += 1
			        print("Finished chain {} ({} of {})".format(chain_ind,ind,total_sims))
			        finished_chains.append(chain_ind)
	    else:
	    	flag_finished = 0
			
    used_threads = 0
    
    # check which columns correspond to rank statistics
    f = open("sbc_output/sbc_output1.out", "r")
    for n in range(0,39):
        f.readline()
        
    col_header = f.readline();
    col_names = col_header.split(',')
    num_cols = len(col_names)
    col_inds = []
    n = 0
    for name in col_names:
        if name[0:6] == "lt_sim":
            col_inds.append(n)
        n += 1
			
    print("Restarting chains with thinning based on ESS...")
    i = 0
    n = 0
    while (used_threads < num_threads):
	f_string = "sbc_output/sbc_output{}.out".format(n+1)
	print("Processing file {}".format(f_string))
	try:
	    lines = np.loadtxt(f_string, comments=["#","lp__"], delimiter=",", unpack=False)
	    min_ESS = num_samples
	    for k in range(7,col_inds[0]):
		samps = lines[:,k]
		var = np.var(samps)
		rho_t = np.ones(num_samples-1)
		t = 1
		while (t < num_samples-1):
		    Vt = np.sum(np.square(samps[t:num_samples-1] - samps[1:num_samples-t]))/(num_samples-t)
		    rho_t[t] = 1 - Vt/(2*var)
		    if (rho_t[t] + rho_t[t-1] < 0):
		    	break;
		    else:
		    	t += 1
		if (t > 3):
		    ESS = num_samples / (1 + 2*np.sum(rho_t[1:t-2]))
		    if ESS < min_ESS:
			min_ESS = ESS
		else:
		    if (t == num_samples-1):
		        # something wrong with autocorrelation, retry with different seed by setting thinning_factor[n]>10
		        min_ESS = num_samples/(max_thin + 1)
	    if (min_ESS > 0):
		thinning_factor[n] = int(round(num_samples/min_ESS))
	    else:
		thinning_factor[n] = max_thin + 1
	    # only attempt for thinning factors < max_thin to prevent glacial computation times
	    while (thinning_factor[n] > max_thin):
	    	print("ESS = {}. Re-running chain with different seed.".format(min_ESS))
		chains[i] = chain(exe, num_samples, num_warmup, 2, 1, data_file, time_out, n+1)
		chains[i].generate_seed_and_run()
		while (chains[i].wait() != 0):
		    time.sleep(0.1)
		chains[i].get_code() # to stop timeout restart
		seeds[n] = copy.deepcopy(chains[i].get_seed())
	        lines = np.loadtxt(f_string, comments=["#","lp__"], delimiter=",", unpack=False)
	        min_ESS = num_samples
	        for k in range(7,col_inds[0]):
		    samps = lines[:,k]
		    var = np.var(samps)
		    rho_t = np.ones(num_samples-1)
		    t = 1
		    while (t < num_samples-1):
		        Vt = np.sum(np.square(samps[t:num_samples-1] - samps[1:num_samples-t]))/(num_samples-t)
		        rho_t[t] = 1 - Vt/(2*var)
		        if (rho_t[t] + rho_t[t-1] < 0):
		        	break;
		        else:
		        	t += 1
		    if (t > 3):
		        ESS = num_samples / (1 + 2*np.sum(rho_t[1:t-2]))
		        if ESS < min_ESS:
			    min_ESS = ESS
		    else:
		        if (t == num_samples-1):
		            # something wrong with autocorrelation, retry with different seed by setting thinning_factor[n]>10
		            min_ESS = num_samples/(max_thin + 1)
		if (min_ESS > 0):
		    thinning_factor[n] = int(round(num_samples/min_ESS))
		else:
		    thinning_factor[n] = max_thin + 1
	    if (thinning_factor[n] < 1):
	    	thinning_factor[n] = 1
	    print("ESS = {}. Running chain for {} sampling iterations with thinning factor of {}.".format(min_ESS, thinning_factor[n]*num_samples, thinning_factor[n]))
	    # start sampling processes
	    chains[i] = chain(exe, thinning_factor[n]*num_samples, num_warmup, 2, thinning_factor[n], data_file, time_out, n+1)
	    chains[i].run_with_seed(seeds[n])
	    n += 1
	    i += 1
	    used_threads += 1
        except:
            print("Error in output file format. Excluding from estimate.".format(f_string))
            n += 1
    
    ind = 0
    
    flag_finished = 0
    finished_chains = []
    
    # wait for processes to finish and start new ones
    while (flag_finished == 0):
        flag_finished = 1
        for i in range(0,num_threads):
            poll = chains[i].get_poll()
	    if poll is not None:
	    	    chain_ind = chains[i].get_ind()
		    if (ind < total_sims + 1):
			    if (n < total_sims):
			        flag_finished = 0
	    	    	        ind += 1
			        print("Finished chain {} ({} of {})".format(chain_ind,ind,total_sims))
			    	# start new process
				f_string = "sbc_output/sbc_output{}.out".format(n+1)
				print("Processing file {}".format(f_string))
				try:
				    lines = np.loadtxt(f_string, comments=["#","lp__"], delimiter=",", unpack=False)
				    min_ESS = num_samples
				    for k in range(7,col_inds[0]):
					samps = lines[:,k]
					var = np.var(samps)
					rho_t = np.ones(num_samples-1)
					t = 1
					while (t < num_samples-1):
					    Vt = np.sum(np.square(samps[t:num_samples-1] - samps[1:num_samples-t]))/(num_samples-t)
					    rho_t[t] = 1 - Vt/(2*var)
					    if (rho_t[t] + rho_t[t-1] < 0):
					    	break;
					    else:
					    	t += 1
					if (t > 3):
					    ESS = num_samples / (1 + 2*np.sum(rho_t[1:t-2]))
					    if ESS < min_ESS:
						min_ESS = ESS
					else:
					    if (t == num_samples-1):
						# something wrong with autocorrelation, retry with different seed by setting thinning_factor[n]>10
						min_ESS = num_samples/(max_thin + 1)
				    if (min_ESS > 0):
					thinning_factor[n] = int(round(num_samples/min_ESS))
				    else:
					thinning_factor[n] = max_thin + 1
	    			    # only attempt for thinning factors < max_thin to prevent glacial computation times
				    while (thinning_factor[n] > max_thin):
				    	print("ESS = {}. Re-running chain with different seed.".format(min_ESS))
					chains[i] = chain(exe, num_samples, num_warmup, 2, 1, data_file, time_out, n+1)
					chains[i].generate_seed_and_run()
					while (chains[i].wait() != 0):
					    time.sleep(0.1)
					chains[i].get_code() # to stop timeout restart
					seeds[n] = copy.deepcopy(chains[i].get_seed())
					lines = np.loadtxt(f_string, comments=["#","lp__"], delimiter=",", unpack=False)
					min_ESS = num_samples
					for k in range(7,col_inds[0]):
					    samps = lines[:,k]
					    var = np.var(samps)
					    rho_t = np.ones(num_samples-1)
					    t = 1
					    while (t < num_samples-1):
						Vt = np.sum(np.square(samps[t:num_samples-1] - samps[1:num_samples-t]))/(num_samples-t)
						rho_t[t] = 1 - Vt/(2*var)
						if (rho_t[t] + rho_t[t-1] < 0):
							break;
						else:
							t += 1
					    if (t > 3):
						ESS = num_samples / (1 + 2*np.sum(rho_t[1:t-2]))
						if ESS < min_ESS:
						    min_ESS = ESS
						else:
						    if (t == num_samples-1):
							# something wrong with autocorrelation, retry with different seed by setting thinning_factor[n]>10
							min_ESS = num_samples/(max_thin + 1)
					if (min_ESS > 0):
					    thinning_factor[n] = int(round(num_samples/min_ESS))
					else:
					    thinning_factor[n] = max_thin + 1
				    if (thinning_factor[n] < 1):
				    	thinning_factor[n] = 1
				    print("ESS = {}. Running chain for {} sampling iterations with thinning factor of {}.".format(min_ESS, thinning_factor[n]*num_samples, thinning_factor[n]))
				    # start sampling processes
				    chains[i] = chain(exe, thinning_factor[n]*num_samples, num_warmup, 2, thinning_factor[n], data_file, time_out, n+1)
				    chains[i].run_with_seed(seeds[n])
				    n += 1
				    i += 1
				except:
				    print("Error in output file format. Excluding from estimate.".format(f_string))
				    n += 1
		    	    else:
				if (chain_ind not in finished_chains):
			    	    flag_finished = 0
		            	    ind += 1
			    	    print("Finished chain {} ({} of {})".format(chain_ind,ind,total_sims))
			    	    finished_chains.append(chain_ind)
	    else:
	        flag_finished = 0
    
    # Bin output into histogram
    print("Generating SBC histogram...")
    bins = np.zeros((num_bins,len(col_inds)))
    for n in range(1,total_sims+1):
        f_string = "sbc_output/sbc_output{}.out".format(n)
        print("Processing file {}".format(f_string))
        try:
            lines = np.loadtxt(f_string, comments=["#","lp__"], delimiter=",", unpack=False)
            samps = lines[:,col_inds]
            r = samps.sum(axis=0)
            col_ind = 0
            for r_m in r:
                bin_ind = int(np.floor(float(r_m/((num_actual_samples+1.0)/num_bins))))
                if (bin_ind > num_bins - 1):
                    bin_ind = num_bins - 1
                bins[bin_ind,col_ind] += 1
                col_ind += 1
        except:
            print("Error in output file format. Excluding from estimate.".format(f_string))
    
    np.savetxt('{}'.format(output_file), bins, delimiter=',')
    print("Finished. Histogram saved to {}.".format(output_file))
    
