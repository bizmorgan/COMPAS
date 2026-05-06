###################################################################
#                                                                 #                                                               
#  Plotting the Spin evolution of a COMPAS run                    #
#                                                                 #
################################################################### 
#when running this file, need to sort that it only runs for BBHs 

import os
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from matplotlib import rcParams
import tempfile
from pathlib import Path
from astropy import units as u
from astropy import constants as const

class determine_bbh(object):
    """
    This class determines which systems in the COMPAS output form Binary Black Holes (BBHs)
    In order to use them for future spin calculations 
    """

    def __init__(self, data_path):
        data_path = Path(data_path)
        self.path = data_path

        if not self.path.is_file():
            raise ValueError(
                f"h5 file not found. Wrong path given?\npath checked = {self.path}")
        else:
            self.h5file = h5.File(self.path, "r")

    def identify_bbh_systems(self):
        """
        Identify systems in the main COMPAS output that evolve to
        binary black holes (BBHs). Stores the corresponding SEEDs and then returns a list of the paths to the detailed output files
        """

        fDCO = self.h5file['BSE_Double_Compact_Objects']

        # SEEDs for all DCO systems
        seed_dco = fDCO['SEED'][...].squeeze()

        # Final compact object stellar types
        st1 = fDCO['Stellar_Type(1)'][...].squeeze()
        st2 = fDCO['Stellar_Type(2)'][...].squeeze()

        # BBH condition: both compact objects are black holes, or neutron stars
        mask_binary_BH = (st1 == 14) & (st2 == 14) | (st1 == 13) & (st2 == 13) | (st1 == 14) & (st2 == 13) | (st1 == 13) & (st2 == 14)

        #mask_binary_BH = (st1 == 14) & (st2 == 14) #just looking at BBHs for now, can add NSs back in later if want to

        # Store BBH seeds for later use with detailed output files
        self.BBH_SEEDS = seed_dco[mask_binary_BH]
        #print(self.BBH_SEEDS.dtype, self.BBH_SEEDS) ######## removing - want to print the overlap. 
        #print(f"Identified {self.BBH_SEEDS.size} BBH systems.")

        if self.BBH_SEEDS.size == 0:
            raise RuntimeError("No BBH systems found in this COMPAS output.")

        # detailed_dir = self.path.parent / 'Detailed_Output'
        # print("Detailed dir exists:", detailed_dir.exists())

        # # removing this for now to use the determine_hmxrb object below 

        # fRLOF = self.h5file['BSE_RLOF']

        # #seed for these systems 
        # seeds_rlof = fRLOF['SEED'][...].squeeze()
        # #store seeds for later use 
        # self.seed_rlof = np.unique(seeds_rlof) 
        # print(f"Identified {self.seed_rlof.size} HMXRB systems.")

        # #Combining requirements for HMXRB and BBH
        # self.combined_seeds = np.intersect1d(self.BBH_SEEDS, self.seed_rlof)
        # print(f"{self.combined_seeds.size, self.combined_seeds} systems are BOTH BBH and HMXRB.")

        # if self.combined_seeds.size == 0:
        #     raise RuntimeError("No matching systems found.")

        return self.BBH_SEEDS

    
    def close(self):
        self.h5file.close()

class determine_hmxrb(object):

    def __init__(self, data_path):
        data_path = Path(data_path)
        self.path = data_path
        if not self.path.is_file():
            raise ValueError(
                f"h5 file not found. Wrong path given?\npath checked = {self.path}")
        else:
            self.h5file = h5.File(self.path, "r")
        
    def identify_hmxrb_systems(self): 
        fRLOF = self.h5file['BSE_RLOF']

        #seed for these systems 
        seeds_rlof = fRLOF['SEED'][...].squeeze()

        st1 = fRLOF['Stellar_Type(1)<MT'][...].squeeze()
        st2 = fRLOF['Stellar_Type(2)<MT'][...].squeeze() 

        pre_rad = fRLOF['Radius(2)|RL<step'][...].squeeze()
        post_rad = fRLOF['Radius(2)|RL>step'][...].squeeze()

        rlof1 = fRLOF['RLOF(1)<MT'][...].squeeze() #want the timestep before mass transfer is occuring, needs to be 1 in the mask to be true 
        rlof2 = fRLOF['RLOF(2)<MT'][...].squeeze()

        time = fRLOF['Time<MT'][...].squeeze() 

        mass1 = fRLOF['Mass(1)<MT'][...].squeeze()
        mass2 = fRLOF['Mass(2)<MT'][...].squeeze()

        sep_all  = fRLOF['SemiMajorAxis<MT'][...].squeeze()
        rad1_all = fRLOF['Radius(1)<MT'][...].squeeze()
        rad2_all = fRLOF['Radius(2)<MT'][...].squeeze()

        mt_rates_all = fRLOF['MassTransferRateDonor'][...].squeeze()

        #Inlcude case A MT here also 
        
        caseA_mask = (
            ((rlof1 == 1) & (st1 <= 2)) |
            ((rlof2 == 1) & (st2 <= 2))
        ) 


        #ratio = post_rad / pre_rad #needs to be >0.8 this is for the secondary star as the first will already be a black hole. Only considering compainion star is star 2

        safe_ratio = np.divide(
            post_rad, pre_rad,
            out=np.full_like(post_rad, np.nan, dtype=float),
            where=pre_rad > 0
        )

        compact1 = (st1 >= 13) & (st2 < 13)
        compact2 = (st2 >= 13) & (st1 < 13)

        hmxrb_mask = (compact1 | compact2) & (safe_ratio > 0.8)
                
        # Case A MT events
        caseA_seeds = seeds_rlof[caseA_mask]
        caseA_times = time[caseA_mask]

        # HMXRB events
        hmxb_seeds_rows = seeds_rlof[hmxrb_mask]
        hmxb_times_rows = time[hmxrb_mask]

        self.hmxrb_seeds = np.unique(hmxb_seeds_rows) 

        #print(self.hmxrb_seeds) #these have the rlof 0.8/0.9 but not the case A MT 

        valid_hmxb_seeds = []
        visible_times = {}
        
        #check hmxrb events had an earlier case A MT event 
        for seed in self.hmxrb_seeds:        #if this loop takes too long look into vectorizing this  *biztodo
            t_hmxb = hmxb_times_rows[hmxb_seeds_rows == seed]
            t_caseA = caseA_times[caseA_seeds == seed]
            # skip if this system never had Case A
            if t_caseA.size == 0:
                continue

            if np.any(t_caseA < np.min(t_hmxb)):
                valid_hmxb_seeds.append(seed)

                visible_times[seed] = {"luminosity": {"time": np.array([]), "Lx": np.array([])}, "disc": {"time": np.array([]), "Lx": np.array([])}}

                #now including extra to be checks for when the high mass xray binary would be observable:
                print("high mass xray binary time for seed", seed ,"from:", np.min(t_hmxb), "to:", np.max(t_hmxb)) 

                # rows belonging to this seed AND HMXB phase
                seed_rows = (hmxb_seeds_rows == seed)

                t_hmxb = hmxb_times_rows[seed_rows]

                # get mass transfer rates for these rows
                mt_rates = mt_rates_all[hmxrb_mask][seed_rows]

                #now must convert this m_dot to luminosity estimation 
                #attaching units with astropy units (easy to convert later) 
                
                mt_rates = mt_rates * u.Msun/u.yr
                mt_rates = mt_rates.to(u.g/u.s) #may need to change this back surely itll be huge ??
                eta = 0.5 #set this in the yaml?
                # Lx = η Mdot c^2
                Lx = eta * mt_rates * const.c**2

                #converting to erg to compare to telescope values 
                Lx = Lx.to(u.erg / u.s)

                #nustar value 
                L_nustar = 1e32 * u.erg / u.s

                #XMM value (biztodo) - look this up 
                

                #detectable if the luminosity is above this value 
                detectable = Lx > L_nustar

                Lx_detectable = Lx[detectable]

                t_detectable = t_hmxb[detectable]
                if t_detectable.size > 0:
                    visible_times[seed]["luminosity"]["time"] = t_detectable
                    visible_times[seed]["luminosity"]["Lx"]   = Lx_detectable
                    print(
                        f"Observable HMXB (by luminosity) for seed {seed} from "
                        f"{np.min(t_detectable)} to {np.max(t_detectable)} Myr"
                    )


                #now following Hirai and Mandel paper - disc formation for observability

                #need R_circ > R_ISCO

                #will need to sort units here *** 

                sep  = sep_all[hmxrb_mask][seed_rows]  * u.Rsun
                rad1 = rad1_all[hmxrb_mask][seed_rows] * u.Rsun
                rad2 = rad2_all[hmxrb_mask][seed_rows] * u.Rsun

                m1 = mass1[hmxrb_mask][seed_rows] * u.Msun
                m2 = mass2[hmxrb_mask][seed_rows] * u.Msun

                st1_rows = st1[hmxrb_mask][seed_rows]
                st2_rows = st2[hmxrb_mask][seed_rows]

                #determining which object is the BH and which is the donor
                bh_is_1 = st1_rows >= 14
                bh_is_2 = st2_rows >= 14

                bh_m    = np.where(bh_is_1, m1, m2)
                donor_m = np.where(bh_is_1, m2, m1)
                donor_r = np.where(bh_is_1, rad2, rad1)
                a_orb   = sep
                
                # bh_m = #need to determine which object is the BH and which is therefore the donor
                # donor_m = # need to determine which object is the BH and which is the donor
                # a_orb = #separation
                # donor_r = #once know which is donor then get this radius from RLOF file  

                #v_orb = np.sqrt(G*(M_donor + M_bh)/a_orb)
                
                R_sch = (2 * const.G * bh_m) / const.c**2 
                R_isco = 3 * R_sch  #R Innermost Stable Circular Orbit

                orb_v = np.sqrt(const.G * (bh_m + donor_m) / a_orb )
                
                v_esc = np.sqrt(2 * const.G * donor_m / donor_r)
                v_inf = 2.6 * v_esc
                rad_wind_v = v_inf * (1 - donor_r/a_orb)

                
                wind_v = np.sqrt(orb_v**2 + rad_wind_v**2) #wind velocity in the vicinity of the accreting object (vector sum of the orbital velocyt and the radial wind velocity) 
                
                R_acc = (2 * const.G * bh_m) / wind_v**2 #capture radius  
                
                omega = np.sqrt(const.G*(donor_m + bh_m)/a_orb**3) #orbital angular frequency 

                # ang mom of captured wind (eta_j = 0.1), assumption in paper, biz check
                eta_j = 0.1
                
                angmom_j = eta_j * omega * R_acc**2
                
                R_circ = (angmom_j**2) / (const.G * bh_m) #R CIRC

                disc_forms = R_circ > R_isco

                Lx_disc_forms = Lx[disc_forms] #using luminosity caluclations defined above - could include additional luminsoity condition?

                t_disc_forms = t_hmxb[disc_forms]
                
                if t_disc_forms.size > 0:
                    visible_times[seed]["disc"]["time"] = t_disc_forms
                    visible_times[seed]["disc"]["Lx"]   = Lx_disc_forms
                    print(
                        f"Observable HMXB (disc formation) for seed {seed} from "
                        f"{np.min(t_disc_forms)} to {np.max(t_disc_forms)} Myr"
                    )
                
                
        #should only be the seeds which are case A then hmxrb
        self.hmxrb_seeds = np.array(valid_hmxb_seeds)
    
        #print(f"Identified {self.hmxrb_seeds.size} HMXRB systems.")
        #if self.hmxrb_seeds.size > 0:
            #seed_lines = "\n".join(f"--random-seed {seed}" for seed in self.hmxrb_seeds.tolist())
            #print(seed_lines)

        return self.hmxrb_seeds, visible_times

    def close(self):
        self.h5file.close()




def get_bbh_hmxrb_overlap(data_path):

    # --- run BBH finder ---
    bbh = determine_bbh(data_path)
    bbh_seeds = bbh.identify_bbh_systems()
    bbh.close()

    # --- run HMXRB finder ---
    hmx = determine_hmxrb(data_path)
    hmxrb_seeds, visible_times = hmx.identify_hmxrb_systems()
    hmx.close()

    # --- intersect seeds ---
    overlapping_seeds = np.intersect1d(bbh_seeds, hmxrb_seeds)

    print(f"\nFound {overlapping_seeds.size} systems that are BOTH BBH and HMXRB.")
    if overlapping_seeds.size > 0:
            seed_lines = "\n".join(f"--random-seed {seed}" for seed in overlapping_seeds.tolist())
            print(seed_lines)

    return overlapping_seeds, visible_times

    



class spin_param_calculations(object):
    """
    This class calculated the time evolution of the dimensionless spin parameter 'a' for black holes  
    """

    def __init__(self, data_path): #could add output directory here too later, and things about spin models etc if needed. data path should be as worked out above
        
        self.path = Path(data_path)

        if not self.path.is_file():
            raise ValueError(f"h5 file not found. Wrong path given?\npath checked = {self.path}")
        self.h5file = h5.File(self.path, "r")

    def calculate_dimensionless_spin(self, Mass, Ang_Mom):
        # Calculate the dimensionless spin parameter chi = cJ/GM^2 
    
        M = Mass * u.Msun
        J = Ang_Mom * u.Msun * u.au**2 / u.yr

        G = const.G  # [gr cm s^2]
        c = const.c    # [cm/s]

        chi = (c * J) / (G * M**2)

        return chi.decompose().value #dimensionless value only 
        
        
    def setData(self):
        # Reads in params needed from the detailed h5 file, then calculated the dimensionless spin parameter 

        Detailed_BSE = self.h5file #h5 file (only one)

        # New masks 
        record_type = Detailed_BSE['Record_Type'][()]
        mask_timestep = record_type == 4

        stellar_type_1 = Detailed_BSE['Stellar_Type(1)'][()]
        stellar_type_2 = Detailed_BSE['Stellar_Type(2)'][()]

        mask_bh1 = mask_timestep & (stellar_type_1 == 14)
        mask_bh2 = mask_timestep & (stellar_type_2 == 14)
        mask_ns1 = mask_timestep & (stellar_type_1 == 13) #Introducing masking for neutron stars also. (same dimensionless spin param)
        mask_ns2 = mask_timestep & (stellar_type_2 == 13)

        mask_1 = mask_bh1 | mask_ns1
        mask_2 = mask_bh2 | mask_ns2

        #assigning params below

        Mass1 = Detailed_BSE['Mass(1)'][()][mask_1] #Msol
        Ang_Mom1 = Detailed_BSE['Ang_Momentum(1)'][()][mask_1] #Msol AU^2 yr^-1
        time1 = Detailed_BSE['Time'][()][mask_1] * u.Myr

        Mass2 = Detailed_BSE['Mass(2)'][()][mask_2] #Msol
        Ang_Mom2 = Detailed_BSE['Ang_Momentum(2)'][()][mask_2] #Msol AU^2 yr^-1
        time2 = Detailed_BSE['Time'][()][mask_2] * u.Myr #Absolute binary evolution time since ZAMS - rather than using age which resets to 0 at each new phase, this is the time that continues to increase throughout the binary evolution.

        chi1 = self.calculate_dimensionless_spin(Mass1, Ang_Mom1)
        chi2 = self.calculate_dimensionless_spin(Mass2, Ang_Mom2)
        

        return chi1, chi2, time1, time2

    def close(self):
        self.h5file.close() 


#removing this function from the calculations object for different directory input (can re-introduce later)
def process_all_detailed_files(detailed_output_dir):
        """
        Loop through all HDF5 files in the Detailed_Output directory (takes detailed output directory as input param)
        and calculate spin evolution for each.
        Returns a dictionary with SEED as key and spin data as value.
        """

        results = {}
        h5_files = detailed_output_dir.glob("BSE_Detailed_Output_*.h5")

        for h5path in h5_files:
            try:
                with h5.File(h5path, "r") as f:
                    seed = int(f["SEED"][0])
                    calc = spin_param_calculations(h5path)
                    chi1, chi2, time1, time2 = calc.setData()
                    results[seed] = {
                        'chi1': chi1,
                        'chi2': chi2,
                        'time1': time1,
                        'time2': time2
                    }
            except (OSError, KeyError, ValueError) as e:
                print(f"Skipping file {h5path} due to error: {e}")
                continue
        print(f"Processed {len(results)} BBH detailed output files.")
        print(results)
        return results


        
def plot_spin_evolution(time, spin1, spin2, outdir='.', show=True, use_latex=True):
    """
    Plots the time evolution of the dimensionless spin parameters for both components in a BBH system.
    """

    if use_latex:
        rcParams['text.usetex'] = True
        rcParams['font.family'] = 'serif'
        rcParams['font.serif'] = ['Computer Modern Roman']
        rcParams['font.size'] = 14

    plt.figure(figsize=(10,6))
    plt.plot(time, spin1, label='Spin Parameter a1', color='blue')
    plt.plot(time, spin2, label='Spin Parameter a2', color='red')
    plt.xlabel('Time (Myr)')
    plt.ylabel('Dimensionless Spin Parameter (a)')
    plt.title('Evolution of Dimensionless Spin Parameters in BBH System')
    plt.legend()
    plt.grid(True)

    if outdir:
        Path(outdir).mkdir(parents=True, exist_ok=True)
        plt.savefig(os.path.join(outdir, 'spin_evolution.png'), dpi=300)

    if show:
        plt.show()
    else:
        plt.close()


def plot_all_spin_evolution(results_dict, outdir='.', show=True, use_latex=False):
    """
    Plots the time evolution of dimensionless spin parameters for all BBH systems in the results dictionary.
    
    Parameters:
    -----------
    results_dict : dict
        Dictionary returned from process_all_detailed_files() with SEED as keys
    outdir : str
        Output directory for saving plots
    show : bool
        Whether to display plots interactively
    use_latex : bool
        Whether to use LaTeX for text rendering
    """
    
    if use_latex:
        rcParams['text.usetex'] = True
        rcParams['font.family'] = 'serif'
        rcParams['font.serif'] = ['Computer Modern Roman']
        rcParams['font.size'] = 14

    if not results_dict:
        print("No results to plot")
        return
    
    # Create individual plots for each SEED
    for seed, data in results_dict.items():
        fig, ax = plt.subplots(figsize=(5, 3))
        
        ax.scatter(data['time1'].value, data['chi1'], label='Spin Parameter $\chi_1$', color='blue', linewidth=2)
        ax.scatter(data['time2'].value, data['chi2'], label='Spin Parameter $\chi_2$', color='red', linewidth=2)
        
        ax.set_xlabel('Time (Myr)', fontsize=9)
        ax.set_ylabel('Dimensionless Spin Parameter (a)', fontsize=9)
        ax.set_title(f'Spin Evolution: SEED {seed}', fontsize=10)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        
        # Save individual plot
        if outdir:
            Path(outdir).mkdir(parents=True, exist_ok=True)
            outfile = os.path.join(outdir, f'spin_evolution_SEED_{seed}.png')
            plt.savefig(outfile, dpi=300, bbox_inches='tight')
            print(f"Saved: {outfile}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    # Optional: Create a combined comparison plot
    fig, axes = plt.subplots(len(results_dict), 1, figsize=(6, 2*len(results_dict)))
    
    # Handle single vs multiple subplots
    if len(results_dict) == 1:
        axes = [axes]
    
    for idx, (seed, data) in enumerate(results_dict.items()):
        ax = axes[idx]
        ax.scatter(data['time1'], data['chi1'], label='χ₁', color='blue', linewidth=2)
        ax.scatter(data['time2'], data['chi2'], label='χ₂', color='red', linewidth=2)
        ax.set_ylabel('Dimensionless Spin (a)', fontsize=9)
        ax.set_title(f'SEED {seed}', fontsize=10)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        
        if idx == len(results_dict) - 1:
            ax.set_xlabel('Time (Myr)', fontsize=10)
    
    plt.tight_layout()
    
    if outdir:
        Path(outdir).mkdir(parents=True, exist_ok=True)
        outfile = os.path.join(outdir, 'spin_evolution_all_systems.png')
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        print(f"Saved combined plot: {outfile}")
    
    if show:
        plt.show()
    else:
        plt.close()


#####################################################################

class SpinAtEventsCalculations:
    #calculates spin at evolutionary timesteps for a single detailed output file. 

    def __init__(self, data_path):
        self.path = Path(data_path)

        if not self.path.is_file():
            raise ValueError(f"H5 file not found:\n{self.path}")

        self.h5file = h5.File(self.path, "r")

    def detect_evolutionary_events(self):

        time_data = self.h5file['Time'][()]
        mt_history = self.h5file['MT_History'][()]
        stype1 = self.h5file['Stellar_Type(1)'][()]
        stype2 = self.h5file['Stellar_Type(2)'][()]

        event_indices = [0]
        event_types = ['ZAMS']

        for ii in range(1, len(time_data)):

            if (mt_history[ii] > 0) and (mt_history[ii] != mt_history[ii-1]):
                event_indices.append(ii)
                event_types.append('MT')

            if stype1[ii] != stype1[ii-1]:
                event_indices.append(ii)
                event_types.append('SN_1' if stype1[ii] in [13,14] else 'Stype_1')

            if stype2[ii] != stype2[ii-1]:
                event_indices.append(ii)
                event_types.append('SN_2' if stype2[ii] in [13,14] else 'Stype_2')

        event_indices.append(-1)
        event_types.append('End')

        return event_indices, event_types

    def calculate_spin_at_events(self):

        spin_calc = spin_param_calculations(self.path)

        time_data = self.h5file['Time'][()]

        indices, _ = self.detect_evolutionary_events()

        indices = np.array([idx if idx >= 0 else len(time_data)+idx for idx in indices])
        event_times = time_data[indices]

        mass1 = self.h5file['Mass(1)'][()][indices]
        mass2 = self.h5file['Mass(2)'][()][indices]
        ang1 = self.h5file['Ang_Momentum(1)'][()][indices]
        ang2 = self.h5file['Ang_Momentum(2)'][()][indices]
        st1 = self.h5file['Stellar_Type(1)'][()][indices]
        st2 = self.h5file['Stellar_Type(2)'][()][indices]

        chi1 = np.full(len(indices), np.nan)
        chi2 = np.full(len(indices), np.nan)

        for i in range(len(indices)):
            if st1[i] in [13,14]:
                chi1[i] = spin_calc.calculate_dimensionless_spin(mass1[i], ang1[i])

            if st2[i] in [13,14]:
                chi2[i] = spin_calc.calculate_dimensionless_spin(mass2[i], ang2[i])

        return event_times, chi1, chi2

    def close(self):
        self.h5file.close()


##batch processing 

def process_event_spin_for_all(detailed_files_dir):
    #pass the directory of the detailed output files. Loop through each file, calculate the spin at events, and store results in a dictionary with SEED as key.

    results = {}
    h5_files = detailed_files_dir.glob("BSE_Detailed_Output_*.h5")

    for h5path in h5_files:
        try:
            with h5.File(h5path, "r") as f:
                seed = int(f["SEED"][0])
                calc = SpinAtEventsCalculations(h5path)
                times, chi1, chi2 = calc.calculate_spin_at_events()
                results[seed] = {
                    'times': times,
                    'chi1': chi1,
                    'chi2': chi2
                }

        except (OSError, KeyError, ValueError) as e:
            print(f"Skipping file {h5path} due to error: {e}")
            continue

    print(f"Processed {len(results)} BBH detailed output files.")
    return results


#plotting 

def plot_spin_at_events(event_results, outdir='.', show=True, use_latex=False):

    if use_latex:
        rcParams['text.usetex'] = True
        rcParams['font.family'] = 'serif'
        rcParams['font.serif'] = ['Computer Modern Roman']
        rcParams['font.size'] = 14

    if not event_results:
        print("No results to plot")
        return

    #plots for individual seeds 
    for seed, data in event_results.items():
        fig, ax = plt.subplots(figsize=(5, 3)) #making them smaller for notebook

        ax.scatter(data['times'], data['chi1'], label='Spin Parameter $\chi_1$', color='blue')
        ax.scatter(data['times'], data['chi2'], label='Spin Parameter $\chi_2$', color='red')

        ax.set_xlabel('Time (Myr)', fontsize=8)
        ax.set_ylabel('Dimensionless Spin Parameter (a)', fontsize=8)
        ax.set_title(f'Spin Evolution: SEED {seed}', fontsize=10)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

        if outdir:
            Path(outdir).mkdir(parents=True, exist_ok=True)
            outfile = os.path.join(outdir, f'spin_at_events_SEED_{seed}.png')
            plt.savefig(outfile, dpi=300, bbox_inches='tight')
            print(f"Saved: {outfile}")

        if show:
            plt.show()

        else:
            plt.close()

    #combined plot for comparison 
    fig, axes = plt.subplots(len(event_results), 1, figsize=(6, 2*len(event_results)))
    if len(event_results) == 1:
        axes = [axes]

    for idx, (seed, data) in enumerate(event_results.items()):
        ax = axes[idx]
        ax.scatter(data['times'], data['chi1'], label='χ₁', color='blue')
        ax.scatter(data['times'], data['chi2'], label='χ₂', color='red')
        ax.set_ylabel('Dimensionless Spin (a)', fontsize=8)
        ax.set_title(f'SEED {seed}', fontsize=10)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

        if idx == len(event_results) - 1:
            ax.set_xlabel('Time (Myr)', fontsize=8)
    plt.tight_layout()

    if outdir:
        Path(outdir).mkdir(parents=True, exist_ok=True)
        outfile = os.path.join(outdir, 'spin_at_events_all_systems.png')
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        print(f"Saved combined plot: {outfile}")

    if show:     
        plt.show()

    else:
        plt.close()


#summary plots for poster 

#####################################################################
# New class for plotting BH vs donor masses for BBH systems
#####################################################################

class BHDonorMassPlot:
    """
    Class to extract black hole masses and donor masses at the time of BH formation
    from COMPAS detailed output files and create scatter plots.
    """

    def __init__(self, detailed_file_path):
        self.path = Path(detailed_file_path)
        if not self.path.is_file():
            raise ValueError(f"H5 file not found:\n{self.path}")
        self.h5file = h5.File(self.path, "r")

    def get_bh_and_donor_masses(self):
        """
        Returns two lists: BH masses and the corresponding donor masses at the time of BH formation.
        Only considers the first BH formation in each binary.
        """

        mass1 = self.h5file['Mass(1)'][()]
        mass2 = self.h5file['Mass(2)'][()]
        stype1 = self.h5file['Stellar_Type(1)'][()]
        stype2 = self.h5file['Stellar_Type(2)'][()]

        record_type = self.h5file['Record_Type'][()]
        mask_timestep = record_type == 4  # only physical timesteps

        # masks for BH formation
        mask_bh1 = mask_timestep & (stype1 == 14)
        mask_bh2 = mask_timestep & (stype2 == 14)

        # find first BH formation step for each component
        bh_masses = []
        donor_masses = []

        if np.any(mask_bh1):
            idx = np.where(mask_bh1)[0][0]
            bh_masses.append(mass1[idx])
            donor_masses.append(mass2[idx])

        if np.any(mask_bh2):
            idx = np.where(mask_bh2)[0][0]
            bh_masses.append(mass2[idx])
            donor_masses.append(mass1[idx])

        return bh_masses, donor_masses

    def close(self):
        self.h5file.close()


# Batch processing for multiple BBH detailed files
def process_all_bh_donor_masses(detailed_files_dir):
    """
    Loops through all BBH detailed output files in the given directory.
    Returns a dictionary with SEED as key and {'bh': [], 'donor': []} as value.
    """

    results = {}
    h5_files = detailed_files_dir.glob("BSE_Detailed_Output_*.h5")

    for h5path in h5_files:
        try:
            with h5.File(h5path, "r") as f:
                seed = int(f['SEED'][0])
                plotter = BHDonorMassPlot(h5path)
                bh_mass, donor_mass = plotter.get_bh_and_donor_masses()
                results[seed] = {'bh': bh_mass, 'donor': donor_mass}
                plotter.close()
        except (OSError, KeyError, ValueError) as e:
            print(f"Skipping {h5path} due to error: {e}")
            continue

    print(f"Processed {len(results)} BBH detailed output files for BH vs donor masses.")
    return results


# Scatter plotting function
def plot_bh_vs_donor_mass(results_dict, outdir='.', show=True, use_latex=False):
    """
    Creates a scatter plot of donor mass vs BH mass for all BBH systems.
    """

    if use_latex:
        rcParams['text.usetex'] = True
        rcParams['font.family'] = 'serif'
        rcParams['font.serif'] = ['Computer Modern Roman']
        rcParams['font.size'] = 14

    all_bh = []
    all_donor = []

    for seed, data in results_dict.items():
        all_bh.extend(data['bh'])
        all_donor.extend(data['donor'])

    plt.figure(figsize=(8,6))
    plt.scatter(all_donor, all_bh, alpha=0.7)
    plt.xlabel('Donor Mass at BH Formation [M$_\odot$]')
    plt.ylabel('Black Hole Mass [M$_\odot$]')
    plt.title('Scatter of Donor Mass vs BH Mass')
    plt.grid(True)

    if outdir:
        Path(outdir).mkdir(parents=True, exist_ok=True)
        outfile = os.path.join(outdir, 'bh_vs_donor_mass.png')
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        print(f"Saved scatter plot: {outfile}")

    if show:
        plt.show()
    else:
        plt.close()

class BHDonorMassSpinPlot(BHDonorMassPlot):
    """
    Extension of BHDonorMassPlot to include dimensionless spin at BH formation
    and create a two-panel plot: 
    (1) Donor mass vs BH mass
    (2) Spin at BH formation vs BH mass
    """

    def get_bh_donor_and_spin(self):
        """
        Returns BH masses, donor masses, and spin of BHs at formation.
        """

        mass1 = self.h5file['Mass(1)'][()]
        mass2 = self.h5file['Mass(2)'][()]
        ang1 = self.h5file['Ang_Momentum(1)'][()]
        ang2 = self.h5file['Ang_Momentum(2)'][()]
        stype1 = self.h5file['Stellar_Type(1)'][()]
        stype2 = self.h5file['Stellar_Type(2)'][()]

        record_type = self.h5file['Record_Type'][()]
        mask_timestep = record_type == 4  # physical timesteps

        mask_bh1 = mask_timestep & (stype1 == 14)
        mask_bh2 = mask_timestep & (stype2 == 14)

        bh_masses = []
        donor_masses = []
        bh_spins = []

        spin_calc = spin_param_calculations(self.path)

        if np.any(mask_bh1):
            idx = np.where(mask_bh1)[0][0]
            bh_masses.append(mass1[idx])
            donor_masses.append(mass2[idx])
            bh_spins.append(spin_calc.calculate_dimensionless_spin(mass1[idx], ang1[idx]))

        if np.any(mask_bh2):
            idx = np.where(mask_bh2)[0][0]
            bh_masses.append(mass2[idx])
            donor_masses.append(mass1[idx])
            bh_spins.append(spin_calc.calculate_dimensionless_spin(mass2[idx], ang2[idx]))

        return bh_masses, donor_masses, bh_spins


def process_all_bh_donor_and_spin(detailed_files_dir):
    """
    Batch processing for multiple BBH detailed files.
    Returns dictionary with SEED as key and dict of BH mass, donor mass, and spin.
    """
    results = {}
    h5_files = detailed_files_dir.glob("BSE_Detailed_Output_*.h5")

    for h5path in h5_files:
        try:
            with h5.File(h5path, "r") as f:
                seed = int(f['SEED'][0])
                plotter = BHDonorMassSpinPlot(h5path)
                bh_mass, donor_mass, bh_spin = plotter.get_bh_donor_and_spin()
                results[seed] = {'bh': bh_mass, 'donor': donor_mass, 'spin': bh_spin}
                plotter.close()
        except (OSError, KeyError, ValueError) as e:
            print(f"Skipping {h5path} due to error: {e}")
            continue

    print(f"Processed {len(results)} BBH detailed output files for BH/donor mass and spin.")
    return results


def plot_bh_donor_and_spin(results_dict, outdir='.', show=True, use_latex=False):
    """
    Creates a two-panel plot:
    Top: Donor mass vs BH mass
    Bottom: BH spin at formation vs BH mass
    """

    if use_latex:
        rcParams['text.usetex'] = True
        rcParams['font.family'] = 'serif'
        rcParams['font.serif'] = ['Computer Modern Roman']
        rcParams['font.size'] = 14

    all_bh = []
    all_donor = []
    all_spin = []

    for seed, data in results_dict.items():
        all_bh.extend(data['bh'])
        all_donor.extend(data['donor'])
        all_spin.extend(data['spin'])

    fig, axes = plt.subplots(2, 1, figsize=(8, 10), sharex=True)

    # Panel 1: donor vs BH mass
    axes[0].scatter(all_donor, all_bh, alpha=0.7)
    axes[0].set_ylabel('Black Hole Mass [M$_\odot$]')
    axes[0].set_title('Donor Mass vs BH Mass')
    axes[0].grid(True)

    # Panel 2: spin vs BH mass
    axes[1].scatter(all_bh, all_spin, alpha=0.7, color='orange')
    axes[1].set_xlabel('Black Hole Mass [M$_\odot$]')
    axes[1].set_ylabel('Dimensionless Spin (a)')
    axes[1].set_title('BH Spin at Formation vs BH Mass')
    axes[1].grid(True)

    plt.tight_layout()

    if outdir:
        Path(outdir).mkdir(parents=True, exist_ok=True)
        outfile = os.path.join(outdir, 'bh_donor_mass_and_spin.png')
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        print(f"Saved combined two-panel plot: {outfile}")

    if show:
        plt.show()
    else:
        plt.close()


# plotting for all simulations 

# will take a list of paths to the detailed output files for each simulation 

# dark baclground plotting: for poster

def plot_multiple_runs_from_dirs(run_dirs, outdir='.', show=True, use_latex=False):
    """
    Plot multiple COMPAS runs on the same 2-panel figure.
    Each run is plotted with a different color and labeled by its directory name.
    """

    #plt.style.use('dark_background')

    if use_latex:
        rcParams['text.usetex'] = True
        rcParams['font.family'] = 'serif'
        rcParams['font.serif'] = ['Computer Modern Roman']
        rcParams['font.size'] = 14

    colours = plt.get_cmap("Dark2").colors

    fig, axes = plt.subplots(2, 1, figsize=(6, 12))

    for i, run_dir in enumerate(run_dirs):

        colour = colours[i % len(colours)]

        label = run_dir.parent.parent.name

        # process the run
        results_dict = process_all_bh_donor_and_spin(run_dir)

        all_bh = []
        all_donor = []
        all_spin = []

        for seed, data in results_dict.items():
            all_bh.extend(data['bh'])
            all_donor.extend(data['donor'])
            all_spin.extend(data['spin'])

        # Panel 1
        axes[0].scatter(
            all_donor,
            all_bh,
            color=colour,
            alpha=0.8,
            s=30,
            linewidth=0,
            label=label
        )

        # Panel 2
        axes[1].scatter(
            all_bh,
            all_spin,
            color=colour,
            alpha=0.8,
            s=30,
            linewidth=0,
            label=label
        )

    # ---- formatting ----
    axes[0].set_ylabel('Black Hole Mass [M$_\\odot$]')
    axes[0].set_xlabel('Donor Mass at BH Formation [M$_\\odot$]')
    axes[0].set_title('Donor Mass vs BH Mass')
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    axes[1].set_xlabel('Compact Object Mass [M$_\\odot$]')
    axes[1].set_ylabel('Dimensionless Spin (a)')
    axes[1].set_title('Spin at Formation vs Mass')
    axes[1].set_yscale("symlog", linthresh=1e-6)
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()

    plt.tight_layout()

    if outdir:
        Path(outdir).mkdir(parents=True, exist_ok=True)
        outfile = os.path.join(outdir, 'bh_donor_spin_multiple_runs.png')
        plt.savefig(outfile, dpi=300, bbox_inches='tight')
        print(f"Saved combined plot: {outfile}")

    if show:
        plt.show()
    else:
        plt.close()