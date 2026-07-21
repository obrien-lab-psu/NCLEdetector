#!/usr/bin/env python3
import sys, getopt, math, os, time, traceback, glob, multiprocessing, copy
import numpy as np
from scipy.stats import norm
import pandas as pd
from statsmodels.stats.multitest import multipletests
import parmed as pmd
import mdtraj as mdt
import pathlib, re
import logging
from EntDetect._logging import setup_logger

class MassSpec:
    """
    Compare ensembles of protein structures with LiPMS and XLMS experimental data.
    Primary author: Yang Jiang
    Secondary author: Ian Sitarik
    """
    #############################################################################################################
    def __init__(self, msm_data_file:str, meta_dist_file:str, LiPMS_exp_file:str, sasa_data_file:str, XLMS_exp_file:str, dist_data_file:str,
                 cluster_data_file:str, OPpath:str, AAdcd_dir:str, native_AA_pdb:str, native_state_idx:int, state_idx_list:list, prot_len:int, n_analysis_frames:int,
                 rm_traj_list:list=[], outdir:str='./', ID:str='', xp_dir:str=None, resid2residueidx_map:dict={},
                 sasa_dir:str=None, n_traj:int=None, sasa_xp_frames_per_traj:int=None, collect_jwalk_npy:bool=False,
                 start:int=0, end:int=999999999999, stride:int=1, verbose:bool=False, num_perm:int=10000, n_boot:int=10000, lag_frame:int=1, nproc:int=1, log_level:int=logging.INFO, logdir:str=None):


        self.msm_data_file = msm_data_file
        self.meta_dist_file = meta_dist_file
        self.LiPMS_exp_file = LiPMS_exp_file
        self.sasa_data_file = sasa_data_file
        self.XLMS_exp_file = XLMS_exp_file
        self.dist_data_file = dist_data_file
        self.xp_dir = xp_dir
        self.cluster_data_file = cluster_data_file
        self.OPpath = OPpath
        self.AAdcd_dir = AAdcd_dir
        self.native_AA_pdb = native_AA_pdb
        self.native_state_idx = native_state_idx
        self.state_idx_list = state_idx_list
        self.rm_traj_list = rm_traj_list
        self.outdir = outdir
        self.ID = ID

        # Initialize logging before any log calls in constructor.
        self.logger = setup_logger('MassSpec', outdir=logdir if logdir is not None else self.outdir, ID=self.ID, log_level=log_level)

        self.resid2residueidx_map = resid2residueidx_map
        if len(self.resid2residueidx_map) == 0:
            self.resid2residueidx_map = {i + 1:i for i in range(prot_len)}
            self.logger.info(f'No resid2residueidx_map provided, using identity mapping for protein length {prot_len} with an offset of -1')
        self.start = start
        self.end = end
        self.stride = stride
        self.verbose = verbose

        self.res_buffer = 5
        self.num_perm = num_perm
        self.prot_len = prot_len
        self.nproc = nproc

        self.if_calc_M = 1

        self.n_analysis_frames = n_analysis_frames # last X ns actually analyzed (trailing-frame window shared by MSM/G/Q/SASA/XP)
        self.lag_frame = lag_frame # down sample trajectories at each #lag_frame frame
        self.n_boot = n_boot
        #self.n_boot = 100

        # make the outdir if it doesnt existgs
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)
            self.logger.debug(f'Creating directory: {self.outdir}')

        # Optionally collect per-trajectory OP files into the dense arrays
        # (SASA.npy / Jwalk.npy) that the consistency test consumes.
        self._maybe_collect_op_arrays(sasa_dir=sasa_dir, n_traj=n_traj,
                                      sasa_xp_frames_per_traj=sasa_xp_frames_per_traj, collect_jwalk_npy=collect_jwalk_npy)
    ##############################################################################

    ##############################################################################
    def _maybe_collect_op_arrays(self, sasa_dir, n_traj, sasa_xp_frames_per_traj, collect_jwalk_npy):
        """Build SASA.npy (and optionally Jwalk.npy) from per-trajectory OP files.

        Runs only when the corresponding cached array is not already available.
        SASA is required as a dense array; XL-MS scoring can stream from
        ``xp_dir`` directly, so Jwalk.npy is only built when explicitly requested
        via ``collect_jwalk_npy``.
        """
        need_sasa = self.sasa_data_file is None or not os.path.exists(self.sasa_data_file)
        if need_sasa and sasa_dir is not None:
            cached = os.path.join(self.outdir, 'SASA.npy')
            if os.path.exists(cached):
                self.sasa_data_file = cached
                self.logger.info(f'Using cached SASA array: {cached}')
            else:
                if n_traj is None or sasa_xp_frames_per_traj is None:
                    raise ValueError('n_traj and sasa_xp_frames_per_traj are required to collect SASA from sasa_dir')
                from EntDetect.order_params import CollectOP
                self.logger.info(f'Collecting SASA from {sasa_dir} into {self.outdir}/SASA.npy')
                collector = CollectOP(sasa_dir=sasa_dir, xp_dir=self.xp_dir, outdir=self.outdir,
                                      ID=self.ID, n_traj=n_traj, sasa_xp_frames_per_traj=sasa_xp_frames_per_traj, prot_len=self.prot_len)
                self.sasa_data_file = collector.collect_SASA()

        if collect_jwalk_npy:
            need_dist = self.dist_data_file is None or not os.path.exists(self.dist_data_file)
            if need_dist and self.xp_dir is not None:
                cached = os.path.join(self.outdir, 'Jwalk.npy')
                if os.path.exists(cached):
                    self.dist_data_file = cached
                    self.logger.info(f'Using cached Jwalk array: {cached}')
                else:
                    if n_traj is None or sasa_xp_frames_per_traj is None:
                        raise ValueError('n_traj and sasa_xp_frames_per_traj are required to collect Jwalk from xp_dir')
                    from EntDetect.order_params import CollectOP
                    self.logger.info(f'Collecting Jwalk from {self.xp_dir} into {self.outdir}/Jwalk.npy')
                    collector = CollectOP(sasa_dir=sasa_dir, xp_dir=self.xp_dir, outdir=self.outdir,
                                          ID=self.ID, n_traj=n_traj, sasa_xp_frames_per_traj=sasa_xp_frames_per_traj, prot_len=self.prot_len)
                    self.dist_data_file = collector.collect_Jwalk()
    ##############################################################################

    ##############################################################################
    def load_LiPMS_data(self, file_path):
        df = pd.read_excel(file_path)
        LiPMS_sig_data = {}
        for index, row in df.iterrows():
            peptide = row['Cut Site']
            if '-' not in peptide:
                resid = int(peptide.strip()[1:])-1
                peptide_range = set(list(np.arange(np.max([0,resid-self.res_buffer]), np.min([self.prot_len,resid+self.res_buffer+1]))))
                if peptide not in LiPMS_sig_data:
                    LiPMS_sig_data[peptide] = {}
                    LiPMS_sig_data[peptide]['peptide_range'] = []
                    LiPMS_sig_data[peptide]['qual_change'] = []
                LiPMS_sig_data[peptide]['peptide_range'] = peptide_range
                if row['Log2 FC'] < 0:
                    LiPMS_sig_data[peptide]['qual_change'].append(-1)
                elif row['Log2 FC'] > 0:
                    LiPMS_sig_data[peptide]['qual_change'].append(1)
                elif row['Log2 FC'] == 0:
                    LiPMS_sig_data[peptide]['qual_change'].append(0)
        removed_key_list = []
        for k, v in LiPMS_sig_data.items():
            a = list(set(v['qual_change']))
            if len(a) == 1:
                v['qual_change'] = a[0]
            else:
                self.logger.info('Site %s has inconsistent changes in abundance: %s'%(k, str(v['qual_change'])))
                removed_key_list.append(k)
        for k in removed_key_list:
            LiPMS_sig_data.pop(k)
        LiPMS_sig_data = dict(sorted(LiPMS_sig_data.items(), key=lambda item: int(item[0].strip()[1:])))
        return LiPMS_sig_data
    ##############################################################################

    ##############################################################################        
    def load_XLMS_data(self, file_path):
        df = pd.read_excel(file_path)
        XLMS_sig_data = {}
        for index, row in df.iterrows():
            XL_site_key = row['Pairs'].strip()
            XL_site = XL_site_key.split('-')
            if XL_site_key not in XLMS_sig_data.keys():
                XLMS_sig_data[XL_site_key] = {}
                XLMS_sig_data[XL_site_key]['qual_change'] = []
                XLMS_sig_data[XL_site_key]['res_pair'] = []
            XLMS_sig_data[XL_site_key]['res_pair'] = [int(XL_site[0][1:])-1, int(XL_site[1][1:])-1]
            if row['log2(heavy/light)'] < 0:
                XLMS_sig_data[XL_site_key]['qual_change'].append(-1)
            elif row['log2(heavy/light)'] > 0:
                XLMS_sig_data[XL_site_key]['qual_change'].append(1)
            elif row['log2(heavy/light)'] == 0:
                XLMS_sig_data[XL_site_key]['qual_change'].append(0)
        removed_key_list = []
        for k, v in XLMS_sig_data.items():
            a = list(set(v['qual_change']))
            if len(a) == 1:
                v['qual_change'] = a[0]
            else:
                self.logger.info('Sites %s has inconsistent changes in abundance: %s'%(k, str(v['qual_change'])))
                removed_key_list.append(k)
        for k in removed_key_list:
            XLMS_sig_data.pop(k)
        XLMS_sig_data = dict(sorted(XLMS_sig_data.items(), key=lambda item: (int(item[0].split('-')[0][1:]), int(item[0].split('-')[1][1:]))))
        return XLMS_sig_data
    ##############################################################################

    ##############################################################################        
    def score_XL(self, pair_AA, JWalk_dist):
        XL_offset = 1.1
        sc_length = {'K': 6.3,
                    'S': 2.5,
                    'T': 2.5,
                    'Y': 6.5,
                    'M': 1.5,}
        KK_mu = 18.6
        KK_sigma = 6.0
        KK_threshold = 33
        
        KK_mu += XL_offset
        KK_sigma = (XL_offset + 3*KK_sigma) / 3
        KK_threshold += XL_offset
        
        mu = KK_mu + (sc_length[pair_AA[0]] + sc_length[pair_AA[1]]) - 2*sc_length['K']
        sigma = (mu - (KK_mu - 3*KK_sigma)) / 3
        threshold = KK_threshold + mu - KK_mu
        
        N = norm(mu, sigma)
        
        if JWalk_dist == -1:
            score = 0
        elif JWalk_dist <= threshold:
            score = N.pdf(JWalk_dist)
        else:
            score = 0
        return score
    ##############################################################################

    ##############################################################################
    def perm_fun(self, perm_idx_list, combined_data, length_1):
        d_1 = combined_data[perm_idx_list[:length_1]]
        d_2 = combined_data[perm_idx_list[length_1:]]
        return self.statistic_fun(d_1, d_2, 0)
    ##############################################################################

    ##############################################################################        
    def permutation_test(self, perm_stat_fun, data_1, data_2, num_perm, side='!='):
        if side not in ['!=', '>', '<']:
            self.logger.info('side parameter is wrong for function "permutation_test". It must be "!=", ">", or "<".')
            sys.exit()
        combined_data = np.array(list(data_1) + list(data_2))
        perm_idx_list_0 = np.arange(len(combined_data))
        t0 = self.perm_fun(perm_idx_list_0, combined_data, len(data_1))
        
        pool = multiprocessing.Pool(self.nproc)
        pool_list = []
        start_time = time.time()
        self.logger.debug('start permutation test')
        for i in range(num_perm):
            perm_idx_list = np.random.permutation(perm_idx_list_0)
            pool_list.append(pool.apply_async(self.perm_fun, (perm_idx_list, combined_data, len(data_1),)))
        pool.close()
        pool.join()
        t_dist = [p.get() for p in pool_list]
        p = 0
        for t in t_dist:
            if side == '!=' and np.abs(t) >= np.abs(t0):
                p += 1
            elif side == '>' and t >= t0:
                p += 1
            elif side == '<' and t <= t0:
                p += 1
        p = (p+1)/(num_perm+1)
        used_time = time.time() - start_time
        self.logger.info('%.2fs'%used_time)
        return p
    ##############################################################################

    ##############################################################################
    def bootstrap(self, boot_fun, data, n_time):
        def fun(boot_fun, sample_idx_list):
            if len(data.shape) == 1:
                new_data = data[sample_idx_list]
                return boot_fun(new_data)
            elif len(data.shape) == 2:
                new_data = data[sample_idx_list, :]
                result = np.array([boot_fun(new_data[:,j]) for j in range(data.shape[1])])
                return result
        idx_list = np.arange(len(data))
        if len(data.shape) == 1:
            boot_stat = np.zeros(n_time)
        elif len(data.shape) == 2:
            boot_stat = np.zeros((n_time, data.shape[1]))
        else:
            self.logger.info('bootstrap: Can only handle 1 or 2 dimentional data')
            sys.exit()
        
        boot_stat = []
        for i in range(n_time):
            sample_idx_list = np.random.choice(idx_list, len(idx_list))
            bs = fun(boot_fun, sample_idx_list)
            boot_stat.append(bs)
        boot_stat = np.array(boot_stat)
        
        # pool = multiprocessing.Pool(nproc)
        # pool_list = []
        # start_time = time.time()
        # print('start bootstrapping')
        # for i in range(n_time):
            # sample_idx_list = np.random.choice(idx_list, len(idx_list))
            # pool_list.append(pool.apply_async(fun, (boot_fun, sample_idx_list)))
        # pool.close()
        # pool.join()
        # boot_stat = np.array([p.get() for p in pool_list])
        # used_time = time.time() - start_time
        # print('%.2fs'%used_time)
        return boot_stat
    ##############################################################################

    ##############################################################################
    def remove_traj_from_frame_list(self, rm_traj_list, frame_list, traj_axis):
        if len(frame_list) == 0:
            result_list = frame_list
        else:
            if traj_axis == 0: # traj ids in 1st row
                sel_idx = []
                for idx, i in enumerate(frame_list[0,:]):
                    if i not in rm_traj_list:
                        sel_idx.append(idx)
                result_list = frame_list[:,sel_idx]
            elif traj_axis == 1: # traj ids in 1st column
                sel_idx = []
                for idx, i in enumerate(frame_list[:,0]):
                    if i not in rm_traj_list:
                        sel_idx.append(idx)
                result_list = frame_list[sel_idx,:]
        return result_list
    ##############################################################################

    ##############################################################################
    def statistic_fun(self, data_1, data_2, ref):
        a = (np.mean(data_1) - np.mean(data_2)) - ref
        b = (np.std(data_1)**2/len(data_1) + np.std(data_2)**2/len(data_2))**0.5
        if a == 0 and b == 0:
            stat = 1.0
        else:
            stat = a/b
        return stat
    ##############################################################################

    ##############################################################################
    def bootstrap_test(self, data_1, data_2, statistic_fun, n_time, side='!='):
        if side not in ['!=', '>', '<']:
            self.logger.info('side parameter is wrong for function "bootstrap_test". It must be "!=", ">", or "<".')
            sys.exit()
        idx_list_1 = np.arange(len(data_1))
        idx_list_2 = np.arange(len(data_2))
        boot_stat = []
        ref_0 = 0
        ref_1 = np.mean(data_1) - np.mean(data_2)
        for i in range(n_time):
            sample_idx_list_1 = np.random.choice(idx_list_1, len(idx_list_1))
            sample_idx_list_2 = np.random.choice(idx_list_2, len(idx_list_2))
            bs = statistic_fun(data_1[sample_idx_list_1], data_2[sample_idx_list_2], ref_1)
            boot_stat.append(bs)
        boot_stat = np.array(boot_stat)
        
        boot_stat_0 = statistic_fun(data_1, data_2, ref_0)
        
        if side == '!=':
            p = (np.min([len(np.where(boot_stat >= boot_stat_0)[0]), len(np.where(boot_stat <= boot_stat_0)[0])])+1) / (n_time+1) * 2
        elif side == '>':
            p = (len(np.where(boot_stat >= boot_stat_0)[0])+1) / (n_time+1)
        elif side == '<':
            p = (len(np.where(boot_stat <= boot_stat_0)[0])+1) / (n_time+1)
        
        if p > 1:
            p = 1.0
        
        return (p, boot_stat)
    ##############################################################################


    ##############################################################################
    def LiP_XL_MS_ConsistencyTest(self,):
        self.logger.info(f'Comparing simulation to experimental data...')
        xlsx_outfile = os.path.join(self.outdir, f'LiPMS_XLMS_consist_pvalues_metastates_v11_down_sample_lag{self.lag_frame}.xlsx')
        npz_outfile = os.path.join(self.outdir, 'LiPMS_XLMS_consist_data_v9.npz')
        self.logger.debug(f'xlsx_outfile: {xlsx_outfile}')
        self.logger.debug(f'npz_outfile: {npz_outfile}')

        if os.path.exists(npz_outfile) and os.path.exists(xlsx_outfile):
            self.logger.info(f'npz_outfile EXISTS: Loading...')
            self.logger.info(f'xlsx_outfile EXISTS: Loading...')
            M_data = np.load(npz_outfile, allow_pickle=True)
            XLSX_df = pd.read_excel(xlsx_outfile)
            #print(f'XLSX_df:\n{XLSX_df}')

        else:

            #################################################################
            # Load MSM data
            MSM_data = pd.read_csv(self.msm_data_file)
            self.logger.info(f'MSM_data\n{MSM_data}')
            meta_states = MSM_data['metastablestate'].unique()
            meta_states = np.array(meta_states, dtype=int)
            self.logger.debug(f'meta_states: {meta_states}')
            num_meta_states = len(meta_states)
            self.logger.debug(f'num_meta_states: {num_meta_states}')
  

            meta_dtrajs_last = []
            traj_idx_to_trajnum = {} # mapping traj_idx to traj number
            for traj_idx, (traj, traj_df) in enumerate(MSM_data.groupby('traj')):
                traj_len = len(traj_df)
                self.logger.debug(f'traj: {traj}, traj_len: {traj_len}\n{traj_df.head()}')

                last = traj_df.iloc[-self.n_analysis_frames:,:]
                last = last.reset_index(drop=True)
                last = last['metastablestate'].values
                #print(f'last: {last}')
                meta_dtrajs_last.append(last)

                self.logger.debug(f'traj_idx: {traj_idx}, traj: {traj}, traj_len: {traj_len}, n_analysis_frames: {self.n_analysis_frames}, last: {last} {len(last)}')
                traj_idx_to_trajnum[traj_idx] = traj

            meta_dtrajs_last = np.array(meta_dtrajs_last)
            self.logger.info(f'meta_dtrajs_last.shape: {meta_dtrajs_last.shape}')
            self.logger.info(f'meta_dtrajs_last\n{meta_dtrajs_last} {meta_dtrajs_last.shape}')
            self.logger.debug(np.unique(meta_dtrajs_last))

            # Keep MSM trajectory indexing aligned with downstream OP arrays.
            rm_traj_set = set(int(t) for t in self.rm_traj_list)
            keep_traj_idx = [
                idx for idx in range(meta_dtrajs_last.shape[0])
                if int(traj_idx_to_trajnum[idx]) not in rm_traj_set
            ]
            meta_dtrajs_last = meta_dtrajs_last[keep_traj_idx, :]
            traj_idx_to_trajnum = {
                new_idx: int(traj_idx_to_trajnum[old_idx])
                for new_idx, old_idx in enumerate(keep_traj_idx)
            }
            # Row indices into the SASA/Jwalk OP arrays. CollectOP stores row
            # (traj_number - 1) for every trajectory 1..n_traj, so OP rows must be
            # selected by trajectory NUMBER rather than by MSM position: the MSM
            # mapping may already exclude mirror trajectories, which would otherwise
            # shift the positional alignment between MSM states and OP values.
            op_row_idx = np.array(
                [traj_idx_to_trajnum[i] - 1 for i in range(meta_dtrajs_last.shape[0])],
                dtype=int,
            )
            self.logger.info(
                f'meta_dtrajs_last.shape after mirror-image removal: {meta_dtrajs_last.shape}'
            )
            self.logger.info(f'op_row_idx n={len(op_row_idx)} first5={op_row_idx[:5]} last5={op_row_idx[-5:]}')
            #################################################################


            #################################################################
            # Create frame_list for each state
            # The result is a list where each element is a 2D array with the first column being the trajectory index and the second column being the frame index
            sel_frame_idx = np.arange(0, self.n_analysis_frames, self.lag_frame)
            self.logger.info(f'sel_frame_idx:\n{sel_frame_idx}')

            frame_list = []
            empty_states = []
            for state_idx in self.state_idx_list:
                self.logger.info(f'\nGetting frames for state {state_idx}...')
                frame_list_0 = np.array(np.where(meta_dtrajs_last[:, sel_frame_idx] == state_idx)).T 
                frame_list_0[:,1] = sel_frame_idx[frame_list_0[:,1]]
                utrajs = np.unique(frame_list_0[:,0])
                self.logger.debug(frame_list_0)
                self.logger.debug(f'utrajs: {utrajs}')

                #frame_list_0 = self.remove_traj_from_frame_list(self.rm_traj_list, frame_list_0, 1)
                if len(frame_list_0) == 0:
                    self.logger.info(f'No frames for state {state_idx} in the last {self.n_analysis_frames} frames. Exitting. The state maybe made entirely of mirror traj and in the self.rm_traj_list!')
                    empty_states.append(state_idx)
                    continue

                frame_list.append(frame_list_0)
            
            ## adjust the frame_list and the state_idx_list depending on what states are populated
            self.logger.debug(f'empty_states: {empty_states}')
            self.state_idx_list = [idx for idx in self.state_idx_list if idx not in empty_states]
            self.logger.info(f'Updated state_idx_list: {self.state_idx_list}')
            for idx, arr in enumerate(frame_list):
                state = self.state_idx_list[idx]
                self.logger.info(f'state_list_idx={idx}, state={state}, shape={arr.shape}')

            native_sel = np.array(np.where(meta_dtrajs_last[:, sel_frame_idx] == self.native_state_idx)).T
            native_sel[:,1] = sel_frame_idx[native_sel[:,1]]
            self.logger.info(f'native_sel:\n{native_sel} {native_sel.shape}')

            #native_sel = self.remove_traj_from_frame_list(self.rm_traj_list, native_sel, 1)
            #print(f'native_sel:\n{native_sel} {native_sel.shape}')
            #################################################################


            #################################################################
            # Load SASA data
            sasa_traj_list = np.load(self.sasa_data_file, allow_pickle=True)[:,-self.n_analysis_frames:,:]
            self.logger.info(f'sasa_traj_list.shape: {sasa_traj_list.shape}')

            # Align OP rows to the MSM trajectory order by trajectory number.
            if len(op_row_idx) and op_row_idx.max() >= sasa_traj_list.shape[0]:
                raise IndexError(
                    f'SASA array has {sasa_traj_list.shape[0]} trajectories but the MSM references '
                    f'trajectory number {int(op_row_idx.max()) + 1}. Ensure --n_traj covers all MSM trajectories.'
                )
            sasa_traj_list = sasa_traj_list[op_row_idx, :, :]
            self.logger.info(f'sasa_traj_list.shape after trajectory alignment: {sasa_traj_list.shape}')
            #################################################################


            #################################################################
            # XL-MS scores are computed either from a pre-built Jwalk.npy
            # (legacy mode) or directly from per-trajectory XP files
            # (memory-friendly streaming mode).
            dist_traj_list = None
            #################################################################


            #################################################################
            # Select the frames without nan residual SASA 
            # NAN SASA indicates a bad backmapped structure and I want to skip them
            nan_frame_sel = np.where(np.isnan(sasa_traj_list))
            nan_frame_sel = np.array([nan_frame_sel[0], nan_frame_sel[1]], dtype=int).T
            self.logger.debug(nan_frame_sel)
            frame_list = [np.array([sel for sel in frame if not sel.tolist() in nan_frame_sel.tolist()]) for frame in frame_list]
            self.logger.info(f'frame_list:')
            for idx, arr in enumerate(frame_list):
                state = self.state_idx_list[idx]
                self.logger.info(f'state_list_idx={idx}, state={state}, shape={arr.shape}')
            native_sel = np.array([sel for sel in native_sel if not sel.tolist() in nan_frame_sel.tolist()])
            self.logger.info(f'native_sel:\n{native_sel} {native_sel.shape}')
            #################################################################


            #################################################################
            # Load LiPMS experimental data
            LiPMS_sig_data = self.load_LiPMS_data(self.LiPMS_exp_file)
            self.logger.info(f'Loaded LiPMS experimental data: {len(LiPMS_sig_data)} peptides')

            # load XLMS experimental data
            XLMS_sig_data = self.load_XLMS_data(self.XLMS_exp_file)
            self.logger.info(f'Loaded XLMS experimental data: {len(XLMS_sig_data)} pairs')
            #################################################################
            

            #################################################################
            # Calculate or load metric matrix
            if self.if_calc_M == 1:
                self.logger.debug('Calculating metric matrix...')
                # calculate metric matrix
                M_LiPMS = np.zeros((*meta_dtrajs_last.shape, len(LiPMS_sig_data)))
                for idx, peptide in enumerate(LiPMS_sig_data.keys()):
                    sel = list(LiPMS_sig_data[peptide]['peptide_range'])
                    SA = np.sum(sasa_traj_list[:,:,sel], axis=-1)
                    M_LiPMS[:,:,idx] = SA
                
                xlms_targets = []
                for idx, key in enumerate(XLMS_sig_data.keys()):
                    pair_AA = [k[0] for k in key.split('-')]
                    key_0 = '-'.join([k[1:]+'|A' for k in key.split('-')])
                    xlms_targets.append((idx, key_0, pair_AA))

                M_XLMS = np.zeros((*meta_dtrajs_last.shape, len(XLMS_sig_data)))
                if self.dist_data_file is not None and os.path.exists(self.dist_data_file):
                    dist_traj_list = np.load(self.dist_data_file, allow_pickle=True)[:,-self.n_analysis_frames:]
                    self.logger.info(f'dist_traj_list.shape: {dist_traj_list.shape}')

                    # Align OP rows to the MSM trajectory order by trajectory number.
                    if len(op_row_idx) and op_row_idx.max() >= dist_traj_list.shape[0]:
                        raise IndexError(
                            f'Jwalk array has {dist_traj_list.shape[0]} trajectories but the MSM references '
                            f'trajectory number {int(op_row_idx.max()) + 1}. Ensure --n_traj covers all MSM trajectories.'
                        )
                    dist_traj_list = dist_traj_list[op_row_idx, :]
                    self.logger.info(f'dist_traj_list.shape after trajectory alignment: {dist_traj_list.shape}')

                    n_traj = min(meta_dtrajs_last.shape[0], dist_traj_list.shape[0])
                    n_frame = min(meta_dtrajs_last.shape[1], dist_traj_list.shape[1])
                    for i in range(n_traj):
                        for j in range(n_frame):
                            frame_data = dist_traj_list[i, j]
                            if frame_data is None:
                                continue
                            for idx, key_0, pair_AA in xlms_targets:
                                if key_0 not in frame_data:
                                    continue
                                JWalk_dist = frame_data[key_0].get('Jwalk', -1)
                                M_XLMS[i, j, idx] = self.score_XL(pair_AA, JWalk_dist)
                elif self.xp_dir is not None and os.path.isdir(self.xp_dir):
                    self.logger.info(f'Computing XL-MS scores in streaming mode from XP files in: {self.xp_dir}')
                    target_lookup = {key_0: (idx, pair_AA) for idx, key_0, pair_AA in xlms_targets}

                    for traj_idx in range(meta_dtrajs_last.shape[0]):
                        traj_num = int(traj_idx_to_trajnum[traj_idx])
                        if traj_num in rm_traj_set:
                            continue

                        fpath = os.path.join(self.xp_dir, f'{self.ID}_Traj{traj_num}.XP')
                        if not os.path.exists(fpath):
                            self.logger.warning(f'Missing XP file for streaming XL-MS scoring: {fpath}')
                            continue

                        df = pd.read_csv(
                            fpath,
                            sep='\t',
                            usecols=['Frame', 'Atom1', 'Atom2', 'SASD'],
                            dtype={'Frame': np.int32, 'SASD': np.float32, 'Atom1': 'string', 'Atom2': 'string'},
                        )
                        if df.empty:
                            continue

                        frame_values = np.sort(df['Frame'].unique())
                        frame_to_idx = {int(f): idx for idx, f in enumerate(frame_values)}

                        atom1_parts = df['Atom1'].str.split('-', expand=True)
                        atom2_parts = df['Atom2'].str.split('-', expand=True)
                        if atom1_parts.shape[1] < 3 or atom2_parts.shape[1] < 3:
                            self.logger.warning(f'Malformed Atom fields in {fpath}; skipping trajectory {traj_num}')
                            continue

                        frame_key_df = pd.DataFrame({
                            'Frame': df['Frame'].values,
                            'SASD': df['SASD'].values,
                            'key': (atom1_parts[1].astype(str) + '|'+ atom1_parts[2].astype(str) + '-' +
                                    atom2_parts[1].astype(str) + '|' + atom2_parts[2].astype(str)).values,
                        })

                        frame_key_df = frame_key_df[frame_key_df['key'].isin(target_lookup.keys())]
                        frame_key_df = frame_key_df.drop_duplicates(subset=['Frame', 'key'], keep='first')
                        if frame_key_df.empty:
                            continue

                        for key_0, key_df in frame_key_df.groupby('key', sort=False):
                            idx, pair_AA = target_lookup[key_0]
                            for frame_num, sasd in zip(key_df['Frame'].values, key_df['SASD'].values):
                                frame_idx = frame_to_idx.get(int(frame_num), None)
                                if frame_idx is None or frame_idx >= meta_dtrajs_last.shape[1]:
                                    continue
                                M_XLMS[traj_idx, frame_idx, idx] = self.score_XL(pair_AA, float(sasd))

                        if (traj_idx + 1) % 50 == 0:
                            self.logger.info(f'Streamed XP scoring progress: {traj_idx + 1}/{meta_dtrajs_last.shape[0]} trajectories')
                else:
                    raise ValueError(
                        'XL-MS scoring requires either dist_data_file (Jwalk.npy) or xp_dir (per-trajectory XP files).'
                    )

                # Save data
                np.savez(npz_outfile, 
                        M_LiPMS = M_LiPMS,
                        M_XLMS = M_XLMS)
                self.logger.info(f'Saved metric matrices to {npz_outfile}')

            else:
                self.logger.debug('Loading metric matrix...')
                M_data = np.load(npz_outfile, allow_pickle=True)
                M_LiPMS = M_data['M_LiPMS'][:,-self.n_analysis_frames:,:]
                M_XLMS = M_data['M_XLMS'][:,-self.n_analysis_frames:,:]
            #################################################################


            #################################################################
            # Consistency tests
            M_str_list = ['LiPMS', 'XLMS']
            exp_data_list = [LiPMS_sig_data, XLMS_sig_data]
            df_list = []
            df_all_state_list = []
            p_list = []
            p_all_state_list = []
            for idx, M in enumerate([M_LiPMS, M_XLMS]):
                self.logger.info(M_str_list[idx]+':')
                exp_data = exp_data_list[idx]
                
                df_list_0 = []
                for idx_0 in range(len(exp_data)):
                    self.logger.info('%s:'%(list(exp_data.keys())[idx_0]))
                    M_0 = M[:,:,idx_0]
                    index_str = list(exp_data.keys())[idx_0]
                    header = ['Near-native state', 'Sample size', '<M>']
                    header += ['p (!=)', 'Adjusted p (!=)']
                    if exp_data[list(exp_data.keys())[idx_0]]['qual_change'] > 0:
                        header += ['p (>)', 'Adjusted p (>)']
                        test_side = '>'
                    else:
                        header += ['p (<)', 'Adjusted p (<)']
                        test_side = '<'

                    M_0_native = M_0[native_sel[:,0], native_sel[:,1]]
                    # bootstrapping to get 95%CI
                    boot_stat_native = self.bootstrap(np.mean, M_0_native, self.n_boot)
                    lb_native = np.percentile(boot_stat_native, 2.5)
                    ub_native = np.percentile(boot_stat_native, 97.5)
                    
                    df_data = []
                    df_all_state_data = []
                    # for near-native states
                    self.logger.info(f'len(self.state_idx_list): {self.state_idx_list} {len(self.state_idx_list)}')
                    for idx_1, state_id in enumerate(self.state_idx_list):
                        self.logger.info(f'state_list_idx={idx_1}, state={state_id}')
                        near_native_sel = np.array(frame_list[idx_1], dtype=int)
                        self.logger.debug(near_native_sel)
                        M_0_near_native = M_0[near_native_sel[:,0], near_native_sel[:,1]]
                        # bootstrapping to get 95%CI
                        boot_stat_near_native = self.bootstrap(np.mean, M_0_near_native, self.n_boot)
                        lb_near_native = np.percentile(boot_stat_near_native, 2.5)
                        ub_near_native = np.percentile(boot_stat_near_native, 97.5)
                        self.logger.info('Near-native state %d vs. Native state %d:'%(state_id+1, num_meta_states))
                        self.logger.info('    Sample size: %d vs. %d'%(len(M_0_near_native), len(M_0_native)))
                        self.logger.info('    <M>: %.4f [%.4f, %.4f] vs. %.4f [%.4f, %.4f]'%(np.mean(M_0_near_native), lb_near_native, ub_near_native, np.mean(M_0_native), lb_native, ub_native))
                        p_value_list_0 = []
                        for ts in ['!=', test_side]:
                            p = self.permutation_test(self.perm_fun, M_0_near_native, M_0_native, self.num_perm, side=ts)
                            # p, _ = bootstrap_test(M_0_near_native, M_0_native, statistic_fun, n_boot, side=ts)
                            self.logger.info('    p-value ("%s") = %.4f'%(ts, p))
                            p_value_list_0.append(p)
                        p_list.append(p_value_list_0)
                        df_data.append([state_id+1, '%d vs. %d'%(len(M_0_near_native), len(M_0_native)), '%.4f [%.4f, %.4f] vs. %.4f [%.4f, %.4f]'%(np.mean(M_0_near_native), lb_near_native, ub_near_native, np.mean(M_0_native), lb_native, ub_native), p_value_list_0[0], 0, p_value_list_0[1], 0])
                        
                    df = pd.DataFrame(df_data, columns=header, index=[index_str]*len(df_data))
                    df_list_0.append(df)

                    # for all states
                    frame_list_1 = np.array(np.where(meta_dtrajs_last < num_meta_states)).T
                    #near_native_sel = self.remove_traj_from_frame_list(self.rm_traj_list, frame_list_1, 1)
                    near_native_sel = np.array([sel for sel in near_native_sel if not sel in nan_frame_sel]) # remove nan SASA
                    
                    # Select only frames seperated by #lag_frame for each state in the trajectory
                    near_native_sel_idx = []
                    for i in np.unique(near_native_sel[:,0]):
                        idx = np.where(near_native_sel[:,0] == i)[0]
                        near_native_sel_idx.append(idx[0])
                        for iidx in idx[1:]:
                            if near_native_sel[iidx,1] - near_native_sel[near_native_sel_idx[-1],1] >= self.lag_frame:
                                near_native_sel_idx.append(iidx)
                    near_native_sel = near_native_sel[near_native_sel_idx,:]

                    M_0_near_native = M_0[near_native_sel[:,0], near_native_sel[:,1]]
                    # bootstrapping to get 95%CI
                    boot_stat_near_native = self.bootstrap(np.mean, M_0_near_native, self.n_boot)
                    lb_near_native = np.percentile(boot_stat_near_native, 2.5)
                    ub_near_native = np.percentile(boot_stat_near_native, 97.5)
                    self.logger.info('All states vs. Native state %d:'%(num_meta_states))
                    self.logger.info('    Sample size: %d vs. %d'%(len(M_0_near_native), len(M_0_native)))
                    self.logger.info('    <M>: %.4f [%.4f, %.4f] vs. %.4f [%.4f, %.4f]'%(np.mean(M_0_near_native), lb_near_native, ub_near_native, np.mean(M_0_native), lb_native, ub_native))
                    p_value_list_0 = []
                    for ts in ['!=', test_side]:
                        p = self.permutation_test(self.perm_fun, M_0_near_native, M_0_native, self.num_perm, side=ts)
                        # p, _ = bootstrap_test(M_0_near_native, M_0_native, statistic_fun, n_boot, side=ts)
                        self.logger.info('    p-value ("%s") = %.4f'%(ts, p))
                        p_value_list_0.append(p)
                    p_all_state_list.append(p_value_list_0)
                    df_all_state_data.append(['All states', '%d vs. %d'%(len(M_0_near_native), len(M_0_native)), '%.4f [%.4f, %.4f] vs. %.4f [%.4f, %.4f]'%(np.mean(M_0_near_native), lb_near_native, ub_near_native, np.mean(M_0_native), lb_native, ub_native), p_value_list_0[0], 0, p_value_list_0[1], 0])
                    df = pd.DataFrame(df_all_state_data, columns=header, index=[index_str]*len(df_all_state_data))
                    df_all_state_list.append(df)
                df_list.append(df_list_0)
            #################################################################


            #################################################################
            # Correct p-value
            p_list = np.array(p_list)
            p_all_state_list = np.array(p_all_state_list)
            adjusted_p_list = []
            for pi in range(p_list.shape[1]):
                pl = p_list[:,pi]
                results = multipletests(pl, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
                adjusted_p_list.append(results[1])
            adjusted_p_list = np.array(adjusted_p_list).T
            adjusted_p_all_state_list = []
            for pi in range(p_all_state_list.shape[1]):
                pl = p_all_state_list[:,pi]
                results = multipletests(pl, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
                adjusted_p_all_state_list.append(results[1])
            adjusted_p_all_state_list = np.array(adjusted_p_all_state_list).T
            #################################################################

            #################################################################
            # Update dataframes
            pi = 0
            for idx in range(len(df_list)):
                for idx_0 in range(len(df_list[idx])):
                    header = list(df_list[idx][idx_0].columns)
                    for i in range(adjusted_p_list.shape[1]):
                        ii = -adjusted_p_list.shape[1]*2 + i*2 + 1
                        df_list[idx][idx_0][header[ii]] = adjusted_p_list[pi:pi+len(self.state_idx_list),i]
                    pi += len(df_list[idx][idx_0])

            pi = 0
            for df in df_all_state_list:
                header = list(df.columns)
                for i in range(adjusted_p_all_state_list.shape[1]):
                    ii = -adjusted_p_all_state_list.shape[1]*2 + i*2 + 1
                    df[header[ii]] = adjusted_p_all_state_list[pi:pi+len(df),i]
                pi += len(df)

            df_list.append(df_all_state_list)
            M_str_list.append('All states')
            #################################################################

            #################################################################
            # Creating Excel Writer Object from Pandas  
            with pd.ExcelWriter(xlsx_outfile, engine="openpyxl") as writer:
                workbook=writer.book
                for idx_0, df_list_0 in enumerate(df_list):
                    worksheet=workbook.create_sheet(M_str_list[idx_0])
                    writer.sheets[M_str_list[idx_0]] = worksheet
                    start_row = 0
                    for idx_1, df in enumerate(df_list_0):
                        df.to_excel(writer, sheet_name=M_str_list[idx_0], index=True, float_format="%.4f", startrow=start_row , startcol=0)
                        start_row += len(df) + 2
            self.logger.debug(f'SAVED: {xlsx_outfile}')
            #################################################################

        self.logger.info(f'LiP_XL_MS_ConsistencyTest DONE!')
        return npz_outfile, xlsx_outfile
    ##############################################################################

    #######################################################################################
    def load_OP(self, start:int=0, end:int=99999999999):
        """
        Loads the GQ values of each trajectory into a 2D array and then appends it to a list
        The list should have Nt = number of trajectories and each array should be n x 2 where n is the number of frames
        """
        self.logger.info(f'Loading G and Q order parameters...')
        Qfiles = glob.glob(os.path.join(self.OPpath, 'Q/*.Q'))
        QTrajs = [int(pathlib.Path(Qf).stem.split('Traj')[-1]) for Qf in Qfiles]

        Gfiles = glob.glob(os.path.join(self.OPpath, 'G/*.G'))
        GTrajs = [int(pathlib.Path(Gf).stem.split('Traj')[-1]) for Gf in Gfiles]

        shared_Trajs = set(QTrajs).intersection(GTrajs)
        shared_Trajs = sorted(shared_Trajs)
        self.logger.debug(f'shared_Trajs: {shared_Trajs}')
        #print(f'Shared Traj between Q and G: {shared_Trajs} {len(shared_Trajs)}')
        self.logger.info(f'Number of Q files found: {len(Qfiles)} | Number of G files found: {len(Gfiles)}')

        assert len(Qfiles) == len(Gfiles), f"The # of Q and G files {len(Qfiles)} and {len(Gfiles)} is not equal"

        ## remove trajectories that are in the rm_traj_list
        if len(self.rm_traj_list) > 0:
            self.logger.info(f'Removing trajectories: {self.rm_traj_list}')
            shared_Trajs = [traj for traj in shared_Trajs if traj not in self.rm_traj_list]
            self.logger.info(f'Number of shared Traj after removing: {len(shared_Trajs)}')

        # loop through the Qfiles and find matching Gfile
        # then load the Q and G time series into a 2D array
        Q_data = []
        G_data = []
        for traj in shared_Trajs:
            #print(f'Traj: {traj}')

            # get the cooresponding G and Q file
            Qf = [f for f in Qfiles if f.endswith(f'Traj{traj}.Q')]
            Gf = [f for f in Gfiles if f.endswith(f'Traj{traj}.G')]
            #print(f'Qf: {Qf}')
            #print(f'Gf: {Gf}')

            ## Quality check to assert that only a single G and Q file were found
            assert len(Qf) == 1, f"the number of Q files {len(Qf)} should equal 1 for Traj {traj}"
            assert len(Gf) == 1, f"the number of G files {len(Gf)} should equal 1 for Traj {traj}"

            # load the G Q data and extract only the time series column
            #Qdata = pd.read_csv(Qf)['Q'].values[self.start:self.end:self.stride]
            Qdata = pd.read_csv(Qf[0], sep=',')
            Gdata = pd.read_csv(Gf[0])
            if start < 0: # start was specified as negative and takens as slicing the end of the arry
                Qdata = Qdata['total'].values[start:]
                Gdata = Gdata['G'].values[start:]
            else:
                Qdata = Qdata[(Qdata['Frame'] >= start) & (Qdata['Frame'] <= end)]
                Qdata = Qdata['total'].values.astype(float)

                Gdata = Gdata[(Gdata['Frame'] >= start) & (Gdata['Frame'] <= end)]
                Gdata = Gdata['G'].values.astype(float)
            # print(f'Shape of OP: Q {Qdata.shape} G {Gdata.shape}')

            ## Quality check that the G and Q data has the same number of frames
            if Qdata.shape != Gdata.shape:
                self.logger.warning(f"WARNING: The number of frames in Q {Qdata.shape} should equal the number of frames in G {Gdata.shape} in Traj {traj}")
                continue
            
            ## Check and ensure that Qdata or Gdata has no nan values
            if np.isnan(Qdata).any():
                raise ValueError(f'There is a NaN value in this Qdata')
             
            if np.isnan(Gdata).any():
                raise ValueError(f'There is a NaN value in this Gdata')
         
            Q_data.append(Qdata)
            G_data.append(Gdata)
        
        Q_data = np.asarray(Q_data)
        G_data = np.asarray(G_data)
        self.logger.debug(f'Q_data: {Q_data.shape}')
        self.logger.debug(f'G_data: {G_data.shape}')
        return Q_data, G_data
    ##############################################################################

    ##############################################################################
    def select_rep_structs(self, consist_data_file:str, consist_result_file:str, total_traj_num_frames:int, n_analysis_frames:int, restart: bool = False):
        """
        After performing the consistency test select representative structures with high consistency.
        When restart=True, an existing viz_rep_struct directory is preserved and any group directory
        that already contains a .done sentinel file is skipped.
        """
        self.logger.info(f'Selecting representative structure')

        ############### Functions ###############
        ##############################################################################
        def parse_consistant_results(consist_result_file, sheet_name, test_type='two-tailed', significant_level=0.05):
            data = {}
            consist_data = pd.read_excel(consist_result_file, sheet_name=sheet_name)
            num_row = consist_data.index.size
            all_key_list = []
            if test_type == 'two-tailed':
                t_idx = -3
            elif test_type == 'one-tailed':
                t_idx = -1
            else:
                self.logger.error('Error: Wrong test_type = %s; can be either "two-tailed" or "one-tailed"'%test_type)
                sys.exit()
            sign = consist_data.columns[t_idx].split()[-1][1:-1]
            for i in range(num_row):
                row = consist_data.loc[i]
                if type(row[-1]) == float and np.isnan(row[-1]):
                    continue
                elif type(row[0]) == float and np.isnan(row[0]):
                    sign = row[t_idx].split()[-1][1:-1]
                elif row[1] == 'All states':
                    continue
                else:
                    if row[0] not in all_key_list:
                        all_key_list.append(row[0])
                    if row[t_idx] < significant_level:
                        if row[0] not in data.keys():
                            mean = float(row[3].strip().split('vs.')[-1].split()[0])
                            words = row[3].strip().split('vs.')[-1].split()[1:]
                            lb = float(words[0][1:-1])
                            ub = float(words[1][:-1])
                            data[row[0]] = [sign, [row[1]-1], mean, [lb, ub]]
                        else:
                            data[row[0]][1].append(row[1]-1)
            
            for key in data.keys():
                data[key].append(all_key_list.index(key))
            
            return data
        ##############################################################################

        ##############################################################################
        def calc_rel_change(traj_idx, frame_idx):
            signal_list = dtrajs_MS[traj_idx, frame_idx]
            rel_change_list = []
            for signal in signal_list:
                if '-' in signal[1:-1]:
                    data = XLMS_consist_data
                    M = M_XLMS
                else:
                    data = LIPMS_consist_data
                    M = M_LiPMS
                mean = data[signal][2]
                v = M[traj_idx, frame_idx, data[signal][4]]
                if mean == 0:
                    rel_change_list.append(np.abs(v/1e-5))
                else:
                    rel_change_list.append(np.abs(v/mean-1))
            return rel_change_list
        ##############################################################################


        ################### MAIN ###############
        ##############################################################################
        self.logger.info(f'Loading consistency test data from {consist_result_file}')
        self.logger.info(f'Loading consistency test data from {consist_data_file}')

        if_backmap = 0
        pulchra_only = True
        significant_level = 0.05
        ##############################################################################

        ##############################################################################
        # Load MSM data
        MSM_data = pd.read_csv(self.msm_data_file)
        self.logger.info(f'MSM_data\n{MSM_data}')
        meta_states = MSM_data['metastablestate'].unique()
        meta_states = np.array(meta_states, dtype=int)
        self.logger.debug(f'meta_states: {meta_states}')
        num_meta_states = len(meta_states)
        self.logger.debug(f'num_meta_states: {num_meta_states}')

        rm_traj_set = set(int(t) for t in self.rm_traj_list)
        meta_dtrajs_last = []
        micro_dtrajs_last = []
        MSM_traj_idx_to_trajnum = {} # mapping traj_idx to traj number (after rm_traj_list filtering)
        for traj, traj_df in MSM_data.groupby('traj'):
            traj = int(traj)
            if traj in rm_traj_set:
                continue
            traj_len = len(traj_df)
            #print(f'traj: {traj}, traj_len: {traj_len}\n{traj_df.head()}')

            last = traj_df.iloc[-n_analysis_frames:,:]
            last = last.reset_index(drop=True)
            meta_last = last['metastablestate'].values
            micro_last = last['microstate'].values
            #print(f'last: {last}')
            meta_dtrajs_last.append(meta_last)
            micro_dtrajs_last.append(micro_last)

            MSM_traj_idx_to_trajnum[len(MSM_traj_idx_to_trajnum)] = traj

        meta_dtrajs_last = np.array(meta_dtrajs_last)
        micro_dtrajs_last = np.array(micro_dtrajs_last)
        self.logger.info(f'meta_dtrajs_last\n{meta_dtrajs_last} {meta_dtrajs_last.shape}')
        self.logger.debug(np.unique(meta_dtrajs_last))
        self.logger.info(f'micro_dtrajs_last\n{micro_dtrajs_last} {micro_dtrajs_last.shape}')
        self.logger.debug(np.unique(micro_dtrajs_last))
        self.logger.debug(f'MSM_traj_idx_to_trajnum: {MSM_traj_idx_to_trajnum}')

        ## load the meta_dist data
        meta_dist = np.load(self.meta_dist_file, allow_pickle=True)
        self.logger.info(f'meta_dist:\n{meta_dist} {meta_dist.shape}')
        ##############################################################################

 
        ##############################################################################
        # Load Consistency Metrics
        consist_data = np.load(consist_data_file, allow_pickle=True)
        M_LiPMS = consist_data['M_LiPMS'][:,-n_analysis_frames:,:]
        M_XLMS = consist_data['M_XLMS'][:,-n_analysis_frames:,:]
        self.logger.info(f'Loaded consistency metrics from {consist_data_file}')
        self.logger.debug(f'M_LiPMS: {M_LiPMS.shape}')
        self.logger.debug(f'M_XLMS: {M_XLMS.shape}')
        ##############################################################################

        ##############################################################################
        # Load consistency test results
        LIPMS_consist_data = parse_consistant_results(consist_result_file, sheet_name='LiPMS', test_type='two-tailed', significant_level=significant_level)
        XLMS_consist_data = parse_consistant_results(consist_result_file, sheet_name='XLMS', test_type='two-tailed', significant_level=significant_level)
        self.logger.info(f'Loaded consistency test results')
        ##############################################################################
     

        ##############################################################################
        # Load cluster data
        # Beware the order of traj in the cluster_data may not be the same as the order in the MSM_data
        # we will correct that here so it is consistent moving forward with the analysis
        cluster_data = np.load(self.cluster_data_file, allow_pickle=True)
        idx2trajfile = cluster_data['idx2trajfile'].tolist()
        idx2traj = np.asarray([int(f.strip().split('/')[-1].split('_')[0]) for f in idx2trajfile])
        # print(f'idx2traj: {idx2traj} {len(idx2traj)}')

        # get the ordering array and reorder all the cluster_data arrays
        # also remove those trajectories that are in the rm_traj_list
        order = np.argsort(idx2traj)
        # print(f'order: {order} {len(order)}')
        idx2traj = idx2traj[order]
        # print(f'idx2traj: {idx2traj} {len(idx2traj)}')

        idx_2_keep = [idx for idx, traj in enumerate(idx2traj) if traj not in self.rm_traj_list]
        # print(f'idx_2_keep: {idx_2_keep} {len(idx_2_keep)}')

        # reorder idx2trajfile
        idx2traj = idx2traj[idx_2_keep]
        # print(f'idx2traj after removal of mirror images: {idx2traj} {len(idx2traj)}')

        # QC to ensure idx2traj matches MSM_traj_idx_to_trajnum
        assert np.array_equal(idx2traj, np.asarray(list(MSM_traj_idx_to_trajnum.values()))), f"Error: idx2traj does not match MSM_traj_idx_to_trajnum: {idx2traj} vs. {np.asarray(list(MSM_traj_idx_to_trajnum.values()))}"

        # load the dtraj data
        dtrajs = cluster_data['dtrajs']
        dtrajs = dtrajs[order]
        dtrajs = dtrajs[idx_2_keep]
        self.logger.info(f'dtrajs after removal of mirror images: {dtrajs.shape}')

        # load the rep_chg_ent_dtrajs
        rep_chg_ent_dtrajs = cluster_data['rep_chg_ent_dtrajs']
        # Map resid to residue idx using resid2residueidx_map
        self.logger.info(f'Mapping of rep_chg_ent_dtrajs resid to residu idx using resid2residueidx_map: {self.resid2residueidx_map}')
        for traj_idx, traj_data in enumerate(rep_chg_ent_dtrajs):
            for frame_idx, frame_data in enumerate(traj_data):
                for fingerprint_id, fingerprint in frame_data.items():
                    # print(f'\n{"#"*50}\nTraj idx: {traj_idx} | Frame idx {frame_idx} | Fingerprint ID: {fingerprint_id}')

                    for fingerprint_key, new_key in {'crossing_resid':'crossing_residx', 'ref_crossing_resid':'ref_crossing_residx'}.items():
                        residx_arr = fingerprint[fingerprint_key]
                        residx_arr = [[self.resid2residueidx_map[x] for x in sublist] for sublist in residx_arr]
                        fingerprint[new_key] = residx_arr
                        # print(f'      mapped {new_key}: {residx_arr}')

                    for fingerprint_key, new_key in {'native_contact':'native_contact_residx', 'ref_native_contact':'ref_native_contact_residx'}.items():
                        residx_arr = fingerprint[fingerprint_key]
                        residx_arr = [self.resid2residueidx_map[x] for x in residx_arr]
                        fingerprint[new_key] = residx_arr
                        # print(f'      mapped {new_key}: {residx_arr}')

                    # for k,v in fingerprint.items():
                    #     print(f'      {k}: {v}')
            
        rep_chg_ent_dtrajs = rep_chg_ent_dtrajs[order]
        rep_chg_ent_dtrajs = rep_chg_ent_dtrajs[idx_2_keep]
        self.logger.info(f'rep_chg_ent_dtrajs after removal of mirror images: {rep_chg_ent_dtrajs.shape}')

        sorted_chg_ent_structure_keyword_list = cluster_data['sorted_chg_ent_structure_keyword_list'].tolist()
        self.logger.debug(f'sorted_chg_ent_structure_keyword_list: {len(sorted_chg_ent_structure_keyword_list)} {sorted_chg_ent_structure_keyword_list[:10]}')

        dtrajs_cluster_idx = np.array([[sorted_chg_ent_structure_keyword_list.index(str(dd)) for dd in d] for d in dtrajs])
        self.logger.debug(f'dtrajs_cluster_idx: {dtrajs_cluster_idx.shape}')
        self.logger.info(f'Loaded cluster_data_file: {self.cluster_data_file}')
        ##############################################################################

        ##############################################################################
        # Load SASA data
        sasa_traj_list = np.load(self.sasa_data_file, allow_pickle=True)[:,-n_analysis_frames:,:]
        self.logger.info(f'Loaded SASA data: {self.sasa_data_file}')

        # remove trajectories that are in the rm_traj_list - 1
        sasa_traj_list = [v for i, v in enumerate(sasa_traj_list) if i not in np.asarray(self.rm_traj_list) - 1]
        sasa_traj_list = np.array(sasa_traj_list)
        self.logger.info(f'sasa_traj_list.shape after removal of mirror images: {sasa_traj_list.shape}')
        ##############################################################################


        ##############################################################################
        # Change the ub and lb in LIPMS_consist_data and XLMS_consist_data
        # (1) Get the native state index (defined by the MSM indexing which is sorted in order of the trajectory number)
        native_sel = np.where(np.isin(meta_dtrajs_last, self.native_state_idx))
        native_sel = np.array([native_sel[0], native_sel[1]], dtype=int).T

        # (2) remove any frames with NAN residual SASA
        nan_frame_sel = np.where(np.isnan(sasa_traj_list))
        nan_frame_sel = np.array([nan_frame_sel[0], nan_frame_sel[1]], dtype=int).T
        native_sel = np.array([sel for sel in native_sel if not sel.tolist() in nan_frame_sel.tolist()])
        self.logger.debug(f'native_sel: {native_sel} {native_sel.shape}')

        # (3) Get the native M_LiPMS and M_XLMS
        native_M_LiPMS = M_LiPMS[native_sel[:,0], native_sel[:,1], :]
        native_M_XLMS = M_XLMS[native_sel[:,0], native_sel[:,1], :]
        native_M_data_outfile = os.path.join(self.outdir, 'native_M_data.npz')
        np.savez(native_M_data_outfile,
                M_LiPMS = native_M_LiPMS,
                M_XLMS = native_M_XLMS,)
        for key in LIPMS_consist_data:
            LIPMS_consist_data[key][3][0] = np.percentile(native_M_LiPMS[:, LIPMS_consist_data[key][4]], 2.5)
            LIPMS_consist_data[key][3][1] = np.percentile(native_M_LiPMS[:, LIPMS_consist_data[key][4]], 97.5)
        for key in XLMS_consist_data:
            XLMS_consist_data[key][3][0] = np.percentile(native_M_XLMS[:, XLMS_consist_data[key][4]], 2.5)
            XLMS_consist_data[key][3][1] = np.percentile(native_M_XLMS[:, XLMS_consist_data[key][4]], 97.5)
        self.logger.debug(f'SAVED: {native_M_data_outfile}')
        ##############################################################################

        ##############################################################################
        # Go through LiPMS data
        LIPMS_struct_data = {}
        dtrajs_MS = np.empty(meta_dtrajs_last.shape, dtype=object)

        # Initialize dtrajs_MS with empty lists
        for i in range(len(dtrajs_MS)):
            for j in range(len(dtrajs_MS[i])):
                dtrajs_MS[i,j] = []

        # Go through each PK site in the LiPMS consistency data
        self.logger.info(f'\nProcessing LiPMS consistency data')
        for pk, data in LIPMS_consist_data.items():
            # print(pk, data)

            LIPMS_struct_data[pk] = {}
            sign = data[0]
            state_list = data[1]
            sasa_ub = data[3][1]
            sasa_lb = data[3][0]
            idx_pk = data[4]

            for idx_0, state_idx in enumerate(state_list):
                LIPMS_struct_data[pk][state_idx] = {}

                idx_list = np.array(np.where(meta_dtrajs_last == state_idx)).T
                idx_list = np.array([sel for sel in idx_list if not sel.tolist() in nan_frame_sel.tolist()]) # Skip frames with NAN residual SASA (bad backmapped structure)

                sasa_list = M_LiPMS[idx_list[:,0],idx_list[:,1], idx_pk]

                for cluster_idx, keyword in enumerate(sorted_chg_ent_structure_keyword_list):
                    idx_list_1 = np.where(dtrajs_cluster_idx[idx_list[:,0],idx_list[:,1]] == cluster_idx)[0]

                    if len(idx_list_1) == 0:
                        continue
                    if sign == '>':
                        idx_list_2 = np.where(sasa_list[idx_list_1] > sasa_ub)[0]
                        idx_rep = np.argmax(sasa_list[idx_list_1])
                    elif sign == '<':
                        idx_list_2 = np.where(sasa_list[idx_list_1] < sasa_lb)[0]
                        idx_rep = np.argmin(sasa_list[idx_list_1])
                    else:
                        idx_list_2 = np.where(np.any([sasa_list[idx_list_1] > sasa_ub, sasa_list[idx_list_1] < sasa_lb], axis=0))[0]
                        idx_rep = np.argmax(np.max([sasa_list[idx_list_1]-sasa_ub, sasa_lb-sasa_list[idx_list_1]], axis=0))
                    if len(idx_list_2) == 0:
                        continue

                    consist_idx_list = idx_list[idx_list_1[idx_list_2], :]
                    rep_idx = idx_list[idx_list_1[idx_rep], :]
                    LIPMS_struct_data[pk][state_idx][cluster_idx] = [rep_idx, consist_idx_list]
                    # print(f'PK: {pk}, state_idx: {state_idx}, cluster_idx: {cluster_idx}, rep_idx: {rep_idx}, consist_idx_list: {len(consist_idx_list)}')

                    for idx in consist_idx_list:
                        dtrajs_MS[idx[0],idx[1]].append(pk)
        ##############################################################################

        ##############################################################################
        # Go through XLMS data
        XLMS_struct_data = {}
        self.logger.info(f'\nProcessing XLMS consistency data')
        for pair, data in XLMS_consist_data.items():
            # print(pair, data)

            # Initialize the structure data for each pair
            XLMS_struct_data[pair] = {}
            sign = data[0]
            state_list = data[1]
            dist_ub = data[3][1]
            dist_lb = data[3][0]
            idx_pair = data[4]

            for state_idx in state_list:
                XLMS_struct_data[pair][state_idx] = {}

                idx_list = np.array(np.where(meta_dtrajs_last == state_idx)).T
                idx_list = np.array([sel for sel in idx_list if not sel.tolist() in nan_frame_sel.tolist()])

                dist_list = M_XLMS[idx_list[:,0],idx_list[:,1],idx_pair]
                for cluster_idx, keyword in enumerate(sorted_chg_ent_structure_keyword_list):
                    idx_list_1 = np.where(dtrajs_cluster_idx[idx_list[:,0],idx_list[:,1]] == cluster_idx)[0]
                    if len(idx_list_1) == 0:
                        continue
                    if sign == '>':
                        idx_list_2 = np.where(dist_list[idx_list_1] > dist_ub)[0]
                        idx_rep = np.argmax(dist_list[idx_list_1])
                    elif sign == '<':
                        idx_list_2 = np.where(dist_list[idx_list_1] < dist_lb)[0]
                        idx_rep = np.argmin(dist_list[idx_list_1])
                    else:
                        idx_list_2 = np.where(np.any([dist_list[idx_list_1] > dist_ub, dist_list[idx_list_1] < dist_lb], axis=0))[0]
                        idx_rep = np.argmax(np.max([dist_list[idx_list_1]-dist_ub, dist_lb-dist_list[idx_list_1]], axis=0))
                    if len(idx_list_2) == 0:
                        continue

                    consist_idx_list = idx_list[idx_list_1[idx_list_2], :]
                    rep_idx = idx_list[idx_list_1[idx_rep], :]
                    XLMS_struct_data[pair][state_idx][cluster_idx] = [rep_idx, consist_idx_list]
                    #print(f'Pair: {pair}, state_idx: {state_idx}, cluster_idx: {cluster_idx}, rep_idx: {rep_idx}, consist_idx_list: {len(consist_idx_list)}')

                    for idx in consist_idx_list:
                        dtrajs_MS[idx[0],idx[1]].append(pair)
        ##############################################################################
        self.logger.debug(f'dtrajs_MS: {dtrajs_MS.shape}')


        ##############################################################################
        # Group consistency data
        self.logger.info(f'\nGrouping consistency data based on consistency signals')
        consist_signal_dict = {}
        for i, d in enumerate(dtrajs_MS): # list of lists of consistency signals in each frame

            for j, dd in enumerate(d): # list of consistency signals in frame idx j

                # print(f'MSM traj idx: {i}, frame idx: {j}, consistency signal: {dd}')
                cluster_idxs = dtrajs_cluster_idx[i,j]
                # print(f'Cluster idx: {cluster_idxs}')

                if str(dd) not in consist_signal_dict.keys():
                    consist_signal_dict[str(dd)] = {cluster_idxs: [[i,j]]}

                elif cluster_idxs not in consist_signal_dict[str(dd)].keys():
                    consist_signal_dict[str(dd)][cluster_idxs] = [[i,j]]

                else:
                    consist_signal_dict[str(dd)][cluster_idxs].append([i,j])

        Num_struct_list = [np.sum(np.array([len(vv) for vv in v.values()])) for k,v in consist_signal_dict.items()]
        sort_idx = np.argsort(-np.array(Num_struct_list, dtype=int))
        sorted_consist_signal_list = [list(consist_signal_dict.keys())[idx] for idx in sort_idx]
        # print(f'sorted_consist_signal_list: {sorted_consist_signal_list}')
        sorted_consist_signal_dict = {}
        for k in sorted_consist_signal_list:
            kk_list = sorted(list(consist_signal_dict[k].keys()))
            sorted_consist_signal_dict[k] = {kk: consist_signal_dict[k][kk] for kk in kk_list}
            
            # for a, b in sorted_consist_signal_dict[k].items():
            #     print(f'\n{k}, {a}, {b}')
        ##############################################################################



        ##############################################################################
        # Group based on metastable states, then consistensy
        self.logger.info(f'\nGrouping consistency data based on metastable states and consistency signals')
        group_dict = {}
        for state_id in self.state_idx_list:
            # print(f'Processing state {state_id}...')
            group_dict[state_id] = {}

            idx_list = np.array(np.where(meta_dtrajs_last == state_id)).T

            for idx_list_0 in idx_list:

                [i, j] = list(idx_list_0)
                dd = dtrajs_MS[i, j]
                # print(f'Processing idx: {i}, {j}, dd: {dd}')
                cluster_idxs = dtrajs_cluster_idx[i,j]
                # print(f'Cluster idx: {cluster_idxs}')

                if str(dd) not in group_dict[state_id].keys():
                    group_dict[state_id][str(dd)] = {cluster_idxs: [[i,j]]}

                elif cluster_idxs not in group_dict[state_id][str(dd)].keys():
                    group_dict[state_id][str(dd)][cluster_idxs] = [[i,j]]

                else:
                    group_dict[state_id][str(dd)][cluster_idxs].append([i,j])

            # sort based on population
            Num_struct_list = [np.sum(np.array([len(vv) for vv in v.values()])) for k,v in group_dict[state_id].items()]
            sort_idx = np.argsort(-np.array(Num_struct_list, dtype=int))
            sorted_list = [list(group_dict[state_id].keys())[idx] for idx in sort_idx]
            new_dict = {}
            for k in sorted_list:
                kk_list = sorted(list(group_dict[state_id][k].keys()))
                new_dict[k] = {kk: group_dict[state_id][k][kk] for kk in kk_list}
            group_dict[state_id] = new_dict

            # print(f'\nGrouped consistency data for state {state_id}:')
            # for k,v in group_dict[state_id].items():
            #     print(f'{k}: {v}')
        ##############################################################################

        ##############################################################################
        # Load Q list
        self.logger.info(f'Loading G and Q data from {self.OPpath}')
        Q_list, G_list = self.load_OP(start=-n_analysis_frames)
        self.logger.info(f'Loaded G and Q data from {self.OPpath}')
        ##############################################################################

        ##############################################################################
        # Get representative structures
        self.logger.info(f'\nGetting representative structures for each group...')
        rep_group_dict = {}
        for state_id in self.state_idx_list:
            # print(f'Processing state {state_id}...')

            rep_group_dict[state_id] = {}
            for k in group_dict[state_id].keys():

                if k == '[]':
                    continue

                # print(f'    Processing consistent signal k={k}...')
                rep_group_dict[state_id][k] = {}
                for kk in group_dict[state_id][k].keys():
                    # print(f'        Processing entanglement ID kk={kk}...')
                    idx_list = np.array(group_dict[state_id][k][kk])
                    # print(f'            idx_list: {idx_list}')

                    # Max micro-states probability
                    micro_prob = meta_dist[state_id][micro_dtrajs_last[idx_list[:,0], idx_list[:,1]]]
                    max_idx = np.where(micro_prob == np.max(micro_prob))[0]
                    max_idx_list = idx_list[max_idx,:]
                    # print(f'            max_idx_list: {max_idx_list}')

                    # Max Q
                    Q_list_0 = Q_list[max_idx_list[:,0], max_idx_list[:,1]]
                    max_idx = np.where(Q_list_0 == np.max(Q_list_0))[0]
                    max_idx_list = max_idx_list[max_idx,:]
                    # print(f'            max_idx_list after Q: {max_idx_list}')
                    
                    # Max G
                    G_list_0 = G_list[max_idx_list[:,0], max_idx_list[:,1]]
                    max_idx = np.where(G_list_0 == np.max(G_list_0))[0]
                    max_idx_list = max_idx_list[max_idx,:]
                    # print(f'            max_idx_list after G: {max_idx_list}')
                    
                    [rep_traj_idx, rep_frame_idx] = max_idx_list[0,:]
                    # print(f'            Representative structure: MSM traj idx {rep_traj_idx} -> {idx2traj[rep_traj_idx]}, frame {rep_frame_idx}')
                    rep_group_dict[state_id][k][kk] = [rep_traj_idx, rep_frame_idx]
        ##############################################################################


        ##############################################################################
        # Save data
        consist_signal_struct_data_outfile = os.path.join(self.outdir, 'consist_signal_struct_data.npz')
        np.savez(consist_signal_struct_data_outfile,
                n_analysis_frames = n_analysis_frames,
                total_traj_num_frames = total_traj_num_frames,
                LIPMS_consist_data=LIPMS_consist_data,
                XLMS_consist_data=XLMS_consist_data,
                LIPMS_struct_data=LIPMS_struct_data,
                XLMS_struct_data=XLMS_struct_data,
                dtrajs_MS=dtrajs_MS,
                sorted_consist_signal_dict=sorted_consist_signal_dict,
                group_dict=group_dict,
                rep_group_dict=rep_group_dict,)
        self.logger.debug(f'SAVED: {consist_signal_struct_data_outfile}')

        # Save info to excel
        self.logger.info(f'Saving info to excel...')
        df_list = []
        sheet_name_list = []
        df_data = []
        for k in sorted_consist_signal_dict.keys():
            k_list = eval(k)
            if k == '[]':
                continue
            k_str = ', '.join(k_list)
            
            for kk in sorted_consist_signal_dict[k].keys():
                kk_list = eval(sorted_chg_ent_structure_keyword_list[kk])
                if len(kk_list) == 0:
                    continue
                kk_str = ', '.join([str(i+1) for i in kk_list])
                num = len(sorted_consist_signal_dict[k][kk])
                df_data.append([k_str, kk_str, len(k_list), num])
        df = pd.DataFrame(df_data, columns=['Consistent signals', 'IDs of Changes in Entanglements', 'Number of consistent signals', 'Number of Structures'])
        df_sorted = df.sort_values(by=['Number of consistent signals', 'Consistent signals', 'Number of Structures'], ascending=[False, True, False])
        df_list.append(df_sorted)
        sheet_name_list.append('Total')

        for state_id in self.state_idx_list:
            sheet_name_list.append('State %d'%(state_id+1))
            df_data = []

            for k in rep_group_dict[state_id].keys():
                k_list = eval(k)
                k_str = ', '.join(k_list)
                #print(f'\nProcessing state {state_id}, consistent signal {k_str}')

                for kk in rep_group_dict[state_id][k].keys():

                    kk_list = eval(sorted_chg_ent_structure_keyword_list[kk])
                    if len(kk_list) == 0:
                        continue
                    kk_str = ', '.join([str(i+1) for i in kk_list])
                    [traj_idx, frame_idx] = rep_group_dict[state_id][k][kk]
                    #print(f'traj_idx: {traj_idx} | frame_idx: {frame_idx}')

                    traj_frame_idx = total_traj_num_frames - n_analysis_frames + frame_idx
                    #print(f'Processing state {state_id}, consistent signal {k_str}, entanglement ID {kk_str}, traj {traj_idx+1}, frame {traj_frame_idx}')

                    num = len(group_dict[state_id][k][kk])
                    micro_prob = meta_dist[state_id][micro_dtrajs_last[traj_idx, frame_idx]]
                    # rel_change_list = calc_rel_change(traj_idx, frame_idx)
                    traj = idx2traj[traj_idx]
                    df_data.append([k_str, kk_str, len(k_list), num, str([traj, traj_frame_idx+1]), micro_prob, Q_list[traj_idx, frame_idx], G_list[traj_idx, frame_idx] ])

            df = pd.DataFrame(df_data, columns=['Consistent signals', 'IDs of Changes in Entanglements', 'Number of consistent signals', 'Number of Structures', 'Representative Structure (Traj #, Frame #)', 'Prob', 'Q', 'G'])
            df_sorted = df.sort_values(by=['Number of consistent signals', 'Prob', 'Q', 'G'], ascending=[False, False, False, False])
            df_list.append(df_sorted)
            
        # Creating Excel Writer Object from Pandas  
        Consistent_structures_v8_outfile = os.path.join(self.outdir, 'Consistent_structures_v8.xlsx')
        with pd.ExcelWriter(Consistent_structures_v8_outfile, engine="openpyxl") as writer:
            workbook=writer.book
            for idx_0, df in enumerate(df_list):
                worksheet=workbook.create_sheet(sheet_name_list[idx_0])
                writer.sheets[sheet_name_list[idx_0]] = worksheet
                df.to_excel(writer, sheet_name=sheet_name_list[idx_0], index=False)
        self.logger.debug(f'SAVED: {Consistent_structures_v8_outfile}')
        ##############################################################################
        

        ##############################################################################
        # Create visualization
        self.logger.info(f'Create visualizations...')
        if self.dist_data_file is not None and os.path.exists(self.dist_data_file):
            self.logger.info(f'Loading dist_traj_list: {self.dist_data_file}')
            dist_traj_list = np.load(self.dist_data_file, allow_pickle=True)[:,-n_analysis_frames:]
            self.logger.info(f'Loaded distance data: {dist_traj_list.shape}')
            dist_traj_list = [v for i, v in enumerate(dist_traj_list) if i not in np.asarray(self.rm_traj_list) - 1]
            dist_traj_list = np.array(dist_traj_list)
            self.logger.info(f'dist_traj_list after removal of mirror images: {dist_traj_list.shape}')
        elif self.xp_dir is not None and os.path.isdir(self.xp_dir):
            self.logger.info(f'Building sparse dist_traj_list from XP files in: {self.xp_dir}')
            dist_traj_list = np.empty((len(idx2traj), n_analysis_frames), dtype=object)
            dist_traj_list[:] = None

            needed_frames_by_traj = {}
            for state_id in self.state_idx_list:
                for k in rep_group_dict[state_id].keys():
                    for kk in rep_group_dict[state_id][k].keys():
                        traj_idx, frame_idx = rep_group_dict[state_id][k][kk]
                        needed_frames_by_traj.setdefault(int(traj_idx), set()).add(int(frame_idx))

            for traj_idx, needed_frame_idx in needed_frames_by_traj.items():
                traj_num = int(idx2traj[traj_idx])
                fpath = os.path.join(self.xp_dir, f'{self.ID}_Traj{traj_num}.XP')
                if not os.path.exists(fpath):
                    self.logger.warning(f'Missing XP file for sparse visualization load: {fpath}')
                    continue

                df = pd.read_csv(
                    fpath,
                    sep='\t',
                    usecols=['Frame', 'Atom1', 'Atom2', 'Euclidean Distance', 'SASD'],
                    dtype={
                        'Frame': np.int32,
                        'SASD': np.float32,
                        'Euclidean Distance': np.float32,
                        'Atom1': 'string',
                        'Atom2': 'string',
                    },
                )
                if df.empty:
                    continue

                frame_values = np.sort(df['Frame'].unique())
                frame_to_idx = {int(f): idx for idx, f in enumerate(frame_values)}
                df['frame_idx'] = df['Frame'].map(frame_to_idx)
                df = df[df['frame_idx'].isin(needed_frame_idx)]
                if df.empty:
                    continue

                atom1_parts = df['Atom1'].str.split('-', expand=True)
                atom2_parts = df['Atom2'].str.split('-', expand=True)
                if atom1_parts.shape[1] < 3 or atom2_parts.shape[1] < 3:
                    self.logger.warning(f'Malformed Atom fields in {fpath}; skipping sparse load for traj {traj_num}')
                    continue

                df['key'] = (
                    atom1_parts[1].astype(str) + '|' + atom1_parts[2].astype(str) + '-' +
                    atom2_parts[1].astype(str) + '|' + atom2_parts[2].astype(str)
                )

                for frame_idx, frame_df in df.groupby('frame_idx', sort=False):
                    frame_dict = {}
                    for _, row in frame_df.iterrows():
                        frame_dict[row['key']] = {
                            'Euclidean': float(row['Euclidean Distance']),
                            'Jwalk': float(row['SASD']),
                        }
                    dist_traj_list[traj_idx, int(frame_idx)] = frame_dict
        else:
            raise ValueError(
                'Representative structure visualization requires either dist_data_file (Jwalk.npy) or xp_dir.'
            )


        # Check if the viz_rep_struct path exists.
        if restart:
            os.makedirs('viz_rep_struct/', exist_ok=True)
            self.logger.info(f'Restart mode: preserving existing viz_rep_struct directory')
        else:
            if os.path.exists('viz_rep_struct'):
                self.logger.info(f'viz_rep_struct exists and will be removed')
                os.system('rm -rf viz_rep_struct/')
            os.system('mkdir viz_rep_struct/')
        os.chdir('viz_rep_struct/')

        if os.path.isdir(self.AAdcd_dir):
            AAtraj_files = glob.glob(os.path.join(self.AAdcd_dir, '*.dcd'))
        else:
            AAtraj_files = glob.glob(self.AAdcd_dir)
        if len(AAtraj_files) == 0:
            raise ValueError(
                f'No AA trajectory files found. AAdcd_dir={self.AAdcd_dir} '
                f'(resolved count={len(AAtraj_files)})'
            )
        self.logger.info(f'AAtraj_files:\n{AAtraj_files[:10]}')
        
        wd = os.getcwd()
        for state_id in self.state_idx_list:
            state_dir = os.path.join(wd, 'State_%d'%(state_id+1))
            os.makedirs(state_dir, exist_ok=True)
            self.logger.info(f'Made {state_dir}')
            # os.chdir(state_dir)
            self.logger.info(f'Length of rep_group_dict[state_id]: {len(rep_group_dict[state_id])}')

            # PHASE 1: Pre-extract per-(traj_idx, frame_idx) data for each (state_id, k) pair
            # to avoid pickling full arrays (~1 GB) to each worker. Extract only the specific
            # slices needed for this task.
            per_frame_data_dict = {}
            for k in rep_group_dict[state_id].keys():
                per_frame_data_dict[k] = {}
                for kk in rep_group_dict[state_id][k].keys():
                    [traj_idx, frame_idx] = rep_group_dict[state_id][k][kk]
                    if (traj_idx, frame_idx) not in per_frame_data_dict[k]:
                        # Extract full rows for M_LiPMS and M_XLMS since they're indexed
                        # by signal-dependent third indices in process_k
                        per_frame_data_dict[k][(traj_idx, frame_idx)] = {
                            'Q': Q_list[traj_idx, frame_idx],
                            'G': G_list[traj_idx, frame_idx],
                            'M_LiPMS_row': M_LiPMS[traj_idx, frame_idx, :],
                            'M_XLMS_row': M_XLMS[traj_idx, frame_idx, :],
                            'rep_chg_ent_dict': rep_chg_ent_dtrajs[traj_idx, frame_idx],
                            'dist_frame_data': dist_traj_list[traj_idx, frame_idx]
                        }

            # PHASE 2: Build args_list with pre-extracted per-frame data instead of full arrays
            args_list = [
                (state_dir, state_id, k, rep_group_dict, sorted_chg_ent_structure_keyword_list, n_analysis_frames, total_traj_num_frames,
                idx2traj, AAtraj_files, self.native_AA_pdb, per_frame_data_dict[k],
                LIPMS_consist_data, XLMS_consist_data, if_backmap, pulchra_only, self.logger.name, restart)
                for k in rep_group_dict[state_id].keys()
            ]
            self.logger.info(f'Processing {len(args_list)} consistent signal groups for state {state_id}...')
            if len(args_list) == 0:
                self.logger.info(f'No consistent signal groups for state {state_id}, skipping...')
                continue
            
            # `spawn` is safer than the Linux default `fork` when worker code touches
            # mdtraj/parmed C extensions and large shared-filesystem DCDs.
            worker_count = min(self.nproc, len(args_list))
            self.logger.info(
                f'Launching representative-structure pool for state {state_id + 1} '
                f'with start_method=spawn, processes={worker_count}, maxtasksperchild=1'
            )
            ctx = multiprocessing.get_context("spawn")
            with ctx.Pool(processes=worker_count, maxtasksperchild=1) as pool:
                for completed_idx, result in enumerate(pool.imap_unordered(process_k, args_list, chunksize=1), start=1):
                    self.logger.info(
                        f'Completed consistent signal group {completed_idx}/{len(args_list)} '
                        f'for state {state_id + 1}: {result}'
                    )
        self.logger.debug('Completion of selecting rep structure')


##################################################################################################
import mdtraj as mdt  # at the top of your file
def process_k(args):
    # PHASE 3: Updated unpacking to receive pre-extracted per-frame data instead of full arrays
    (state_dir, state_id, k, rep_group_dict, sorted_chg_ent_structure_keyword_list, n_analysis_frames, total_traj_num_frames,
     idx2traj, AAtraj_files, native_AA_pdb, per_frame_data_dict,
     LIPMS_consist_data, XLMS_consist_data, if_backmap, pulchra_only, logger_name, restart) = args

    logger = logging.getLogger(logger_name)

    k_list = eval(k)
    k_str = '_'.join(k_list)
    k_str_dir = os.path.join(state_dir, k_str)
    os.makedirs(k_str_dir, exist_ok=True)
    logger.info(f'Made {k_str_dir}')
    key_order = ['type', 'code', 'native_contact', 'native_contact_residx',  'linking_value', 'crossing_resid', 'crossing_residx', 'crossing_pattern', 'gauss_linking_number', 'topoly_linking_number', 
                 'ref_native_contact', 'ref_native_contact_residx', 'ref_linking_value', 'ref_crossing_resid', 'ref_crossing_residx', 'ref_crossing_pattern', 'ref_gauss_linking_number', 'ref_topoly_linking_number']
    completed_groups = 0

    try:
        for kk in rep_group_dict[state_id][k].keys():
            kk_list = eval(sorted_chg_ent_structure_keyword_list[kk])
            if len(kk_list) == 0:
                continue

            kk_str = '_'.join([str(i+1) for i in kk_list])
            kk_str_dir = os.path.join(k_str_dir, kk_str)
            if restart and os.path.isfile(os.path.join(kk_str_dir, '.done')):
                logger.info(f'Restart: skipping completed group {kk_str_dir}')
                continue
            os.makedirs(kk_str_dir, exist_ok=True)
            logger.info(f'Made {kk_str_dir}')

            [traj_idx, frame_idx] = rep_group_dict[state_id][k][kk]
            traj = idx2traj[traj_idx]
            traj_frame_idx = total_traj_num_frames - n_analysis_frames + frame_idx

            AAtraj_file = match_pattern(AAtraj_files, f'{traj}')
            if len(AAtraj_file) != 1:
                raise ValueError(f'Found {len(AAtraj_file)} AA traj files for traj {traj}, expected 1.')

            load_start = time.time()
            logger.info(
                f'Loading AA trajectory for state {state_id + 1}, signals {k_str}, group {kk_str}: '
                f'traj={traj}, frame={traj_frame_idx + 1}, file={AAtraj_file[0]}'
            )
            state_cor = mdt.load(AAtraj_file[0], top=native_AA_pdb)[traj_frame_idx].center_coordinates().xyz * 10
            logger.info(
                f'Loaded AA trajectory for state {state_id + 1}, signals {k_str}, group {kk_str} '
                f'in {time.time() - load_start:.2f} s'
            )

            # PHASE 3: Fetch pre-extracted frame data from per_frame_data_dict
            frame_info = per_frame_data_dict[(traj_idx, frame_idx)]
            rep_chg_ent_dict = frame_info['rep_chg_ent_dict']

            rep_ent_dict = {tuple(v['code']): [] for kkk, v in rep_chg_ent_dict.items()}
            for kkk, v in rep_chg_ent_dict.items():
                v['chg_index'] = kkk
                rep_ent_dict[tuple(v['code'])].append(v)

            viz_start = time.time()
            logger.info(
                f'Generating representative visualization for state {state_id + 1}, '
                f'signals {k_str}, group {kk_str}'
            )
            gen_state_visualizion(state_id, kk_str, kk_str_dir, native_AA_pdb, state_cor, native_AA_pdb, rep_ent_dict,
                          logger,
                                  if_backmap=if_backmap, pulchra_only=pulchra_only, exp_signal_str=k_str)
            logger.info(
                f'Finished representative visualization for state {state_id + 1}, '
                f'signals {k_str}, group {kk_str} in {time.time() - viz_start:.2f} s'
            )

            # PHASE 3: Use pre-extracted Q and G from frame_info
            Q = frame_info['Q']
            G = frame_info['G']
            info_file = os.path.join(kk_str_dir, 'info.txt')
            with open(info_file, 'w') as f:
                f.write('State #%d\n' % (state_id + 1))
                f.write('Q: %f\n' % (Q))
                f.write('G: %f\n' % (G))
                f.write('Consistent experiment signals: %s\n' % k)
                f.write('Changes in entanglement cluster IDs: %s\n' % [i+1 for i in kk_list])
                f.write('Trajectory number: %d\n' % (traj))
                f.write('Frame number: %d\n' % (traj_frame_idx + 1))
                f.write('%s\n' % ('-' * 64))
                for signal in k_list:
                    if signal in LIPMS_consist_data.keys():
                        # PHASE 3: Use pre-extracted M_LiPMS row, index by signal-dependent third dimension
                        M = frame_info['M_LiPMS_row'][LIPMS_consist_data[signal][4]]
                        f.write('%s: M: %.4f\n' % (signal, M))
                    elif signal in XLMS_consist_data.keys():
                        # PHASE 3: Use pre-extracted M_XLMS row, index by signal-dependent third dimension
                        M = frame_info['M_XLMS_row'][XLMS_consist_data[signal][4]]
                        requested_key = '%s|A-%s|A' % (signal.split('-')[0][1:], signal.split('-')[1][1:])
                        # PHASE 3: Use pre-extracted dist frame data
                        frame_data = frame_info['dist_frame_data']
                        if frame_data is None or requested_key not in frame_data:
                            logger.warning(
                                f'Missing XL-MS key for signal {signal}. '
                                f'Traj {idx2traj[traj_idx]}, frame {frame_idx}, '
                                f'MSM indices [{traj_idx}, {frame_idx}]. '
                                f'Requested key: {requested_key}. '
                                f'Available keys in frame_dict: {list(frame_data.keys()) if frame_data is not None else "None"}. '
                                f'Skipping this signal in info file.'
                            )
                            f.write('%s: M: %.4f Jwalk: N/A Euclidean: N/A (key not in XP data)\n' % (signal, M))
                        else:
                            d_data = frame_data[requested_key]
                            f.write('%s: M: %.4f Jwalk: %.4f Euclidean: %.4f\n' % (signal, M, d_data['Jwalk'], d_data['Euclidean']))
                f.write('%s\n' % ('-' * 64))
                for kkk, v in rep_chg_ent_dict.items():
                    f.write('%d:\n' % (kkk + 1))
                    for kkkk in key_order:
                        f.write(' ' * 4 + '%s: %s\n' % (kkkk, v[kkkk]))
            logger.debug(f'SAVED: {info_file}')
            with open(os.path.join(kk_str_dir, '.done'), 'w'):
                pass
            completed_groups += 1
    except Exception:
        logger.exception(
            f'Representative-structure worker failed for state {state_id + 1}, signals {k_str}'
        )
        raise

    return f'{k_str} ({completed_groups} subgroup(s))'
##############################################################################

##############################################################################
def match_pattern(strings, user_substring):
    # Build the regex pattern: non-digit (\D), user substring, then underscore
    pattern = re.compile(rf"\D{re.escape(user_substring)}_")

    # Filter the list to only those matching the pattern
    matches = [s for s in strings if pattern.search(s)]
    return matches
##############################################################################

##############################################################################
def gen_state_visualizion(state_id, ent_id, kk_str_dir, psf, state_cor, native_AA_pdb, rep_ent_dict, logger, if_backmap=True, pulchra_only=False, exp_signal_str=''):
    def idx2sel(idx_list):
        if len(idx_list) == 0:
            return ''
        else:
            sel = 'index'
            idx_0 = idx_list[0]
            idx_1 = idx_list[0]
            sel_0 = ' %d'%idx_0
            for i in range(1, len(idx_list)):
                if idx_list[i] == idx_list[i-1] + 1:
                    idx_1 = idx_list[i]
                else:
                    if idx_1 > idx_0:
                        sel_0 += ' to %d'%idx_1
                    sel += sel_0
                    idx_0 = idx_list[i]
                    idx_1 = idx_list[i]
                    sel_0 = ' %d'%idx_0
            if idx_1 > idx_0:
                sel_0 += ' to %d'%idx_1
            sel += sel_0
            return sel

    AA_name_list = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                    'HIE', 'HID', 'HIP']
    protein_colorid_list = [12, 15]
    loop_colorid_list = [4, 1]
    thread_colorid_list = [7, 0]
    nc_colorid_list = [3, 3]
    crossing_colorid_list = [8, 8]
    LIP_colorid_list = [13, 10]
    XL_colorid_list = [13, 10]
    thread_cutoff=3
    terminal_cutoff=3
    
    logger.info('Generate visualization of state %d'%(state_id + 1))
    logger.debug(f'ent_id: {ent_id}')
    logger.debug(f'kk_str_dir: {kk_str_dir}')
    logger.debug(f'psf: {psf}')
    logger.debug(f'state_cor: {state_cor}')
    logger.debug(f'native_AA_pdb: {native_AA_pdb}')
    logger.debug(f'if_backmap: {if_backmap}')
    logger.debug(f'pulchra_only: {pulchra_only}')
    logger.debug(f'exp_signal_str: {exp_signal_str}')
    logger.info(f'rep_ent_dict:')
    for k, v in rep_ent_dict.items():
        logger.info(f'k:\n{k}')
        logger.info(f'v:\n{v} {len(v)}')

    logger.info(f'Loading structure template from {psf}')
    struct = pmd.load_file(psf)
    struct.coordinates = state_cor

    # backmap
    pdb_path = os.path.join(kk_str_dir, f'state_{state_id + 1}.pdb')
    if if_backmap:
        if pulchra_only:
            pulchra_only = '1'
        else:
            pulchra_only = '0'
        temp_pdb_path = os.path.join(kk_str_dir, 'tmp.pdb')
        logger.info(f'Writing temporary backmapping PDB {temp_pdb_path}')
        struct.save(temp_pdb_path, overwrite=True)
        os.system('backmap.py -i '+native_AA_pdb+' -c '+temp_pdb_path+' -p '+pulchra_only)
        os.system('mv tmp_rebuilt.pdb '+pdb_path)
        os.system('rm -f '+temp_pdb_path)
        os.system('rm -rf ./rebuild_tmp/')
    else:
        logger.info(f'Writing representative PDB {pdb_path}')
        struct.save(pdb_path, overwrite=True)

    logger.info(f'Loading reference structure from {native_AA_pdb}')
    ref_struct = pmd.load_file(native_AA_pdb)
    logger.info(f'Loading representative structure from {pdb_path}')
    current_struct = pmd.load_file(pdb_path)
    
    # parse exp_signal_str
    if exp_signal_str == '':
        exp_signal_list = []
    else:
        exp_signal_list = exp_signal_str.strip().split('_')
        exp_signal_list = [es.strip().split('-') for es in exp_signal_list]
        exp_signal_list = [[int(ees[1:]) for ees in es] for es in exp_signal_list]

    ##############################################
    ## no change of entaglement
    if len(list(rep_ent_dict.keys())) == 0:
        vmd_outfile = os.path.join(kk_str_dir, f'vmd_s{state_id}_none.tcl')
        f = open(vmd_outfile, 'w')
        f.write('# Entanglement type: no change\n')
        f.write('''display rendermode GLSL
display projection Orthographic
axes location off

color Display {Background} white

mol new ./'''+(pdb_path)+''' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID '''+str(protein_colorid_list[1])+'''
mol selection {all}
mol material AOChalky
mol addrep top
''')
        f.close()
        logger.debug(f'SAVED: {vmd_outfile}')

    ##############################################
    ## Create vmd script for each type of change
    for ent_code, rep_ent_list in rep_ent_dict.items():
        # print(f'ent_code:\n{ent_code}')
        # print(f'rep_ent_list:\n{rep_ent_list}')
        pmd_struct_list = [ref_struct, current_struct]
        struct_dir_list = [native_AA_pdb, pdb_path]
        key_prefix_list = ['ref_', '']
        repres_list = ['', '']
        align_sel_list = ['', '']
        

        for chg_ent_fingerprint_idx, chg_ent_fingerprint in enumerate(rep_ent_list):
            # print(f'\nProcessing change of entanglement fingerprint {chg_ent_fingerprint_idx}:\n{chg_ent_fingerprint}')
            vmd_script = '''# Entanglement type: '''+str(chg_ent_fingerprint['type'])+'''
package require topotools
display rendermode GLSL
display projection Orthographic
axes location off

color Display {Background} white
color Labels {Bonds} black

label textsize 0.000001

'''
            for struct_idx, pmd_struct in enumerate(pmd_struct_list):
                struct_dir = struct_dir_list[struct_idx]
                protein_colorid = protein_colorid_list[struct_idx]
                loop_colorid = loop_colorid_list[struct_idx]
                thread_colorid = thread_colorid_list[struct_idx]
                nc_colorid = nc_colorid_list[struct_idx]
                crossing_colorid = crossing_colorid_list[struct_idx]
                LIP_colorid = LIP_colorid_list[struct_idx]
                XL_colorid = XL_colorid_list[struct_idx]
                key_prefix = key_prefix_list[struct_idx]
                # print(f'Processing structure {struct_idx}: {struct_dir}')
                # print(f'  Protein color ID: {protein_colorid}')
                # print(f'  Loop color ID: {loop_colorid}')
                # print(f'  Thread color ID: {thread_colorid}')
                # print(f'  NC color ID: {nc_colorid}')
                # print(f'  Crossing color ID: {crossing_colorid}')
                # print(f'  LIP color ID: {LIP_colorid}')
                # print(f'  XL color ID: {XL_colorid}')

                # Clean ligands
                clean_sel_idx = np.zeros(len(pmd_struct.atoms))
                for res in pmd_struct.residues:
                    if res.name in AA_name_list:
                        for atm in res.atoms:
                            clean_sel_idx[atm.idx] = 1
                pmd_clean_struct = pmd_struct[clean_sel_idx]
                clean_idx_to_idx = np.where(clean_sel_idx == 1)[0]

                # vmd selection string for protein
                idx_list = []
                for res in pmd_struct.residues:
                    if res.name in AA_name_list:
                        idx_list += [atm.idx for atm in res.atoms]
                vmd_sel = idx2sel(idx_list)

                repres = '''mol new '''+struct_dir+''' type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
    mol delrep 0 top
    mol representation NewCartoon 0.300000 10.000000 4.100000 0
    mol color ColorID '''+str(protein_colorid)+'''
    mol selection {'''+vmd_sel+'''}
    mol material GlassBubble
    mol addrep top
    '''
                align_sel = vmd_sel
  
                nc = chg_ent_fingerprint[key_prefix+'native_contact_residx']
                chg_idx = chg_ent_fingerprint['chg_index']

                idx_list = []
                for res in pmd_clean_struct.residues:
                    if res.idx in nc:
                        idx_list += [atm.idx for atm in res.atoms if atm.name == 'CA']
                nc_sel = idx2sel(clean_idx_to_idx[idx_list])
                # print(f'  Native contact residx: {nc}, selection: {nc_sel}')

                idx_list = []
                for res in pmd_clean_struct.residues:
                    if res.idx >= nc[0] and res.idx <= nc[1]:
                        idx_list += [atm.idx for atm in res.atoms]
                loop_sel = idx2sel(clean_idx_to_idx[idx_list])
                # print(f'  Loop residx: {list(range(nc[0], nc[1]+1))}, selection: {loop_sel}')
          

                align_sel += ' and not (%s)'%loop_sel
                ref_coss_resid = chg_ent_fingerprint['ref_crossing_residx']
                cross_resid = chg_ent_fingerprint['crossing_residx']
                thread = []
                thread_sel_list = []
                for ter_idx in range(len(ref_coss_resid)):
                    thread_0 = []
                    resid_list = ref_coss_resid[ter_idx] + cross_resid[ter_idx]
                    if len(resid_list) > 0:
                        thread_0 = [np.min(resid_list)-5, np.max(resid_list)+5]
                        if ter_idx == 0:
                            thread_0[0] = np.max([thread_0[0], terminal_cutoff])
                            thread_0[1] = np.min([thread_0[1], nc[0]-thread_cutoff])
                        else:
                            thread_0[0] = np.max([thread_0[0], nc[1]+thread_cutoff])
                            thread_0[1] = np.min([thread_0[1], len(struct.atoms)-1-terminal_cutoff])
                        idx_list = []
                        for res in pmd_clean_struct.residues:
                            if res.idx >= thread_0[0] and res.idx <= thread_0[1]:
                                idx_list += [atm.idx for atm in res.atoms]
                        thread_0_sel = idx2sel(clean_idx_to_idx[idx_list])
                        thread_sel_list.append(thread_0_sel)
                        align_sel += ' and not (%s)'%thread_0_sel
                    else:
                        thread_sel_list.append('')
                    thread.append(thread_0)

                ln = chg_ent_fingerprint[key_prefix+'topoly_linking_number']
                cross = []
                for i in range(len(chg_ent_fingerprint[key_prefix+'crossing_residx'])):
                    cross.append([])
                    for j in range(len(chg_ent_fingerprint[key_prefix+'crossing_residx'][i])):
                        cross[-1].append(chg_ent_fingerprint[key_prefix+'crossing_pattern'][i][j]+str(chg_ent_fingerprint[key_prefix+'crossing_residx'][i][j]))
                repres += '# idx: native_contact_residx %s, topoly_linking_number %s, crossing_residx %s.\n'%(str(nc), str(ln), str(cross))
                repres +='''    mol representation NewCartoon 0.350000 10.000000 4.100000 0
    mol color ColorID '''+str(loop_colorid)+'''
    mol selection {'''+loop_sel+'''}
    mol material Opaque
    mol addrep top
    mol representation VDW 1.000000 12.000000
    mol color ColorID '''+str(nc_colorid)+'''
    mol selection {'''+nc_sel+'''}
    mol material Opaque
    mol addrep top
    set sel [atomselect top "'''+nc_sel+'''"]
    set idx [$sel get index]
    topo addbond [lindex $idx 0] [lindex $idx 1]
    mol representation Bonds 0.300000 12.000000
    mol color ColorID '''+str(nc_colorid)+'''
    mol selection {'''+nc_sel+'''}
    mol material Opaque
    mol addrep top
    '''
                for ter_idx, thread_resid in enumerate(thread):
                    if len(thread_resid) == 0:
                        continue
                    repres += '''mol representation NewCartoon 0.350000 10.000000 4.100000 0
    mol color ColorID '''+str(thread_colorid)+'''
    mol selection {'''+thread_sel_list[ter_idx]+'''}
    mol material Opaque
    mol addrep top
    '''
                    if len(chg_ent_fingerprint[key_prefix+'crossing_residx'][ter_idx]) > 0:
                        idx_list = []
                        for res in pmd_clean_struct.residues:
                            if res.idx in chg_ent_fingerprint[key_prefix+'crossing_residx'][ter_idx]:
                                idx_list += [atm.idx for atm in res.atoms if atm.name == 'CA']
                        crossing_sel = idx2sel(clean_idx_to_idx[idx_list])
                        repres += '''mol representation VDW 1.000000 12.000000
mol color ColorID '''+str(crossing_colorid)+'''
mol selection {'''+crossing_sel+'''}
mol material Opaque
mol addrep top
'''
                #######################################
                ## showing experimenal signal residues
                for es in exp_signal_list:
                    idx_list = []
                    for res in pmd_clean_struct.residues:
                        if res.idx+1 in es:
                            idx_list += [atm.idx for atm in res.atoms if atm.name == 'CA']
                    es_sel = idx2sel(clean_idx_to_idx[idx_list])
                    if len(es) == 1: # LiP-MS signal
                        repres += '''# LiP-MS PK site '''+str(es[0])+'''
    mol representation VDW 1.000000 12.000000
    mol color ColorID '''+str(LIP_colorid)+'''
    mol selection {'''+es_sel+'''}
    mol material Opaque
    mol addrep top
    '''
                    elif len(es) == 2: # XL-MS signal
                        repres += '''# XL-MS pair ('''+str(es[0])+''', '''+str(es[1])+''')
    mol representation VDW 1.000000 12.000000
    mol color ColorID '''+str(XL_colorid)+'''
    mol selection {'''+es_sel+'''}
    mol material Opaque
    mol addrep top
    label add Bonds '''+('%d/%d %d/%d'%(struct_idx, idx_list[0], struct_idx, idx_list[1]))+'''
    '''

                    if struct_idx == 0:
                        repres += '''mol representation VDW 1.000000 12.000000
        mol color Name
        mol selection {not ('''+vmd_sel+''') and not water}
        mol material Opaque
        mol addrep top
        '''
                repres_list[struct_idx] = repres
                # print(f'representation:\n{repres}')
                align_sel_list[struct_idx] = align_sel
                # print(f'align selection:\n{align_sel}')

            vmd_script += '\n'.join(repres_list)
            vmd_script += '''
set sel1 [atomselect 0 "'''+align_sel_list[0]+''' and name CA"]
set sel2 [atomselect 1 "'''+align_sel_list[1]+''' and name CA"]
set trans_mat [measure fit $sel1 $sel2]
set move_sel [atomselect 0 "all"]
$move_sel move $trans_mat
    '''
            vmd_outfile = os.path.join(kk_str_dir, f'vmd_s{state_id + 1}_e{chg_idx + 1}_n{ent_code[0]}_c{ent_code[1]}.tcl')
            f = open(vmd_outfile, 'w')
            f.write(vmd_script)
            f.close()
            logger.debug(f'SAVED: {vmd_outfile}')
           
##############################################################################
