import time, sys
from tqdm import tqdm
import multiprocessing as mp
from scipy.stats import bootstrap, kstest
import math, random
import logging
from EntDetect._logging import setup_logger
import argparse
import glob
import numpy as np
import pandas as pd
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.model_selection import train_test_split, StratifiedKFold, KFold, cross_validate, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, average_precision_score, f1_score, recall_score, precision_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import euclidean
import matplotlib.pyplot as plt
import os

try:
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter
except ImportError:
    robjects = None
    pandas2ri = None
    importr = None
    localconverter = None
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import poisson, binom, fisher_exact, chi2, norm
import scipy.stats as st
from matplotlib.ticker import MultipleLocator
from scipy.special import expit
from scipy.stats import entropy

#pd.set_option('display.max_rows', 4000)

class ProteomeLogisticRegression:
    """
    A class to handle the data analysis process including encoding, regression, and statistical tests.
    """

    ################################################################################################
    def __init__(self, dataframe_files:str, outdir:str, gene_list:str, ID:str, reg_formula:str, log_level:int=logging.INFO, logdir:str=None):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        - dataframe_files (str): Path to residue feature files. 
                                The file contains columns for each regression variable and 1 column for the response variable. 
                                The rows should equal the number of samples in the model. 
                                The uniprot ID should be in the file name and match what is in the gene list file
        - outdir (str): Path to the output directory.
        - gene_lists (str): Path to gene lists to use. 
                            The file is a single column uniprot ID with no header                            
        - ID (str): ID for output filenames.
        - reg_formula (str): Regression formula.
        """
        self.dataframe_files = dataframe_files
        self.outdir = outdir
        self.gene_list = np.atleast_1d(np.loadtxt(gene_list, dtype=str))
        self.ID = ID
        self.logger = setup_logger('ProteomeLogisticRegression', outdir=logdir if logdir is not None else outdir, ID=ID, log_level=log_level)
        self.logger.info('Initializing ProteomeLogisticRegression')
        self.logger.debug(f'gene_list: {self.gene_list} {len(self.gene_list)}')
        self.reg_formula = reg_formula
        #self.logger = self.setup_logging()
        #self.gene_list_files = glob.glob(self.gene_list)

        if not os.path.exists(f'{self.outdir}'):
            os.makedirs(f'{self.outdir}')
            self.logger.info(f'Made output directories {self.outdir}')
    ################################################################################################

    ################################################################################################
    def encode_boolean_columns(self, df: pd.DataFrame, boolean_columns: list) -> pd.DataFrame:
        """
        Encodes boolean-like columns in a DataFrame to binary 0 and 1.

        Parameters:
        - df (pd.DataFrame): The input DataFrame.
        - boolean_columns (list): A list of column names to be encoded.

        Returns:
        - pd.DataFrame: The DataFrame with encoded columns.
        """
        label_encoder = LabelEncoder()
        
        for column in boolean_columns:
            if column in df.columns:
                df[column] = label_encoder.fit_transform(df[column])
            else:
                self.logger.info(f"Column '{column}' does not exist in the DataFrame.")
        
        return df
    ################################################################################################

    ################################################################################################
    def regression(self, df, formula):
        """
        Performs quasi-binomial regression analysis on the provided DataFrame.

        Parameters:
        - df (pd.DataFrame): DataFrame containing the data for regression.
        - formula (str): The formula specifying the regression model.

        Returns:
        - table_1_df (pd.DataFrame): DataFrame containing the regression results with p-values.
        """
        model = sm.GLM.from_formula(formula, family=sm.families.Binomial(), data=df)
        #model = smf.logit(formula=formula, data=df)
        result = model.fit()

        # Get the cov_params
        cov_matrix = result.cov_params()
        self.logger.info(f'cov_matrix:\n{cov_matrix}')
    
        # Get the coefficients
        self.logger.info("Coefficients:")
        coefficients = {'A': 0}
        for k,v in result.params.items():
            if 'AA' in k:
                k = k.replace('AA[T.', '').replace(']', '')
            coefficients[k] = v
        for k,v in coefficients.items():
            self.logger.info(f'{k}: {v}')
        self.coefficients = coefficients

        ## recalculate the pvalue to add more digits as statsmodels truncates it to 0 if it is below 0.0001 for some reason. 
        self.logger.debug(result.summary())
        table = result.summary().tables[1]
        table_df = pd.DataFrame(table.data[1:], columns=table.data[0])
        pvalues = []
        for z in table_df['z']:
            z = float(z)
            if z < 0:
                p = st.norm.cdf(z)
            else:
                p = 1 - st.norm.cdf(z)
            pvalues += [p*2]

        table_df['P>|z|'] = pvalues
        table_df = table_df.rename(columns={"": "var"})
        return table_df, cov_matrix
    ################################################################################################

    ################################################################################################
    def load_data(self, sep:str, reg_var:list, response_var:str, var2binarize:list, mask_column:str):
        """
        Loads the residue feature files and filters the data for analysis.
        
        Parameters:
        - sep (str): The separator used in the CSV files.
        - reg_var (list): List of regression variables to include in the analysis.
        - response_var (str): The response variable for the regression.
        - var2binarize (list): List of variables to binarize. Best not to use booleans and convert to 0/1
        - mask_column (list): Column header to use for masking the data. Should be a column containing null values for samples to exclude. 
        """

        self.data = pd.DataFrame()
        self.n = 0

        if os.path.isfile(self.dataframe_files):
            self.logger.info(f'Loading single design matrix file: {self.dataframe_files}')
            self.data = pd.read_csv(self.dataframe_files, sep=sep, usecols=reg_var+[response_var]+[mask_column, 'gene'])
            self.data = self.data[self.data['gene'].isin(self.gene_list)]
            self.n = len(self.data['gene'].dropna().unique())
        else:
            files = glob.glob(os.path.join(self.dataframe_files, '*'))
            #files = [f for f in files if any(s in f for s in self.gene_list)] # get only those files in the gene list
            self.logger.info(f'Number of files: {len(files)}')

            for i, gene in enumerate(self.gene_list):
                gene_resFeat = [f for f in files if gene in f]
                if len(gene_resFeat) == 0:
                    self.logger.warning(f"WARNING: No residue feature file found for gene {gene}")
                    continue
                elif len(gene_resFeat) > 1:
                    self.logger.warning(f"WARNING: More than 1 residue feature file found for gene {gene}")
                    continue
                gene_resFeat_file = gene_resFeat[0]
                #print(f'gene_resFeat_file: {gene_resFeat_file} {i}')
                if len(self.data) == 0:
                    self.data = pd.read_csv(gene_resFeat_file, sep=sep, usecols=reg_var+[response_var]+[mask_column, 'gene'])
                    self.n += 1
                else:
                    self.data = pd.concat((self.data, pd.read_csv(gene_resFeat_file, sep=sep, usecols=reg_var+[response_var]+[mask_column, 'gene'])))
                    self.n += 1

        #self.data = self.data[self.data['gene'].isin(self.reg_genes)]
        self.data = self.data[self.data['AA'] != 'NC']
        self.data = self.data[self.data[mask_column].notna()]
        self.data = self.data[self.data['AA'].notna()]
        self.data = self.data.reset_index()
        self.logger.info(f'Loaded Regression DataFrame:\n{self.data}')
        self.data = self.encode_boolean_columns(self.data, boolean_columns=var2binarize)
        self.data = self.data[[response_var]+reg_var]
        self.logger.info(f'Loaded Regression DataFrame:\n{self.data}')
        #print(f"Data loaded and filtered. Number of unique genes: {len(self.data['gene'].unique())}")
    ################################################################################################

    ################################################################################################
    def run(self, ):
        """
        Orchestrates the workflow by loading data, performing regression, and saving results.
        """

        # Perform regression
        reg, cov_matrix = self.regression(self.data, self.reg_formula)
        reg['coef'] = reg['coef'].astype(float)
        reg['OR'] = np.exp(reg['coef'].astype(float))
        reg['std err'] = reg['std err'].astype(float)
        reg['z'] = reg['z'].astype(float)
        reg['P>|z|'] = reg['P>|z|'].astype(float)
        reg['[0.025'] = np.exp(reg['[0.025'].astype(float))
        reg['0.975]'] = np.exp(reg['0.975]'].astype(float))
        reg['ID'] = self.ID
        reg['n'] = self.n 
        self.logger.info(f'Regression Results:\n{reg.to_string()}')

        return reg
    ################################################################################################

    # ################################################################################################
    # def plot_regression_results(self, data:pd.DataFrame, reg_df: pd.DataFrame, outfile: str, title: str, reg_var:str, response_var:str):
    #     """
    #     Plot the regression results as a single single plot of the probability of the response variable vs the regression variables.
    #     P(x) = 1 / (1 + exp(-(B0 + B1*X1 + B2*X2 + ... + Bn*Xn)))
    #     """
    #     import matplotlib.pyplot as plt
    #     print(data)
    #     print(reg_df)

    #     ## get the coefficients into an easy to read dictionary
    #     coefficients = {}
    #     for rowi, row in reg_df.iterrows():
    #         if 'AA' in row['var']:
    #             row['var'] = row['var'].replace('AA[T.', '').replace(']', '')
    #         coefficients[row['var']] = row['coef']
    #     print(f'coefficients: {coefficients}')

    #     ## plot cut_C_Rall vs region and then fit a curve to it where cut_C_Rall ~ 1 / (1 + exp(-(B0 + B1*X1 + B2*X2 + ... + Bn*Xn)))
    #     ## where X1 is the region and X2..Xn is the AA terms
    #     df = data.copy()
    #     # One-hot encode the AA column (exclude AA not in coefficients)
    #     aa_keys = set(coefficients.keys()) - {'Intercept', 'region'}
    #     df = df[df['AA'].isin(aa_keys)]  # filter only known AA
    #     aa_dummies = pd.get_dummies(df['AA'])
    #     df = pd.concat([df, aa_dummies], axis=1)

    #     # Add coefficient-weighted terms
    #     df['linear_predictor'] = coefficients['Intercept'] + coefficients['region'] * df['region']
    #     for aa in aa_keys:
    #         if aa in df.columns:
    #             df['linear_predictor'] += coefficients[aa] * df[aa]

    #     # Compute predicted probability
    #     df['predicted'] = expit(df['linear_predictor'])
    #     print(df)

    #     # Aggregate by region for plotting
    #     agg = df.groupby('region')[['cut_C_Rall', 'predicted']].mean().reset_index()

    #     # Plot
    #     plt.figure(figsize=(8, 6))
    #     plt.plot(agg['region'], agg['cut_C_Rall'], label='Observed', marker='o')
    #     plt.plot(agg['region'], agg['predicted'], label='Predicted (logistic fit)', marker='x')
    #     plt.xlabel('Region')
    #     plt.ylabel('Probability of cut_C_Rall = 1')
    #     plt.title('Observed vs Predicted cut_C_Rall by Region')
    #     plt.legend()
    #     plt.grid(True)
    #     plt.tight_layout()


    #     plt.savefig(outfile)
    #     print(f'SAVED: {outfile}')
    # ################################################################################################

#########################################################################################################################
#########################################################################################################################
class MonteCarlo:

    """
    A class to handle the data analysis process including encoding, regression, and statistical tests.
    """

    def __init__(self, dataframe_files:str, outdir:str, gene_list:str, ID:str, reg_formula:str, response_var:str, test_var:str, random:bool, n_groups:int, steps:int, C1:float, C2:float, beta:float, linearT:bool, log_level:int=logging.INFO, logdir:str=None):
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        - dataframe_files (str): Path to residue feature files.
        - outdir (str): Path to the output directory.
        - gene_lists (str): Path to gene lists to use.
        - ID (str): ID for output filenames.
        - reg_formula (str): Regression formula.
        - response_var (str): The response variable for the regression.
        - test_var (str): The test variable for the regression.
        - random (bool): Whether to randomize the data.
        - n_groups (int): Number of groups for the Monte Carlo simulation.
        - steps (int): Number of steps for the Monte Carlo simulation.
        - C1 (float): Coefficient for the energy function.
        - C2 (float): Coefficient for the energy function.
        - beta (float): Initial temperature for the simulated annealing.
        - linearT (bool): Whether to use linear temperature scaling.
        """
        self.dataframe_files = dataframe_files
        self.outdir = outdir
        self.gene_list = np.atleast_1d(np.loadtxt(gene_list, dtype=str))
        self.num_genes = len(self.gene_list)
        self.logger = setup_logger('MonteCarlo', outdir=logdir if logdir is not None else outdir, ID=ID, log_level=log_level)
        self.logger.debug(f'gene_list: {self.gene_list} {self.num_genes}')
        
        self.ID = ID
        self.reg_formula = reg_formula
        self.response_var = response_var
        self.test_var = test_var
        self.steps = steps
        self.C1 = C1
        self.C2 = C2
        self.beta = beta
        self.linearT = linearT
        self.data = {}
        self.n_groups = n_groups
        self.random = random

        if not os.path.exists(f'{self.outdir}'):
            os.makedirs(f'{self.outdir}')
            self.logger.info(f'Made output directories {self.outdir}')

        self.logger.info(f'{"#"*100}\nNEW RUN')

        # store the parameters in the log file
        self.logger.info(f'dataframe_files: {self.dataframe_files}')
        self.logger.info(f'outdir: {self.outdir}')
        self.logger.info(f'gene_list: {self.gene_list} {self.num_genes}')
        self.logger.info(f'ID: {self.ID}')
        self.logger.info(f'response_var: {self.response_var}')
        self.logger.info(f'test_var: {self.test_var}')
        self.logger.info(f'steps: {self.steps}')
        self.logger.info(f'C1: {self.C1}')
        self.logger.info(f'C2: {self.C2}')
        self.logger.info(f'linearT: {self.linearT}')
        self.logger.info(f'beta: {self.beta}')
        self.logger.info(f'n_groups: {self.n_groups}')
        
    ################################################################################################

    ################################################################################################
    def setup_logging(self):
        """
        Sets up the logging configuration.

        Returns:
        - logger (logging.Logger): Configured logger.
        """
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)
        return logger
    ################################################################################################

    ################################################################################################
    def encode_boolean_columns(self, df: pd.DataFrame, boolean_columns: list) -> pd.DataFrame:
        """
        Encodes boolean-like columns in a DataFrame to binary 0 and 1.

        Parameters:
        - df (pd.DataFrame): The input DataFrame.
        - boolean_columns (list): A list of column names to be encoded.

        Returns:
        - pd.DataFrame: The DataFrame with encoded columns.
        """
        label_encoder = LabelEncoder()

        for column in boolean_columns:
            if column in df.columns:
                df[column] = label_encoder.fit_transform(df[column])
            else:
                self.logger.info(f"Column '{column}' does not exist in the DataFrame.")

        return df
    ################################################################################################

    ################################################################################################
    def regression(self, df, formula, genes):
        """
        Performs quasi-binomial regression analysis on the provided DataFrame.

        Parameters:
        - df (pd.DataFrame): DataFrame containing the data for regression.
        - formula (str): The formula specifying the regression model.

        Returns:
        - table_1_df (pd.DataFrame): DataFrame containing the regression results with p-values.
        """

        model = sm.GLM.from_formula(formula, family=sm.families.Binomial(), data=df)
        result = model.fit()

        ## recalculate the pvalue to add more digits as statsmodels truncates it to 0 if it is below 0.0001 for some reason.
        table = result.summary().tables[1]
        table_df = pd.DataFrame(table.data[1:], columns=table.data[0])
        pvalues = []
        for z in table_df['z']:
            z = float(z)
            if z < 0:
                p = st.norm.cdf(z)
            else:
                p = 1 - st.norm.cdf(z)
            pvalues += [p*2]
        table_df['P>|z|'] = pvalues
        table_df = table_df.rename(columns={"": "var"})

        coef = table_df[table_df['var'] == 'region']['coef'].values[0]
        std = table_df[table_df['var'] == 'region']['std err'].values[0]

        # get size dist
        size_dist = self.prot_size[self.prot_size['gene'].isin(genes)]['prot_size'].values
        return table_df, float(coef), float(std), table_df[table_df['var'] == 'region'], size_dist
    ################################################################################################

    ################################################################################################
    def metrics(self, df, genes:list):

        ctable = pd.crosstab(df[self.response_var], df[self.test_var])
        res = fisher_exact(ctable)
        OR, pvalue = res.statistic, res.pvalue

        # get size dist
        size_dist = self.prot_size[self.prot_size['gene'].isin(genes)]['prot_size'].values

        return OR, pvalue, size_dist
    ################################################################################################

    ################################################################################################
    def load_data(self, sep:str, reg_var:list, var2binarize:list, mask_column:str, ID_column:str, Length_column:str):
        """
        Loads the residue feature files and filters the data for analysis.
        
        Parameters:
        - sep (str): The separator used in the CSV files.
        - reg_var (list): List of regression variables to include in the analysis.
        - var2binarize (list): List of variables to binarize. Best not to use booleans and convert to 0/1
        - mask_column (list): Column header to use for masking the data. Should be a column containing null values for samples to exclude. 
        - ID_column (str): Column header for the gene ID.
        """

        reg_var = list(reg_var)
        if self.test_var not in reg_var:
            self.logger.info(f'Appending test_var {self.test_var} to reg_var for MonteCarlo loading')
            reg_var.append(self.test_var)

        self.data = pd.DataFrame()
        self.n = 0
        self.prot_size = {'gene':[], 'prot_size':[]}
        if os.path.isfile(self.dataframe_files):
            self.logger.info(f'Loading single design matrix file: {self.dataframe_files}')
            self.data = pd.read_csv(self.dataframe_files, sep=sep, usecols=reg_var+[self.response_var]+[mask_column, ID_column, Length_column])
            self.data = self.data[self.data[ID_column].isin(self.gene_list)]
            self.n = len(self.data[ID_column].dropna().unique())

            size_df = self.data[[ID_column, Length_column]].dropna().drop_duplicates(subset=[ID_column])
            self.prot_size = {
                'gene': size_df[ID_column].tolist(),
                'prot_size': size_df[Length_column].tolist(),
            }
        else:
            files = glob.glob(os.path.join(self.dataframe_files, '*'))
            #files = [f for f in files if any(s in f for s in self.gene_list)] # get only those files in the gene list
            self.logger.info(f'Number of files: {len(files)}')

            for i, gene in enumerate(self.gene_list):
                gene_resFeat = [f for f in files if gene in f]
                if len(gene_resFeat) == 0:
                    self.logger.warning(f"WARNING: No residue feature file found for gene {gene}")
                    continue
                elif len(gene_resFeat) > 1:
                    self.logger.warning(f"WARNING: More than 1 residue feature file found for gene {gene}")
                    continue
                gene_resFeat_file = gene_resFeat[0]
                #print(f'gene_resFeat_file: {gene_resFeat_file} {i}')
                if len(self.data) == 0:
                    self.data = pd.read_csv(gene_resFeat_file, sep=sep, usecols=reg_var+[self.response_var]+[mask_column, ID_column, Length_column])
                    self.n += 1
                    self.prot_size['gene'] += [gene]
                    self.prot_size['prot_size'] += [self.data[Length_column].values[0]]
                else:
                    df = pd.read_csv(gene_resFeat_file, sep=sep, usecols=reg_var+[self.response_var]+[mask_column, ID_column, Length_column])
                    self.data = pd.concat((self.data, df))
                    self.n += 1
                    self.prot_size['gene'] += [gene]
                    self.prot_size['prot_size'] += [df[Length_column].values[0]]

        #self.data = self.data[self.data['gene'].isin(self.reg_genes)]
        self.data = self.data[self.data['AA'] != 'NC']
        self.data = self.data[self.data[mask_column].notna()]
        self.data = self.data[self.data['AA'].notna()]
        self.data = self.data.reset_index()
        self.logger.info(f'Loaded Regression DataFrame:\n{self.data}')
        self.logger.info(f'number of genes: {len(self.data[ID_column].unique())}')

        self.data = self.encode_boolean_columns(self.data, boolean_columns=var2binarize)
        self.data = self.data[[ID_column]+[self.response_var]+reg_var]
        self.logger.info(f'Loaded Regression DataFrame:\n{self.data}')

        self.prot_size = pd.DataFrame(self.prot_size)
        self.logger.info(f'self.prot_size: {self.prot_size}')
        #print(f"Data loaded and filtered. Number of unique genes: {len(self.data['gene'].unique())}")

        if self.random:
            self.logger.info(f'Randomizing {self.response_var} column')
            self.data[self.response_var] = self.data[self.response_var].sample(frac=1).values
            self.logger.info(f'Randomized Regression DataFrame:\n{self.data}')

    ################################################################################################

    ################################################################################################
    def run(self, encoded_df, ID_column:str):
        """
        Orchestrates the workflow by loading data, performing regression, and saving results.
        """

        ## load reference size dist
        self.ref_sizes = self.prot_size['prot_size'].values
        self.logger.info(f'ref_sizes:\n{self.ref_sizes}')

        #####################################################################################################
        ## do it randomly untill the groups ceof are in an assending order
        groups = {}
        subgroups = self.create_unique_subgroups(np.arange(len(self.gene_list)), self.n_groups)
        coefs = []
        last_step = -1
        for n, subgroup in enumerate(subgroups):
            #logging.info(n, len(subgroup), subgroup, type(subgroup))

            subgroup = np.array(subgroup, dtype=int)
            groups[n] = {}
            groups[n]['genes'] = [self.gene_list[s] for s in subgroup]

            # get the inital OR, number of cuts, and regression row
            #_, coef, std, cuts, reg_row, size_dist = self.regression(encoded_df[encoded_df['gene'].isin(groups_temp[n]['genes'])][reg_vars], self.reg_formula, groups_temp[n]['genes'])
            #response_var:str, test_var:str, genes:list
            OR, pvalue, size_dist = self.metrics(encoded_df[encoded_df[ID_column].isin(groups[n]['genes'])], groups[n]['genes'])
            self.logger.debug(f'OR: {OR} pvalue: {pvalue} size_dist: {size_dist}')

            state_size_dist_bootres = bootstrap((size_dist,) , np.mean)
            state_size_mean, state_size_lb, state_size_ub = np.mean(size_dist), state_size_dist_bootres.confidence_interval.low, state_size_dist_bootres.confidence_interval.high
            ks_stat_size = kstest(self.ref_sizes, size_dist).statistic

            E = -1*self.C1*np.log(OR) + self.C2*(ks_stat_size) 

            groups[n]['OR'] = [OR]
            groups[n]['pvalue'] = [pvalue]
            groups[n]['size_dist'] = [size_dist]
            groups[n]['psize_mean'] = [state_size_mean]
            groups[n]['psize_lb'] = [state_size_lb]
            groups[n]['psize_ub'] = [state_size_ub]
            groups[n]['step'] = [last_step]
            groups[n]['E'] = [E]
            groups[n]['ks_stat_size'] = [ks_stat_size]
            groups[n]['beta'] = [self.beta]

        ## log the starting information 
        for state in range(len(groups)):
            state_data = groups[state]
            state_OR = state_data['OR'][-1]
            state_pvalue = state_data['pvalue'][-1]
            num_genes = len(state_data['genes'])
            state_step = state_data['step'][-1]
            state_E = state_data['E'][-1]
            ks_stat_size = state_data['ks_stat_size'][-1]
            state_beta = state_data['beta'][-1]

            self.logger.info(f'STEP: {state_step} | State: {state} | OR: {state_OR} | pvalue: {state_pvalue} | num_genes: {num_genes} | E: {state_E} | ks_stat_size: {ks_stat_size} | beta: {state_beta}')
            self.logger.debug(f'STEP: {state_step} | State: {state} | OR: {state_OR} | pvalue: {state_pvalue} | num_genes: {num_genes} | E: {state_E} | ks_stat_size: {ks_stat_size} | beta: {state_beta}')


        #####################################################################################################
        #####################################################################################################
        ## Start optimizer and run for X steps using simulated annealing where 
        # the energy function for each state is defined as E(i) = -C1*OR(i) + C2*ks-test(ref_sizes, sizes))|
        # for each step we swap non-overlapping states. For exampline in a 5 state system
        # step 1 we swap 1-2, 3-4
        # step 2 we swap 2-3, 4-5
        # this constitutes one full montecarlo step
        # 
        # after every 1000 MC steps move to the next beta where beta ranges from (1/15 to 1000)
        #
        # The objective function is then M = exp(-beta*deltaE)
        # if deltaE is <=0 Accept the step
        # else
        #   get random float [0,1]
        #   M > random_float and M < 1 accept, else reject
        #
        # beta should start at 50 times the starting energy scale and decrease by half every 1000 steps till it reaches 1 (?)

        def swap_n(list1, list2, n):
            # Ensure n is not larger than the smaller list size
            n = min(n, len(list1), len(list2))
            list1 = np.array(list1)
            list2 = np.array(list2)
            
            # Randomly select n indices from both lists
            indices1 = random.sample(range(len(list1)), n)
            indices2 = random.sample(range(len(list2)), n)

            # Swap the elements at the selected indices
            for i in range(n):
                idx1 = indices1[i]
                idx2 = indices2[i]
                gene1 = list1[indices1[i]]
                gene2 = list2[indices2[i]]
                #print(i, idx1, idx2, gene1, gene2)
                list1[idx1] = gene2
                list2[idx2] = gene1

            return list1, list2

        ## create swapping scheme for each monte carlo step
        self.logger.info(f'Making paris for MC swaps of {self.n_groups} groups in {self.steps} steps')
        pairs = [(i, i+1) for i in range(0, self.n_groups - 1 )]
        self.logger.debug(f'pairs: {pairs}')
      
        reps = self.steps
        beta = self.beta
        beta_i = 0
        if self.linearT == False:
            betas = np.linspace(beta, 1000, 75)
        elif self.linearT == True:
            T = 1/beta
            Ts = np.linspace(T, 0.001, 75)
            self.logger.debug(f'Ts: {Ts}')
            betas = 1/Ts
        self.logger.debug(f'betas: {betas} wiht start beta: {betas[beta_i]}')
        logging.info(f'Starting simulations with {reps} steps and beta = {self.beta}')
        
        for step in tqdm(range(last_step + 1, reps + last_step + 1)):
            #for step in range(last_step + 1, reps + last_step + 1):
            #logging.info(f'{"#"*100}\n{"#"*100}\nSTEP: {step}')
            #print(f'{"#"*100}\n{"#"*100}\nSTEP: {step}')

            ## For each pair in the pairs to test
            for pair in pairs:
                #print(f'{"#"*100}\npair: {pair}')
                
                # Get previous energy 
                Eold = groups[pair[0]]['E'][-1] + groups[pair[1]]['E'][-1]
                #print(f'Eold: {Eold}')

                # Get old state genes
                p0_genes = groups[pair[0]]['genes']
                p1_genes = groups[pair[1]]['genes']

                # Swap n genes 
                p0_genes_prime, p1_genes_prime = swap_n(p0_genes, p1_genes, 5)

                # get new regression info
                #_, p0_coef, p0_std, p0_cuts, p0_reg_row, p0_size_dist = self.regression(encoded_df[encoded_df['gene'].isin(p0_genes_prime)][reg_vars], self.reg_formula, p0_genes_prime)
                OR, pvalue, size_dist = self.metrics(encoded_df[encoded_df[ID_column].isin(p0_genes_prime)], p0_genes_prime)
                p0_OR, p0_pvalue, p0_size_dist = self.metrics(encoded_df[encoded_df[ID_column].isin(p0_genes_prime)], p0_genes_prime)
                p0_ks_stat_size = kstest(self.ref_sizes, p0_size_dist).statistic
                Ep0 = -1*self.C1*np.log(p0_OR) + self.C2*(p0_ks_stat_size) 

                #_, p1_coef, p1_std, p1_cuts, p1_reg_row, p1_size_dist = self.regression(encoded_df[encoded_df['gene'].isin(p1_genes_prime)][reg_vars], self.reg_formula, p1_genes_prime)
                p1_OR, p1_pvalue, p1_size_dist = self.metrics(encoded_df[encoded_df[ID_column].isin(p1_genes_prime)], p1_genes_prime)
                p1_ks_stat_size = kstest(self.ref_sizes, p1_size_dist).statistic
                Ep1 = -1*self.C1*np.log(p1_OR) + self.C2*(p1_ks_stat_size)

                # Calculate new E and deltaE
                Enew = Ep0 + Ep1
                #print(f'Enew: {Enew}')
                deltaE = Enew - Eold
                #print(f'deltaE: {deltaE}')

                # Apply metropolis critera
                rand_float = random.uniform(0, 1)
                M = np.exp(-1*beta*deltaE)
                if deltaE <= 0:
                    accept_M = True
                else:
                    ## Apply metropolis critera
                    if M < 1 and M > rand_float:
                        accept_M = True
                    else:
                        accept_M = False
                #print(f'M: {M} with accept_M: {accept_M}')
                
                if accept_M:
                    groups[pair[0]]['genes'] = p0_genes_prime
                    groups[pair[1]]['genes'] = p1_genes_prime
            
            ## Get step summary after all pairs have been 
            #print(f'{"#"*100}\nStep summary for {step}')
            logstr = [f'{"#"*50}']
            for state, state_data in groups.items():
                state_genes = state_data['genes']
                #_, step_coef, step_std, step_cuts, step_reg_row, step_size_dist = self.regression(encoded_df[encoded_df['gene'].isin(state_genes)][reg_vars], self.reg_formula, state_genes)
                state_OR, state_pvalue, state_size_dist = self.metrics(encoded_df[encoded_df[ID_column].isin(state_genes)], state_genes)
                state_size_dist_bootres = bootstrap((state_size_dist,) , np.mean)
                state_size_mean, state_size_lb, state_size_ub = np.mean(state_size_dist), state_size_dist_bootres.confidence_interval.low, state_size_dist_bootres.confidence_interval.high

                ks_stat_size = kstest(self.ref_sizes, state_size_dist).statistic

                E = -1*self.C1*np.log(state_OR) + self.C2*(ks_stat_size) 
                #print(f'STEP: {step} | state {state} OR: {state_OR} pvalue: {state_pvalue} cuts: {state_cuts} size_mean: {state_size_mean} | ks_stat: {ks_stat} | E: {E}')
                logstr += [f'STEP: {step} | state {state} OR: {state_OR} pvalue: {state_pvalue} size_mean: {state_size_mean} | ks_stat_size: {ks_stat_size} | E: {E} | beta: {beta}']

                groups[state]['OR'] += [state_OR]
                groups[state]['pvalue'] += [state_pvalue]
                groups[state]['size_dist'] += [state_size_dist]
                groups[state]['psize_mean'] += [state_size_mean]
                groups[state]['psize_lb'] += [state_size_lb]
                groups[state]['psize_ub'] += [state_size_ub]
                groups[state]['step'] += [step]
                groups[state]['ks_stat_size'] += [ks_stat_size]
                groups[state]['E'] += [E]
                groups[state]['beta'] += [beta]

 
            ## check ranks
            #old_ranks = sorted(range(len(old_coefs)), key=lambda i: old_coefs[i], reverse=True)
            #new_ranks = sorted(range(len(new_coefs)), key=lambda i: new_coefs[i], reverse=True)
            #rank_cond = new_ranks == old_ranks

            # logging.info status of step to log file
            if step % 100 == 0:
                logstr += [f'{"#"*50}']
                logstr = '\n'.join(logstr)
                logging.info(logstr)


            # update beta
            if step % 750 == 0 and step > 10:
                beta_i += 1
                if beta_i < len(betas):
                    beta = betas[beta_i]
                else:
                    beta = 1000
                logging.info(f'Beta update: {beta}')


        logging.info(f'{"#*50"}Simulation complete')
        logging.info(f'{"#*50"}Final state stats')
        for state in range(len(groups)):
            state_data = groups[state]
            state_OR = state_data['OR'][-1]
            state_pvalue = state_data['pvalue'][-1]
            state_size_dist = state_data['size_dist'][-1]
            state_size_mean = state_data['psize_mean'][-1]
            state_size_lb = state_data['psize_lb'][-1]
            state_size_ub = state_data['psize_ub'][-1]
            state_ks_stat_size = state_data['ks_stat_size'][-1]
            state_E = state_data['E'][-1]
            state_beta = state_data['beta'][-1]
            logging.info(f'State: {state} | OR: {state_OR} | pvalue: {state_pvalue} | state_size_mean: {state_size_mean:.0f} ({state_size_lb:.0f}, {state_size_ub:.0f}) | state_ks_stat_size: {state_ks_stat_size} | E: {state_E} | beta: {state_beta}')
            self.logger.debug(f'State: {state} | OR: {state_OR} | pvalue: {state_pvalue} | state_size_mean: {state_size_mean:.0f} ({state_size_lb:.0f}, {state_size_ub:.0f}) | state_ks_stat_size: {state_ks_stat_size} | E: {state_E} | beta: {state_beta}')
        #####################################################################################################

        #####################################################################################################
        ## Save final results
        dfs = []
        for state, state_data in groups.items():
            #logging.info(n, len(subgroup), subgroup, type(subgroup))

            # get the inital OR, number of cuts, and regression row
            _, coef, std, reg_row, size_dist = self.regression(encoded_df[encoded_df[ID_column].isin(groups[state]['genes'])], self.reg_formula, groups[state]['genes'])
            state_size_dist_bootres = bootstrap((size_dist,) , np.mean)
            state_size_mean, state_size_lb, state_size_ub = np.mean(size_dist), state_size_dist_bootres.confidence_interval.low, state_size_dist_bootres.confidence_interval.high

            ks_stat_size = kstest(self.ref_sizes, size_dist).statistic

            E = -1*self.C1*coef + self.C2*(ks_stat_size) 

            reg_row['state'] = state
            reg_row['beta'] = beta
            reg_row['ks_stat_size'] = ks_stat_size
            reg_row['E'] = E
            reg_row['psize_mean'] = state_size_mean
            reg_row['psize_lb'] = state_size_lb
            reg_row['psize_ub'] = state_size_ub
            reg_row['step'] = step
            reg_row['OR'] = np.exp(coef)
            reg_row['OR_lb'] = np.exp(reg_row['[0.025'].astype(float))
            reg_row['OR_ub'] = np.exp(reg_row['0.975]'].astype(float))
            reg_row['ID'] = self.ID
            reg_row['n'] = len(state_data['genes'])
            dfs += [reg_row]

            ## save the final gene list for the state
            state_final_genelist_outfile = os.path.join(self.outdir, f'State{state}_final_genelist_{self.ID}.txt')
            logging.info(state_data['genes'])
            np.savetxt(state_final_genelist_outfile, list(state_data['genes']), fmt='%s')
            logging.info(f'SAVED: {state_final_genelist_outfile}')

            ## save the state data for this state
            state_df = {'state':[], 'step':[], 'OR':[], 'pvalue':[], 'psize_mean':[], 'psize_lb':[], 'psize_ub':[], 'ks_stat_size':[], 'E':[], 'beta':[]}
            for i in range(len(state_data['step'])):
                state_df['state'] += [state]
                state_df['step'] += [state_data['step'][i]]
                state_df['OR'] += [state_data['OR'][i]]
                state_df['pvalue'] += [state_data['pvalue'][i]]
                state_df['psize_mean'] += [state_data['psize_mean'][i]]
                state_df['psize_lb'] += [state_data['psize_lb'][i]]
                state_df['psize_ub'] += [state_data['psize_ub'][i]]
                state_df['ks_stat_size'] += [state_data['ks_stat_size'][i]]
                state_df['E'] += [state_data['E'][i]]
                state_df['beta'] += [state_data['beta'][i]]
            state_df = pd.DataFrame(state_df)
            #print(state_df)
            state_final_traj_outfile = os.path.join(self.outdir, f'State{state}_final_traj_{self.ID}.csv')
            state_df.to_csv(state_final_traj_outfile, index=False)
            logging.info(f'SAVED: {state_final_traj_outfile}')

        outdf = pd.concat(dfs)
        self.logger.debug(outdf)

        ## Save the final step regression data
        final_step_outfile = os.path.join(self.outdir, f'Final_step_reg_{self.ID}.csv')
        outdf.to_csv(final_step_outfile, index=False)
        self.logger.debug(f'SAVED: {final_step_outfile}')
        logging.info(f'SAVED: {final_step_outfile}')

        #####################################################################################################
    ################################################################################################

    ################################################################################################
    def create_unique_subgroups(self, array, m):
        # Shuffle the array
        np.random.shuffle(array)
        
        # Calculate the size of each subgroup
        # If you want each subgroup to have equal size:
        subgroup_size = len(array) // m
        
        # Initialize the list of subgroups
        subgroups = []
        
        # Split the array into subgroups
        for i in range(m):
            start_index = i * subgroup_size
            end_index = start_index + subgroup_size
            
            # Handle the case where the last subgroup might be larger
            if i == m - 1:
                subgroups.append(array[start_index:])
            else:
                subgroups.append(array[start_index:end_index])
        
        return subgroups
    ################################################################################################

    ################################################################################################
    def is_decending(self, array):
        # Iterate through the array and compare adjacent elements
        for i in range(len(array) - 1):
            if array[i] < array[i + 1]:
                return False
        return True
    ################################################################################################

    ################################################################################################
    def get_per_res_cuts(self, df, cutkey):
        per_res_cuts_df = {'gene':[], 'per_res_cuts':[]}
        for gene, gene_df in df.groupby('gene'):
            per_res_cuts_df['gene'] += [gene]
            per_res_cuts_df['per_res_cuts'] += [np.sum(gene_df[cutkey])/len(gene_df)]
        per_res_cuts_df = pd.DataFrame(per_res_cuts_df)
        self.logger.debug(per_res_cuts_df)
        return per_res_cuts_df
    ################################################################################################
    
#########################################################################################################################
#########################################################################################################################
class FoldingPathwayStats:

    """
    A class to analyze the folding pathway statistics resulting from the Markov State Model (MSM) of folding generated by EntDetect.clustering.MSMNonNativeEntanglementClustering
    """

    ####################################################################
    def __init__(self, outdir:str='./FoldingPathwayStats', 
                 n_window:int=200, 
                 n_traj:int=1000, 
                 state_type:str='metastablestate',
                 rm_traj_list:list=[], log_level:int=logging.INFO, logdir:str=None): 
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        - outdir (str): Path to the output directory.
        """
        self.outdir = outdir
        self.n_window = n_window
        self.n_traj = n_traj
        self.state_type = state_type
        self.rm_traj_list = [str(t) for t in rm_traj_list]  # Ensure rm_traj_list is a list of strings
        if self.state_type not in ['microstate', 'metastablestate']:
            raise ValueError(f'state_type must be [microstate | metastablestate]')

        self.logger = setup_logger('FoldingPathwayStats', outdir=logdir if logdir is not None else self.outdir, log_level=log_level)
        if not os.path.exists(f'{self.outdir}'):
            os.makedirs(f'{self.outdir}')
            self.logger.debug(f'Made output directories {self.outdir}')

        ## set delta time between frames in experimental seconds
        dt = 0.015/1000 # time step 
        nsave = 5000 # number of steps between frames
        alpha = 4331293.0 
        self.dt = dt*nsave*alpha/1e9 # in seconds
    ####################################################################

    ####################################################################
    def post_trans(self, msm_data_file:str, traj_type_col:str, traj_type_list:list, outdir:str=None, outfile_name:str=None):
        """
        The folding pathways are identified as following:

        (1) For each discrete trajectory, put the starting state at the first frame into the pathway;

        (2) Move forward along the trajectory and find the state that is different with the last state recoded in the pathway. 
            If the state has not yet been recorded in the pathway, then put it into the pathway. 
            Otherwise, cut the pathway at the first place where this state is recorded and then move forward;

        (3) Repeat step (2) until reach the end of the trajectory.

        This will yield a pathway that has no loop on the route and only records the on-pathway states for each discrete trajectory.
        """
        # print(f'Loading MSM data from {self.msm_data}')
        # msm_data = pd.read_csv(self.msm_data)
        # print(f'msm_data:\n{msm_data}')
        if not msm_data_file:
            raise ValueError('msm_data_file must be provided')
        if not traj_type_list:
            raise ValueError('traj_type_list must contain at least one trajectory type label')

        target_outdir = outdir if outdir is not None else self.outdir
        os.makedirs(target_outdir, exist_ok=True)

        self.logger.info(f'Loading MSM data from {msm_data_file}')
        msm_data = pd.read_csv(msm_data_file)
        msm_data = msm_data[~msm_data['traj'].isin(self.rm_traj_list)]
        self.logger.info(f'msm_data:\n{msm_data}')

        folding_pathways = {}
        folding_pathways_df = {traj_type_col:[], 'pathway':[], 'probability':[]}
        for traj_type in traj_type_list:
            self.logger.info(f'Processing {traj_type} trajectories')
            traj_type_msm_data = msm_data[msm_data[traj_type_col] == traj_type]
            self.logger.info(f'traj_type_msm_data:\n{traj_type_msm_data}')
                
            # Quality check that there is data
            if traj_type_msm_data.empty:
                raise ValueError(f"No data found for trajectory type: {traj_type}")
            
            # start_states = [0]
            end_states = traj_type_msm_data[self.state_type].unique()
            # print(f'Start states: {start_states}, End states: {end_states}')

            pathways = {}
            states_on_pathway = []
            # start_states = [str(s+1) for s in start_states]
            end_states = [str(s+1) for s in end_states]

            for traj, traj_df in traj_type_msm_data.groupby('traj'):
                # print(f'Analyzing folding pathway of traj: {traj}')
                # print(traj_df)

                path = []
                md = traj_df[self.state_type].values ## get the state labels

                # (1) For each discrete trajectory, put the starting state at the first frame into the pathway;
                    
                # (2) Move forward along the trajectory and find the state that is different with the last state recoded in the pathway. 
                # If the state has not yet been recorded in the pathway, then put it into the pathway. 
                # Otherwise, cut the pathway at the first place where this state is recorded and then move forward;
                start_state = str(md[0]+1) # +1 because the states are 0-indexed in the data
                
                path.append(str(md[0]+1))
                for mdi in md[1:]:
                    tag_find = False
                    for pi in range(len(path)):
                        if path[pi] == str(mdi+1):
                            path = path[0:pi+1]
                            tag_find = True
                            break
                    if not tag_find:
                        path.append(str(mdi+1))

                
                if path[-1] not in end_states:
                    continue      
                
                for p in path:
                    if not int(p) in states_on_pathway:
                        states_on_pathway.append(int(p))
                
                path = ' -> '.join(path)
                # print(f'path: {path}')
                if not path in pathways.keys():
                    pathways[path] = 1
                else:
                    pathways[path] += 1
            
            for k, v in pathways.items():
                self.logger.info(f'{k}: {v}')

            tot_num = 0
            for path in pathways.keys():
                tot_num += pathways[path]
            for path in pathways.keys():
                pathways[path] /= tot_num

            sort_pathways = sorted(pathways.items(), key=lambda x: x[1], reverse=True)
            
            states_on_pathway.sort()
            folding_pathways[traj_type] = {'pathways': sort_pathways, 'states_on_pathway': states_on_pathway}

            for path, prob in sort_pathways:
                folding_pathways_df[traj_type_col].append(traj_type)
                folding_pathways_df['pathway'].append(path)
                folding_pathways_df['probability'].append(prob)

        # print(folding_pathways)
        folding_pathways_df = pd.DataFrame(folding_pathways_df)
        self.logger.debug(folding_pathways_df)
        if outfile_name is None:
            outfile_name = f'FoldingPathways_{self.state_type}_{"-".join(traj_type_list)}.csv'
        outfile = os.path.join(target_outdir, outfile_name)
        folding_pathways_df.to_csv(outfile, index=False)
        self.logger.info(f'Saved folding pathways to {outfile}')

        return folding_pathways_df
    ####################################################################

    ####################################################################
    def JS_divergence(self, msm_data_file:str, traj_type_col:str, traj_type_list:list, meta_set_file:str=None, outdir:str=None, outfile_name:str=None):

        # Load MSM data
        # print(f'Loading MSM data from {self.msm_data}')
        # msm_data = pd.read_csv(self.msm_data)
        # print(f'msm_data:\n{msm_data}')
        # msm_data = msm_data[~msm_data['traj'].isin(self.rm_traj_list)]
        if not msm_data_file:
            raise ValueError('msm_data_file must be provided')
        if not traj_type_list:
            raise ValueError('traj_type_list must contain at least one trajectory type label')

        target_outdir = outdir if outdir is not None else self.outdir
        os.makedirs(target_outdir, exist_ok=True)

        self.logger.info(f'Loading MSM data from {msm_data_file}')
        msm_data = pd.read_csv(msm_data_file)
        msm_data = msm_data[~msm_data['traj'].isin(self.rm_traj_list)]
        self.logger.info(f'msm_data:\n{msm_data}')

        dtrajs = np.asarray([msm_data[msm_data['traj'] == t][self.state_type].values for t in msm_data['traj'].unique()])
        self.logger.info(f'dtrajs shape: {dtrajs.shape}\n{dtrajs}')

        # Load the meta set
        meta_set = []
        if self.state_type == 'microstate':
            if not meta_set_file:
                raise ValueError('meta_set_file must be provided when state_type is microstate')
            meta_set_df = pd.read_csv(meta_set_file)
            self.logger.info(f'meta_set_df:\n{meta_set_df}')
            meta_set = [s['microstates'].values for i, s in meta_set_df.groupby('metastable_state')]
            self.logger.debug(f'meta_set: {meta_set}')

        # Get the number of states present for the state type
        n_states = len(np.unique(dtrajs))
        self.logger.info('Number of %sstates: %d'%(self.state_type, n_states))

        # Get the max number of frames
        max_T_len = dtrajs.shape[1]
        self.logger.debug(f'max_T_len: {max_T_len}')
            
        # make list of traj_idx for each mutant type
        mtype2trajid = {traj_type: [] for traj_type in traj_type_list}
        for i_ax, (traj, traj_df) in enumerate(msm_data.groupby('traj')):
            traj_type = traj_df[traj_type_col].values[0]
            if traj_type in mtype2trajid:
                mtype2trajid[traj_type].append(i_ax)
        self.logger.debug(f'mtype2trajid: {mtype2trajid}')

        # analysis MSM for each mutant
        P_list = []
        for i_ax, traj_type in enumerate(traj_type_list):
            dtrajs_0 = dtrajs[mtype2trajid[traj_type]]
            self.logger.info(f'Processing {traj_type} trajectories with {len(dtrajs_0)} trajectories')
            P_list_0 = np.zeros((max_T_len, n_states))
            for i in range(len(dtrajs_0)):
                # print(f'Processing trajectory {dtrajs_0[i][-self.n_window:]}')
                (N, be) = np.histogram(dtrajs_0[i][-self.n_window:], bins=np.arange(-0.5, n_states, 1))
                dtraj_last = np.argwhere(N == np.max(N))[0][0]
                # print(f'Traj {i} last state histogram: {N} {np.max(N)},\n be: {be},\n dtraj_last: {dtraj_last}')
    
                for j in range(max_T_len):
                    if j >= len(dtrajs_0[i]): 
                        state_0 = dtraj_last
                    else:
                        state_0 = dtrajs_0[i][j]
                    P_list_0[j,state_0] += 1
            P_list.append(P_list_0)
            self.logger.debug(P_list_0)

        # Jensen-Shannon divergence
        JS_list = []
        if self.state_type == 'microstate':
            for ms in meta_set:
                P_list_0 = []
                for i_ax, traj_type in enumerate(traj_type_list):
                    P = np.copy(P_list[i_ax][:,ms])
                    for i in range(len(P)):
                        if np.sum(P[i,:]) != 0:
                            P[i,:] = P[i,:] / np.sum(P[i,:])
                        else:
                            P[i,0] = 1
                    P_list_0.append(P)
                M_0 = 0.5 * (P_list_0[0] + P_list_0[1])
                JS_list.append(0.5 * (entropy(P_list_0[0], M_0, axis=1) + entropy(P_list_0[1], M_0, axis=1)))

        P = P_list[0]
        for i in range(len(P)):
            if np.sum(P[i,:]) != 0:
                P[i,:] = P[i,:] / np.sum(P[i,:])
        self.logger.debug(f'P: {P} {P.shape}')

        Q = P_list[1]
        for i in range(len(Q)):
            if np.sum(Q[i,:]) != 0:
                Q[i,:] = Q[i,:] / np.sum(Q[i,:])
        self.logger.debug(f'Q: {Q} {Q.shape}')

        M = 0.5 * (P + Q)
        self.logger.debug(M)
        entropy_arr = 0.5 * (entropy(P, M, axis=1) + entropy(Q, M, axis=1))
        self.logger.info(f'entropy_arr (n={len(entropy_arr)}): {entropy_arr}')
        # entropy_arr = 0.5 * (entropy(P_list[0], M, axis=1) + entropy(P_list[1], M, axis=1))

        JS_list.append(entropy_arr)
        JS_list = np.array(JS_list).T

        ## write the output file
        if outfile_name is None:
            outfile_name = f'JS_div_{self.state_type}_{"-".join(traj_type_list)}.dat'
        outfile = os.path.join(target_outdir, outfile_name)
        fo = open(outfile, 'w')
        fo.write('%10s '%('Time(s)'))
        for j in range(JS_list.shape[1]-1):
            fo.write('%10s '%('P%d'%(j+1)))
        fo.write('%10s\n'%('JSD'))
        for i in range(max_T_len):
            fo.write('%10.4f '%((i+1)*self.dt))
            for j in range(JS_list.shape[1]):
                fo.write('%10.4f '%(JS_list[i,j]))
            fo.write('\n')
        fo.close()
        self.logger.debug(f'SAVED: {outfile}')
        return JS_list
    ####################################################################

#########################################################################################################################
#########################################################################################################################
class MSMStats:

    """
    A class to analyze various statistical properties of the Markov State Model (MSM) of folding generated by EntDetect.clustering.MSMNonNativeEntanglementClustering
    """

    ####################################################################
    def __init__(self, outdir:str='./MSMStats', 
                 n_window:int=200, 
                 n_traj:int=1000, 
                 state_type:str='metastablestate',
                 rm_traj_list:list=[], log_level:int=logging.INFO, logdir:str=None): 
        """
        Initializes the DataAnalysis class with necessary paths and parameters.

        Parameters:
        - outdir (str): Path to the output directory.
        """
        self.outdir = outdir
        self.n_window = n_window
        self.n_traj = n_traj
        self.state_type = state_type
        self.rm_traj_list = [str(t) for t in rm_traj_list]  # Ensure rm_traj_list is a list of strings
        if self.state_type not in ['microstate', 'metastablestate']:
            raise ValueError(f'state_type must be [microstate | metastablestate]')

        self.logger = setup_logger('MSMStats', outdir=logdir if logdir is not None else self.outdir, log_level=log_level)
        if not os.path.exists(f'{self.outdir}'):
            os.makedirs(f'{self.outdir}')
            self.logger.debug(f'Made output directories {self.outdir}')

        ## set delta time between frames in experimental seconds
        self.dt = 0.015/1000 # time step 
        nsave = 5000 # number of steps between frames
        alpha = 4331293.0
        self.dt = self.dt*nsave*alpha/1e9 # in seconds
        self.end_t = 60 # in seconds
        self.n_boot = 100
        self.num_proc = 20
        self.num_points_plot = 1000
        self.if_boot = True
    ####################################################################

    ####################################################################
    def StateProbabilityStats(self, msm_data_file:str, traj_type_col:str, traj_type_list:list, outfile_name:str='MSTS.csv', outdir:str=None):
        """
        This function calculates the state probabilities from the MSM data and saves them to a file.
        It also calculates the bootstrap statistics for the state probabilities.
        """
        if not traj_type_list:
            raise ValueError('traj_type_list must contain at least one trajectory type label')

        target_outdir = outdir if outdir is not None else self.outdir
        os.makedirs(target_outdir, exist_ok=True)
        outfile = os.path.join(target_outdir, outfile_name)
        if os.path.exists(outfile):
            self.logger.info(f'File {outfile} already exists. Skipping calculation.')
            df = pd.read_csv(outfile)
            self.logger.info(f'Loaded existing data from {outfile}')
            self.logger.info(f'df:\n{df}')
            return df
        
        else:
            # Load MSM data
            self.logger.info(f'Loading MSM data from {msm_data_file}')
            msm_data = pd.read_csv(msm_data_file)
            self.logger.info(f'msm_data:\n{msm_data}')
            msm_data = msm_data[~msm_data['traj'].isin(self.rm_traj_list)]
            self.logger.info(f'msm_data:\n{msm_data}')
            dtrajs = np.asarray([msm_data[msm_data['traj'] == t][self.state_type].values for t in msm_data['traj'].unique()])
            self.logger.info(f'dtrajs shape: {dtrajs.shape}\n{dtrajs}')

            # npzfile = np.load(msm_data_file, allow_pickle=True)

            max_T_len = int(np.ceil(self.end_t/self.dt))
            self.logger.debug(f'max_T_len: {max_T_len}')
            interval = int(max_T_len / self.num_points_plot)
            self.logger.debug(f'interval: {interval}')
            sample_idx = [max_T_len-1-i*interval for i in range(int(max_T_len/interval), -1, -1)]
            if sample_idx[0] != 0:
                sample_idx = [0] + sample_idx
            self.logger.debug(f'sample_idx: {sample_idx}')

            # dtrajs = npzfile['meta_dtrajs']
            n_states = 0
            tag_error = False
            for i, md in enumerate(dtrajs):
                if n_states < np.max(md):
                    n_states = np.max(md)
                if len(md) < max_T_len:
                    #print(f"WARNING: Traj #{i+1} stopped early")
                    tag_error = True
            n_states += 1
            self.logger.debug(f'n_states: {n_states}')

            MSTS_list = []
            boot_stat_list = []
            df = {traj_type_col:[], 'Time(s)':[], 'State':[], 'Probability':[], 'Lower CI':[], 'Upper CI':[]}
            for i_ax, traj_type in enumerate(traj_type_list):
                self.logger.info(f'Processing {traj_type} trajectories')

                ## Get the dtrajs for this trajectory type
                meta_dtrajs = msm_data[msm_data[traj_type_col] == traj_type]
                self.logger.info(f'meta_dtrajs:\n{meta_dtrajs}')
                meta_dtrajs = np.asarray([meta_dtrajs[meta_dtrajs['traj'] == t][self.state_type].values for t in meta_dtrajs['traj'].unique()])
                self.logger.info(f'meta_dtrajs:\n{meta_dtrajs} shape: {meta_dtrajs.shape}')

                
                # MSTS
                PPT = np.zeros((meta_dtrajs.shape[1], n_states))
                self.logger.info(f'PPT shape: {PPT.shape}')
                t_span = (np.arange(meta_dtrajs.shape[1])+1)*self.dt
                self.logger.debug(f't_span: {t_span} shape: {t_span.shape}')

                for md in meta_dtrajs:
                    for i in range(len(t_span)):
                        PPT[i,md[i]]+=1
                PPT /= len(meta_dtrajs)
                self.logger.info(f'PPT:\n{PPT} {PPT.shape}')


                ## Bootstrapping
                self.logger.info(f'Bootstrapping 95% ci...')
                boot_arrs = []
                if self.if_boot:
                    for booti in range(self.n_boot):
                        sample_idx = np.random.choice(np.arange(len(meta_dtrajs)), len(meta_dtrajs), replace=True)
                        boot_meta_dtrajs = meta_dtrajs[sample_idx,:]

                        boot_PPT = np.zeros((boot_meta_dtrajs.shape[1], n_states))
                        for md in boot_meta_dtrajs:
                            for i in range(len(t_span)):
                                boot_PPT[i,md[i]]+=1
                        boot_PPT /= len(boot_meta_dtrajs)
                        # print(f'{booti} boot_PPT:\n{boot_PPT} {boot_PPT.shape}')

                        boot_arrs.append(boot_PPT)

                boot_arrs = np.array(boot_arrs)
                self.logger.info(f'boot_arrs shape: {boot_arrs.shape}')
                lower = np.percentile(boot_arrs, 2.5, axis=0)
                upper = np.percentile(boot_arrs, 97.5, axis=0)
                self.logger.info(f'lower:\n{lower} shape: {lower.shape}')
                self.logger.info(f'upper:\n{upper} shape: {upper.shape}')

                ## make the output dataframe 
                for i in range(PPT.shape[0]):
                    for j in range(n_states):
                        df[traj_type_col].append(traj_type)
                        df['Time(s)'].append(t_span[i])
                        df['State'].append(j+1)  # +1 because states are 0-indexed in the data
                        df['Probability'].append(PPT[i,j])
                        df['Lower CI'].append(lower[i,j])
                        df['Upper CI'].append(upper[i,j])
            df = pd.DataFrame(df)
            df = df.sort_values(by=[traj_type_col, 'State', 'Time(s)'])
            self.logger.debug(df)

            ## save the dataframe to a csv file
            df.to_csv(outfile, index=False)
            self.logger.debug(f'SAVED: {outfile}')
            return df
    ####################################################################

    ####################################################################
    def Plot_StateProbabilityStats(self, df:pd.DataFrame, traj_type_col:str, traj_type_list:list, outdir:str=None, filename_prefix:str=''):
        """
        This function plots the state probabilities from the MSM data.
        Makes one plot for each traj_type.
        It also adds the confidence intervals.
        """
        if df is None:
            raise ValueError('df must be provided. Generate it first with StateProbabilityStats(...).')

        target_outdir = outdir if outdir is not None else self.outdir
        os.makedirs(target_outdir, exist_ok=True)

        for traj_type in traj_type_list:
            prefix = f'{filename_prefix}_' if filename_prefix else ''
            outfile = os.path.join(target_outdir, f'{prefix}{traj_type}_MSTS_plot.png')
            self.logger.info(f'Plotting state probabilities to {outfile}')
            plt.figure(figsize=(10, 6))
            
            traj_df = df[df[traj_type_col] == traj_type]
            for state in traj_df['State'].unique():
                state_df = traj_df[traj_df['State'] == state]
                plt.plot(state_df['Time(s)'], state_df['Probability'], label=f'{traj_type} State {state}')
                plt.fill_between(state_df['Time(s)'], state_df['Lower CI'], state_df['Upper CI'], alpha=0.2)

            plt.xlabel('Time (s)')
            plt.ylabel('Probability')
            plt.title('State Probabilities Over Time')
            plt.legend()
            plt.grid()
            plt.savefig(outfile)
            plt.close()
            self.logger.debug(f'SAVED: {outfile}')
    ####################################################################
#########################################################################################################################
#########################################################################################################################
