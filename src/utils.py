import os
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from conorm import tmm,mrn,getmm
from functools import partial

PROJECT_ROOT = str(Path(__file__).resolve().parents[1]) + "/"

def read_yaml(path: str) -> dict:
    """Returns YAML file as dict"""
    # function written by Olivier Collier olivier.collier@parisnanterre.fr 
    with open(path, 'r') as file_in:
        config = yaml.safe_load(file_in)
    return config

def _rename_columns(name, cohort):
    # function written by Olivier Collier olivier.collier@parisnanterre.fr 
    new_name = name + '@' + cohort
    return new_name

class DataSet:
    # class written by Olivier Collier olivier.collier@parisnanterre.fr 
    """Gathers data from RNA-seq or microarrays.

    Attributes
    ----------
    dataset_name_list: list of str
        A list of subfolder names in the data folder corresponding to the datasets to be used.
    features: pandas.DataFrame
        A dataframe gathering RNA expression data for every patient. Rows correspond to patients and columns to genes or probes.
    metadata: pandas.DataFrame
        A dataframe gathering information about every patient in the dataset. Rows correspond to patients and columns to clinical features.
    common_predictors: list of str
        The list of common predictors (genes or probes) for which expression data is available in every dataset.
    
    Methods
    -------
    load_data()
        Loads data from all datasets in dataset_name_list.
    get_smaller_subset(restrictions)
        Returns smaller subsets with restrictions on patients and predictors.
    """

    def __init__(self, dataset_name_list, dictionary_list=None):
        """
        Parameters:
        -----------
        dataset_name_list: list of str
            A list of subfolder names in the data folder corresponding to the datasets to be used.
        dictionary_list: list of str
            A list of dictionaries to use to convert predictor names for each dataset.
        features: pandas.DataFrame
            A dataframe gathering RNA expression data for every patient.
        metadata: pandas.DataFrame
            A dataframe gathering information about every patient in the dataset. 
        """

        self.dataset_name_list = dataset_name_list
        if dictionary_list:
            self.dictionary_list = dictionary_list
        else:
            self.dictionary_list = [None for _ in dataset_name_list]
        self.features = None
        self.log2_normalized_features = None
        self.anndata = None
        self.group = None
        self.dds = None
        self.stat_res = None
        self.results_df = None
        self.results_df_random = None
        self.results_df_winsorisation_pydeseq2 = None
        self.results_df_Wilcoxon = None
        self.genes_up = None
        self.genes_down = None
        self.mat = None
        self.mat_random = None
        self.mat_winsorisation_pydeseq2 = None
        self.mat_Wilcoxon = None
        self.tf_acts = None
        self.gsea_Hallmark = None
        self.gsea_Kegg = None
        self.gsea_Reactome = None
        self.gsea_myPathways = None

        self.metadata = None
        self.common_predictors = None

    def _find_common_predictors(self): 
        """Computes the list of all predictors that are common to all datasets in dataset_name_list.    

        Returns:
        --------
        common_predictors: list of str
            The list of all predictors that are common to all datasets in dataset_name_list. 
        """

        common_predictors = self._load_predictors(self.dataset_name_list[0])
        common_predictors = self._translate(common_predictors, self.dictionary_list[0])
        for dataset_name, dictionary_name in zip(self.dataset_name_list, self.dictionary_list):
            predictors = self._load_predictors(dataset_name)
            translated_predictors = self._translate(predictors, dictionary_name)
            common_predictors = self._intersection(common_predictors, translated_predictors)
        common_predictors.sort()

        return common_predictors
    
    def _translate(self, predictors, dictionary_name):
        """"Translates predictor names using prescribed dictionary.
        
        Parameters:
        -----------
        predictors: list of str
            A list of predictor names to be translated.
        dictionary_name: str or None
            The name of the file containing a mapping table from old to new predictor names.

        Returns:
        --------
        translated_predictors: list of str
            The list of translated predictor names.
        """

        if dictionary_name:
            dictionary = self._load_dictionary(dictionary_name)
            translated_predictors = [ dictionary.get(p, p) for p in predictors]
            return translated_predictors
        else:
            return predictors

    def _load_dictionary(self, dictionary_name):
        """Loads required dictionary.

        Parameters:
        -----------
        dictionary_name: str
            The name of the file containing the required dictionary in PROJECT_ROOT/src/gene_translation/translators/.

        Returns:
        --------
        dictionary: pd.DataFrame
            A dictionary whose keys are the old predictors, and whose values are the new predictors.
        """

        path_to_dictionary = PROJECT_ROOT + "src/gene_translation/translators/" + dictionary_name
        dictionary_as_dataframe = pd.read_csv(filepath_or_buffer=path_to_dictionary,
                                 sep = ' ',
                                 header = None,
                                 names = ['old_predictors', 'new_predictors'])
        dictionary = {}
        for index, row in dictionary_as_dataframe.iterrows():
            old_predictor = row['old_predictors']
            new_predictor = row['new_predictors']
            dictionary[old_predictor] = new_predictor

        return dictionary

    def load_data(self):
        """Loads data from all datasets in dataset_name_list.
        """

        self.common_predictors = self._find_common_predictors()
        self.features = self._load_features()
        self.metadata = self._load_metadata()

    def _load_predictors(self, dataset_name):
        """Loads list of predictors corresponding to a particular dataset.

        Parameters:
        -----------
        dataset_name: str
            The name of the subfolder in the data folder corresponding to the considered dataset.

        Returns:
        --------
        predictors: list of str
            The list of predictors corresponding to the specified dataset.
        """

        if dataset_name.endswith('/'):
            dataset_name = dataset_name[:-1]
        if self.dataset_type == 'micro_array':
            path_to_data = PROJECT_ROOT + 'data/' + dataset_name + '/treated_files/probes.txt' 
        elif self.dataset_type.startswith('RNA-seq'):
            path_to_data = PROJECT_ROOT + 'data/' + dataset_name + '/treated_files/genes.txt'
        predictors = pd.read_csv(path_to_data, sep = '\t', header=None)[0].tolist()
        predictors.sort()

        return predictors

    def _load_features(self):
        """Loads features corresponding to all considered datasets.
        
        Returns:
        --------
        total_features: pandas.DataFrame
            A dataframe containing features corresponding to all considered datasets.
        """
        
        total_features = pd.DataFrame([])
        for dataset_name, dictionary_name in zip(self.dataset_name_list, self.dictionary_list):
            if dataset_name.endswith('/'):
                dataset_name = dataset_name[:-1]
            features = self._load_features_from_dataset(dataset_name, dictionary_name)
            if total_features.shape == (0, 0):
                total_features = features
            else:
                total_features = pd.concat([total_features, features], axis=0)
        total_features.sort_index(axis=0, inplace=True)
        total_features.sort_index(axis=1, inplace=True)

        return total_features
    
    def _load_features_from_dataset(self, dataset_name, dictionary_name):
        """Loads features corresponding to a particular dataset.

        Returns:
        --------
        features: pandas.DataFrame
            A dataframe containing features corresponding to a particular dataset.
        """
        
        if dataset_name.endswith('/'):
            dataset_name = dataset_name[:-1]
        path_to_data = PROJECT_ROOT + 'data/' + dataset_name + '/treated_files/counts.txt'
        features = pd.read_csv(path_to_data, sep = '\t')
        if dictionary_name:
            dictionary = self._load_dictionary(dictionary_name)
            features = features.rename(index=dictionary)
            features = features.loc[~features.index.duplicated(keep='first'), ]
        features = features.transpose()
        features = features[self.common_predictors]
        features.sort_index(axis=0, inplace=True)
        features.sort_index(axis=1, inplace=True)

        return features
    
    def _load_metadata(self):
        """Loads additional information on every patient in the considered datasets.

        Returns:
        --------
        total_metadata: pandas.DataFrame
            A dataframe containing additional information on every patient in the considered datasets.
        """

        total_metadata = pd.DataFrame([])
        for dataset_name in self.dataset_name_list:
            metadata = self._load_metadata_from_dataset(dataset_name)
            if total_metadata.shape == (0, 0):
                total_metadata = metadata
            else:
                total_metadata = pd.concat([total_metadata, metadata], axis=0)
        total_metadata.sort_index(inplace=True)
        total_metadata = total_metadata.fillna('NA')

        return total_metadata

    def _load_metadata_from_dataset(self, dataset_name):
        """Loads additional information on every patient in a particular dataset.

        The names of the features in the dataset will be changed according to the config file included in the dataset_name subfolder.

        Parameters:
        -----------
        dataset_name: str
            The name of the subfolder in the data folder containing the dataset of interest.

        Returns:
        --------
        total_metadata: pandas.DataFrame
            A dataframe containing additional information on every patient in a particular dataset.
        """

        path_to_data = PROJECT_ROOT + 'data/' + dataset_name + '/treated_files/metadata.txt'
        metadata = pd.read_csv(path_to_data, sep = '\t')
        metadata = metadata.transpose()
        config = self._load_config_from_dataset(dataset_name)
        metadata = self._rename_metadata(metadata, config)
        metadata.sort_index(inplace=True)

        return metadata

    def _load_config_from_dataset(self, dataset_name):
        """Loads config file giving information on a particular dataset.

        Parameters:
        -----------
        dataset_name: str
            The name of the subfolder in the data folder containing the dataset of interest.

        Returns:
        --------
        total_metadata: pandas.DataFrame
            A dataframe containing additional information on every patient in a particular dataset.
        """
        
        path_to_data = PROJECT_ROOT + 'data/' + dataset_name + '/treated_files/config.yml'
        with open(path_to_data, 'r') as file_in:
            config = yaml.safe_load(file_in)
        return config

    def _rename_metadata(self, metadata, config):

        for predictor_name in list(config.keys()):
            old_predictor_name = config[predictor_name]['name']
            
            # rename columns
            if old_predictor_name in metadata.columns:
                metadata = metadata.rename(columns={old_predictor_name: predictor_name})
            else:
                # skip silently if column not in metadata
                continue

            # remap categorical values only if found in config
            if config[predictor_name]['type'] == 'categorical':
                dictionary = config[predictor_name]['values']

                # safely map values without KeyError
                metadata[predictor_name] = metadata[predictor_name].map(dictionary).fillna(metadata[predictor_name])

        return metadata


    def _intersection(self, list1, list2):
        """Computes the intersection between two lists.
        
        Parameters:
        -----------
        list1, list2: list
        """

        return [elt for elt in list1 if elt in list2]

    def get_smaller_subset(self, restrictions):
        """Returns a smaller subset with restrictions on patients and predictors.

        Parameters:
        -----------
        restrictions: dict
            A dictionary with two keys: 'probes'/'genes' (depending on the type of data) and 'patients'. 
            The value for 'probes'/'genes' is a list of probes/genes.
            The value for 'patients' is a dictionary whose keys are feature names corresponding to the metadata, and whose values are the values of these features to select.
        
        Returns:
        --------
        new_dataset: pandas.DataFrame
            Object DataSet with restrictions on patients and predictors.
        """

        if self.dataset_type == 'micro_array':
            new_dataset = MicroArray_DataSet(self.dataset_name_list)
            predictor_restrictions = restrictions['probes']
        elif self.dataset_type.startswith('RNA-seq'):
            new_dataset = RNAseq_DataSet(self.dataset_name_list)
            predictor_restrictions = restrictions['genes']
        patient_restrictions = restrictions['patients']
        
        new_dataset.common_predictors = self.common_predictors
        relevant_predictors = self._get_relevant_predictors(predictor_restrictions) 
        relevant_patients = self._get_relevant_patients(patient_restrictions)        
        new_dataset.common_predictors = relevant_predictors
        new_dataset.features = self._apply_restrictions_for_features(relevant_patients, relevant_predictors)
        new_dataset.metadata = self._apply_restrictions_for_metadata(relevant_patients)

        return new_dataset
    
    def _apply_restrictions_for_features(self, relevant_patients, relevant_predictors):
        """Transforms feature dataframe to take restrictions into account.

        Parameters:
        -----------
        relevant_patients: list
            List of patients to select.
        relevant_predictors: list
            List of predictors to select.

        Returns:
        --------
        new_features: pandas.DataFrame
            A dataframe of features restricted to prescribed patients and predictors.
        """
        
        new_features = self.features.loc[relevant_patients]
        new_features = new_features[relevant_predictors]
        new_features.sort_index(inplace=True)

        return new_features

    def _apply_restrictions_for_metadata(self, relevant_patients):
        """Transforms metadata dataframe to take restrictions into account.

        Parameters:
        -----------
        relevant_patients: list
            List of patients to select.
        
        Returns:
        --------
        new_features: pandas.DataFrame
            A dataframe of metadata restricted to prescribed patients.
        """

        new_metadata = self.metadata.loc[relevant_patients]
        new_metadata.sort_index(inplace=True)

        return new_metadata

    def _get_relevant_patients(self, patient_restrictions):
        """Returns the list of patients corresponding to given restrictions.
        
        Parameters:
        -----------
        patient_restrictions: dict
            A dictionary whose keys are feature names corresponding to the metadata, and whose values are the values of these features to select.
            For example, {'FEATURE': FEATURE_VALUE} corresponds to selecting patients for which FEATURE takes the value FEATURE_VALUE.

        Returns:
        --------
        relevant_patients: list
            The list of patients corresponding to prescribed restrictions on patients.
        """

        relevant_patients = self.metadata.index
        if patient_restrictions != None:
            for predictor_name in patient_restrictions:
                authorized_values = patient_restrictions[predictor_name]
                relevant_patients = [patient for patient in relevant_patients if self.metadata[predictor_name][patient] in authorized_values]

        return relevant_patients
    
    def _get_relevant_predictors(self, predictor_restrictions):
        """Returns the list of predictors corresponding to given restrictions.
        
        Parameters:
        -----------
        predictor_restrictions: list
            A list of predictors (probes or genes) to select, if they exist in the dataset.

        Returns:
        --------
        relevant_predictors: list
            The list of predictors to select.
        """

        relevant_predictors = self.common_predictors
        if predictor_restrictions != None:
            relevant_predictors = self._intersection(relevant_predictors, predictor_restrictions) 

        return relevant_predictors
    
    def _get_group_indices(self):
        """Returns a list with the number of the corresponding group for each predictor.

        A group is made of the indices corresponding to the same biological predictor accross different features. 
        
        Returns:
        --------
        group_indices: list

        """
    
        number_predictors = len(self.common_predictors)
        group_indices = [i for i in range(1, number_predictors+1) for _ in self.dataset_name_list]

        return group_indices

class MicroArray_DataSet(DataSet):
    # class written by Olivier Collier olivier.collier@parisnanterre.fr 
    """Gathers data from microarray.

    Child of the DataSet class.
    """

    def __init__(self, dataset_name_list, dictionary_list=None):
        self.dataset_type = 'micro_array'
        DataSet.__init__(self, dataset_name_list, dictionary_list)

class RNAseq_DataSet(DataSet):
    # Adapeted from class written by Olivier Collier olivier.collier@parisnanterre.fr 
    """Gathers data from RNA-seq.
    
    Child of the DataSet class.
    """

    def __init__(self, dataset_name_list, dictionary_list=None):
        self.dataset_type = 'RNA-seq:counts'
        DataSet.__init__(self, dataset_name_list, dictionary_list)


    def transform_counts_to_TMM_normalized_counts(self):
        """Transforms counts to TMM-normalized counts.
        """

        if self.dataset_type == 'RNA-seq:TMM-normalized_counts':
            pass
        elif self.dataset_type == 'RNA-seq:counts':
            tmm_normalized_counts = tmm(self.features.transpose()).transpose()
            self.features = tmm_normalized_counts
            self.dataset_type = 'RNA-seq:TMM-normalized_counts'
        elif self.dataset_type == 'RNA-seq:TMM-normalized_log2counts':
            tmm_normalized_log2counts = np.power(2, self.features) - 1
            self.features = tmm_normalized_log2counts
            self.dataset_type = 'RNA-seq:TMM-normalized_counts'

    def transform_counts_to_TMM_normalized_log2counts(self):
        """Transforms counts to TMM-normalized and log2-transformed cuonts.
        """

        if self.dataset_type == 'RNA-seq:TMM-normalized_log2counts':
            pass
        elif self.dataset_type == 'RNA-seq:counts':
            tmm_normalized_counts = tmm(self.features.transpose()).transpose()
            tmm_normalized_log2counts = np.log2(tmm_normalized_counts + 1)
            self.features = tmm_normalized_log2counts
            self.dataset_type = 'RNA-seq:TMM-normalized_log2counts'
        elif self.dataset_type == 'RNA-seq:TMM-normalized_counts':
            tmm_normalized_log2counts = np.log2(self.features + 1)
            self.features = tmm_normalized_log2counts
            self.dataset_type = 'RNA-seq:TMM-normalized_log2counts'

    def transform_counts_to_GeTMM_normalized_log2counts(self):
        """Transforms counts to GeTMM-normalized and log2-transformed cuonts.
        """

        if self.dataset_type == 'RNA-seq:GeTMM-normalized_log2counts':
            pass
        elif self.dataset_type == 'RNA-seq:counts':
            # gene_length
            gene_lengths = pd.read_csv(PROJECT_ROOT + 'data/gene_lenght/Homo_sapiens.GRCh38.112_gene_lenght.csv')
            gene_lengths.set_index('gene_name', inplace=True)
            expression_data = self.features.transpose()
            common_genes = common_genes = expression_data.index.intersection(gene_lengths.index)
            expression_data = expression_data.loc[common_genes]
            gene_lengths= gene_lengths.loc[common_genes]
            # normalization
            getmm_normalized_counts = getmm(expression_data, gene_lengths).transpose()
            # Update genes.txt if needed.
            # log2 transformation
            getmm_normalized_log2counts = np.log2(getmm_normalized_counts + 1)
            self.features = getmm_normalized_log2counts
            self.dataset_type = 'RNA-seq:GeTMM-normalized_log2counts'

        elif self.dataset_type == 'RNA-seq:GeTMM-normalized_counts':
            getmm_normalized_log2counts = np.log2(self.features + 1)
            self.features = getmm_normalized_log2counts
            self.dataset_type = 'RNA-seq:GeTMM-normalized_log2counts'

def pool_dataset(dict_dataset, pool_dataset_name):
    """
    Pool several datasets while retaining all metadata fields.
    - Concatenate count matrices by rows.
    - Concatenate metadata tables while automatically aligning columns.
    - Add a 'dataset' column to record the source dataset for each sample.
    """

    dst = RNAseq_DataSet([pool_dataset_name]) 
    Raw_counts = None
    metadata_ = pd.DataFrame()

    for dataset_name, dataset in dict_dataset.items():

        # --- 1) Counts ---
        counts = dataset.features.copy()
        counts = counts.round().astype(int)

        # Keep only genes shared across datasets.
        if Raw_counts is None:
            Raw_counts = counts
        else:
            common_genes = Raw_counts.columns.intersection(counts.columns)
            Raw_counts = Raw_counts[common_genes]
            counts = counts[common_genes]
            Raw_counts = pd.concat([Raw_counts, counts], axis=0)

        # --- 2) Metadata ---
        meta = dataset.metadata.copy()

        # Add source dataset information.
        meta["dataset"] = dataset_name

        # Concatenate while automatically aligning columns.
        metadata_ = pd.concat([metadata_, meta], axis=0, sort=False)

    # finalisation
    dst.features = Raw_counts
    dst.metadata = metadata_

    dst.common_predictors = list(dst.features.columns)
    dst.dataset_type = "RNA-seq:counts"

    return dst

  
def load_tcga_primary(project_id, patient_filter=None):
    """
    project_id: str, ex 'TCGA-BRCA'
    patient_filter: optional dictionary passed to get_smaller_subset.
    """
    ds = RNAseq_DataSet([project_id])
    ds.load_data()
    base_filter = {'genes': None, 'patients': {'sample_type': ['Primary Tumor']}}

    if patient_filter is not None:
        # Merge filtering criteria.
        base_filter['patients'].update(patient_filter)

    ds = ds.get_smaller_subset(base_filter)
    ds.group = ds.metadata
    return ds

def label_by_terciles(values: pd.Series,
                      lower_q=1/3,
                      upper_q=2/3,
                      low_label="LOW",
                      high_label="HIGH",
                      unknown_label="unknown") -> pd.Series:
    """
    Assign low_label to values <= q_low,
    high_label to values >= q_high,
    and unknown_label to the remaining values.
    """
    q_low = values.quantile(lower_q)
    q_high = values.quantile(upper_q)

    lab = pd.Series(unknown_label, index=values.index, dtype=object)
    lab[values <= q_low] = low_label
    lab[values >= q_high] = high_label

    return lab

def label_terciles(vst_df, merged_meta, condition_a_tester, name_gene):
    expr = vst_df.loc[mask_test, name_gene]
    low_lab, high_lab = condition_a_tester[1], condition_a_tester[2]
    return label_by_terciles(expr, 1/3, 2/3, low_lab, high_lab, 'unknown').rename('condition_transfer')
#if __name__ == "__main__":

