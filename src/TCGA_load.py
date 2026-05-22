import requests
import pandas as pd
import os
import numpy as np
from lxml import etree
import argparse
import shutil
import requests
import time
import urllib3
from pathlib import Path

# GDC API URL
GDC_API_URL = "https://api.gdc.cancer.gov/"
PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _path_with_trailing_slash(path):
    path = str(path)
    return path if path.endswith("/") else path + "/"


def search_files(endpoint, filters, fields, file_size=None):
    """ 
    Search for files in the GDC Data Portal using the GDC API.

    :param endpoint: The endpoint to search (e.g., "files", "cases", "projects")
    :type endpoint: str
    :param filters: A dictionary containing the filters to apply to the search
    :type filters: dict
    :param fields: A string containing the fields to return in the search results
    :type fields: str
    :param file_size: The maximum file size to return (optional)
    :type file_size: str
    :return: A list of files matching the search criteria
    :rtype: list
    :note: The filters and fields should be formatted according to the GDC API documentation
    """
    params = {
        "filters": filters,
        "fields": fields,
        "format": "JSON",
        "size": file_size if file_size else "10000"
    }
    
    response = requests.post(GDC_API_URL + endpoint, json=params)
    response.raise_for_status()
    return response.json()["data"]["hits"]



def download_file(file_id, path_name_file, max_retries=5):
    """
    Download a file from the GDC Data Portal using the GDC API.
    
    :param file_id: The GDC UUID of the file to download
    :type file_id: str
    :param path_name_file: The path and name of the file to save the data to
    :type path_name_file: str
    :param max_retries: The maximum number of times to retry downloading the file in case of connection errors (default: 5)
    :type max_retries: int, optional
    :return: None
    :rtype: None
    :note: If the request fails due to network issues, it will retry up to max_retries times before raising an exception.
    """

        
    GDC_API_URL = "https://api.gdc.cancer.gov/"
    data_endpoint = f"{GDC_API_URL}data/{file_id}"
    
    for attempt in range(max_retries):
        try:
            with requests.get(data_endpoint, stream=True, timeout=(10, 30)) as r:
                r.raise_for_status()
                with open(path_name_file, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            #print(f"Download successful: {path_name_file}")
            return
        except requests.exceptions.ConnectionError as e:
            #print(f"ConnectionError: {e}, retrying {attempt+1}/{max_retries}")
            time.sleep(5)  
        except requests.exceptions.HTTPError as e:
            #print(f"HTTPError: {e}, stopping retries.")
            break
        except requests.exceptions.Timeout:
            #print(f"Timeout error, retrying {attempt+1}/{max_retries}")
            time.sleep(5)
    
    print("Failed to download file after multiple attempts.")


def download_file_urllib3(file_id, path_name_file, max_retries=5):
    """
    Download a file from the GDC Data Portal using the GDC API.
    
    :param file_id: The GDC UUID of the file to download
    :type file_id: str
    :param path_name_file: The path and name of the file to save the data to
    :type path_name_file: str
    :param max_retries: The maximum number of times to retry downloading the file in case of connection errors (default: 5)
    :type max_retries: int, optional
    :return: None
    :rtype: None
    :note: If the request fails due to network issues, it will retry up to max_retries times before raising an exception.
    """
    http = urllib3.PoolManager()
    data_endpoint = f"https://api.gdc.cancer.gov/data/{file_id}"

    try:
        response = http.request("GET", data_endpoint, preload_content=False)
        with open(path_name_file, 'wb') as f:
            f.write(response.data)
        response.release_conn()
        #print(f"Download successful: {path_name_file}")
    except urllib3.exceptions.HTTPError as e:
        print(f"HTTP Error: {e}")

# DATA RNA-SEQ

def get_rnaseq_data(project_name):
    """
    Get RNA-seq data for a specific TCGA project from the GDC Data Portal.
    
    :param project_name: The name of the TCGA project to search for (e.g., "TCGA-BRCA")
    :type project_name: str
    :return: A list of RNA-seq files for the specified project
    :rtype: list
    """
    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": project_name}},
            {"op": "=", "content": {"field": "files.data_category", "value": "Transcriptome Profiling"}},
            {"op": "=", "content": {"field": "files.data_type", "value": "Gene Expression Quantification"}},
            {"op": "=", "content": {"field": "files.analysis.workflow_type", "value": "STAR - Counts"}}
        ]
    }

    fields = "file_id,file_name,cases.submitter_id,cases.samples.sample_id,cases.samples.sample_type"
    
    return search_files("files", filters, fields)

def get_counts_data( rna_seq_files, path_to_download = '../data/download_data/'):
    """
    Get RNA-seq counts data for a specific TCGA project from the GDC Data Portal.

    :param rna_seq_files: A list of RNA-seq files for the specified project
    :type rna_seq_files: list

    :param path_to_download: The path to save the downloaded files
    :type path_to_download: str
    :return: a tuple containing the combined counts data and the sample metadata in  DataFrames
    :rtype: tuple of pd.DataFrame
    """
    # Create a directory to store the downloaded files
    if not os.path.exists(path_to_download+'RNA_seq'):
        os.makedirs(path_to_download+'RNA_seq')
    # Loop through the RNA-seq files and download each one
    for file in rna_seq_files:
        if file['file_name'] not in os.listdir(path_to_download+'RNA_seq'):
            download_file_urllib3(file['file_id'], path_to_download+ f'RNA_seq/{file["file_name"]}')
    downloaded_files_RNA = os.listdir(path_to_download+'RNA_seq')
    dataframes =[]
    downloaded_files_RNA_2 = []
    for file in downloaded_files_RNA:
        try:
            dataframes.append(pd.read_csv(path_to_download+ f'RNA_seq/{file}', sep='\t', header=1))
            downloaded_files_RNA_2.append(file)
        except:
            print(f'error for file {file}')

    #dataframes = [pd.read_csv(path_to_download+ f'RNA_seq/{file}', sep='\t', header=1) for file in downloaded_files_RNA ]

    # Drop rows with missing values
    dataframes = [df.dropna(axis=0) for df in dataframes] 
    # Merge the dataframes to combine the counts
    combined_counts = pd.concat([df['unstranded'] for df in dataframes], axis=1)

    # Adjust the column names
    combined_counts.index = dataframes[0]['gene_name']
    combined_counts.index.name = None
    combined_counts.columns = downloaded_files_RNA_2

    # Extract the patient and sample IDs from the RNA-seq files
    data = []
    for file in rna_seq_files:
        if file['file_name'] in downloaded_files_RNA_2:
            file_id = file['file_id']
            file_name = file['file_name']
            submitter_id = file['cases'][0]['submitter_id']  # patient ID
            sample_id = file['cases'][0]['samples'][0]['sample_id']  # sample ID
            sample_type = file['cases'][0]['samples'][0]['sample_type']  # sample type (e.g., Primary Tumor)
            
            # Add the data to the list
            data.append({
                "file_id": file_id,
                "file_name": file_name,
                "patient_id": submitter_id,
                "sample_id": sample_id,
                "sample_type": sample_type
            })

    # Convert the list to a DataFrame
    df_ = pd.DataFrame(data)
    df_["patient_id_2"] = df_["patient_id"].str.split("-").str[2]

    # ----------------------------------------------------------
    #  BATCH EFFECT EXTRACTION (TSS, CENTER, PLATE, ALIQUOT)
    # ----------------------------------------------------------

    # We request technical metadata from the GDC "cases" endpoint
    batch_filters = {
        "op": "in",
        "content": {
            "field": "samples.sample_id",
            "value": df_["sample_id"].tolist()
        }
    }

    batch_fields = ",".join([
        "submitter_id",
        "samples.sample_id",
        "samples.tissue_source_site",                                 # TSS batch
        "samples.portions.analytes.aliquots.submitter_id",            # aliquot barcode
        "samples.portions.analytes.aliquots.center",                  # sequencing center
        "samples.portions.plate"                                      # plate ID
    ])

    batch_info = search_files("cases", batch_filters, batch_fields)

    # Parse the batch metadata returned by GDC
    batch_records = []
    for case in batch_info:
        case_submitter = case.get("submitter_id", "unknown")

        for sample in case.get("samples", []):
            sample_id = sample.get("sample_id", "unknown")
            tss = sample.get("tissue_source_site", "unknown")

            # Default values
            plate = "unknown"
            center = "unknown"
            aliquot = "unknown"

            if "portions" in sample and len(sample["portions"]) > 0:
                portion = sample["portions"][0]

                # Plate
                if "plate" in portion:
                    plate = portion["plate"]

                # Aliquot + center
                if "analytes" in portion and len(portion["analytes"]) > 0:
                    analyte = portion["analytes"][0]
                    if "aliquots" in analyte and len(analyte["aliquots"]) > 0:
                        aliquot_data = analyte["aliquots"][0]
                        aliquot = aliquot_data.get("submitter_id", "unknown")
                        center = aliquot_data.get("center", "unknown")

            batch_records.append({
                "sample_id": sample_id,
                "TSS": tss,
                "plate": plate,
                "center": center,
                "aliquot": aliquot
            })

    df_batch = pd.DataFrame(batch_records)

    # Merge batch info into df_
    df_ = df_.merge(df_batch, on="sample_id", how="left")

    # Replace missing values
    df_[["TSS", "plate", "center", "aliquot"]] = df_[["TSS", "plate", "center", "aliquot"]].fillna("unknown")


    # Set the sample ID as the index
    combined_counts_2 = combined_counts[df_['file_name']]
    combined_counts_2.columns = df_['sample_id']
    combined_counts_2.columns.name = None
    

    return combined_counts_2, df_

# CLINICAL DATA
def download_clinical_data(project_name, path_to_download = '../data/download_data/'):
    """
    Download clinical data for a specific TCGA project from the GDC Data Portal.

    :param project_name: The name of the TCGA project to search for (e.g., "TCGA-BRCA")
    :type project_name: str
    :param path_to_download: The path to save the downloaded files
    :type path_to_download: str
    :return: None
    :rtype: None
    """

    clinical_filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": project_name}},
            {"op": "=", "content": {"field": "files.data_category", "value": "Clinical"}},
            {"op": "=", "content": {"field": "files.data_format", "value": "BCR XML"}}
        ]
    }

    survival_fields = "file_id,file_name,cases.submitter_id,cases.diagnoses.days_to_death,cases.diagnoses.vital_status,cases.diagnoses.days_to_last_follow_up"
    survival_files = search_files("files", clinical_filters, survival_fields)
    print(f"Number of clinical files to download: {len(survival_files)}")

    if not os.path.exists(path_to_download+'clinical'):
        os.makedirs(path_to_download+'clinical')

    for file in survival_files:
        if file['file_name'] not in os.listdir(path_to_download+'clinical'):
            download_file(file['file_id'], path_to_download+ f'clinical/{file["file_name"]}' )

        
def load_clinical_data(project_name, df_, path_to_load, path_to_download= '../data/download_data/'):
    project_name_2 = project_name.split('-')[1].lower()

    clinical_data = {
        'patient_id': [],
        'gender': [],
        'vital_status': [],
        'tumor_tissue_site': [],
        'race': [],
        'ethnicity': [],
        'age_at_diagnosis': [],
        'days_to_death': [],                # for Kaplan-Meier analysis
        'days_to_last_follow_up': [],       # for living patients
        # --- batch effect ---
        'batch_number': [],
        'center': [],
        'file_uuid': [],
        'project_code': [],
        'disease_code': [],
    }

    # Uncomment these keys if additional clinical fields are needed.
    # clinical_data['clinical_stage'] = []
    # clinical_data['pathologic_stage'] = []
    # clinical_data['tumor_grade'] = []
    # clinical_data['tumor_status'] = []
    # clinical_data['residual_tumor'] = []
    # clinical_data['primary_therapy_outcome'] = []

    if project_name == 'TCGA-BRCA':
        clinical_data['er_status_by_ihc'] = []
        clinical_data['pr_status_by_ihc'] = []
        clinical_data['her2_status_by_ihc'] = []

    downloaded_files_clinical = [f for f in os.listdir(path_to_download+'clinical') if f.endswith('.xml')]

    for xml_file in downloaded_files_clinical:
        try:
            tree = etree.parse(path_to_download + '/clinical/{}'.format(xml_file))
            root = tree.getroot()
            nsmap = root.nsmap

            # -------------------------
            # admin:admin (batch information) - once per file
            # -------------------------
            admin_block = root.find('.//admin:admin', namespaces=nsmap)
            if admin_block is not None:
                batch_number = admin_block.findtext('admin:batch_number', namespaces=nsmap)
                center = admin_block.findtext('admin:bcr', namespaces=nsmap)
                file_uuid = admin_block.findtext('admin:file_uuid', namespaces=nsmap)
                project_code = admin_block.findtext('admin:project_code', namespaces=nsmap)
                disease_code = admin_block.findtext('admin:disease_code', namespaces=nsmap)

                # Use "unknown" when the field is missing.
                batch_number = batch_number or "unknown"
                center = center or "unknown"
                file_uuid = file_uuid or "unknown"
                project_code = project_code or "unknown"
                disease_code = disease_code or "unknown"
            else:
                batch_number = center = file_uuid = project_code = disease_code = "unknown"

            # -------------------------
            # patients
            # -------------------------
            for patient in root.findall(f'.//{project_name_2}:patient', namespaces=nsmap):

                patient_id = patient.findtext('.//shared:bcr_patient_barcode', namespaces=nsmap)
                gender = patient.findtext('.//shared:gender', namespaces=nsmap)
                vital_status = patient.findtext('.//clin_shared:vital_status', namespaces=nsmap)
                tumor_tissue_site = patient.findtext('.//clin_shared:tumor_tissue_site', namespaces=nsmap)
                race = patient.findtext('.//clin_shared:race', namespaces=nsmap)
                ethnicity = patient.findtext('.//clin_shared:ethnicity', namespaces=nsmap)
                age_at_diagnosis = patient.findtext('.//clin_shared:age_at_initial_pathologic_diagnosis', namespaces=nsmap)

                days_to_death = patient.findtext('.//clin_shared:days_to_death', namespaces=nsmap)
                days_to_last_follow_up = patient.findtext('.//clin_shared:days_to_last_follow_up', namespaces=nsmap)

                # Batch information is identical for all patients in this file.
                clinical_data['batch_number'].append(batch_number)
                clinical_data['center'].append(center)
                clinical_data['file_uuid'].append(file_uuid)
                clinical_data['project_code'].append(project_code)
                clinical_data['disease_code'].append(disease_code)

                if project_name == 'TCGA-BRCA':
                    er_status_by_ihc = patient.findtext('.//brca_shared:breast_carcinoma_estrogen_receptor_status', namespaces=nsmap)
                    pr_status_by_ihc = patient.findtext('.//brca_shared:breast_carcinoma_progesterone_receptor_status', namespaces=nsmap)
                    her2_status_by_ihc = patient.findtext('.//brca_shared:lab_proc_her2_neu_immunohistochemistry_receptor_status', namespaces=nsmap)

                clinical_data['patient_id'].append(patient_id)
                clinical_data['gender'].append(gender)
                clinical_data['vital_status'].append(vital_status)
                clinical_data['tumor_tissue_site'].append(tumor_tissue_site)
                clinical_data['race'].append(race)
                clinical_data['ethnicity'].append(ethnicity)
                clinical_data['age_at_diagnosis'].append(age_at_diagnosis)
                clinical_data['days_to_death'].append(days_to_death)
                clinical_data['days_to_last_follow_up'].append(days_to_last_follow_up)

                # If these fields are needed, also uncomment their keys in the dictionary above.
                # try:
                #     clinical_stage = patient.findtext('.//shared_stage:clinical_stage', namespaces=nsmap)
                #     pathologic_stage = patient.findtext('.//shared_stage:pathologic_stage', namespaces=nsmap)
                #     tumor_grade = patient.findtext('.//shared:neoplasm_histologic_grade', namespaces=nsmap)
                #     tumor_status = patient.findtext('.//clin_shared:person_neoplasm_cancer_status', namespaces=nsmap)
                #     residual_tumor = patient.findtext('.//clin_shared:residual_tumor', namespaces=nsmap)
                #     primary_therapy_outcome = patient.findtext('.//clin_shared:primary_therapy_outcome_success', namespaces=nsmap)
                #     clinical_data['clinical_stage'].append(clinical_stage)
                #     clinical_data['pathologic_stage'].append(pathologic_stage)
                #     clinical_data['tumor_grade'].append(tumor_grade)
                #     clinical_data['tumor_status'].append(tumor_status)
                #     clinical_data['residual_tumor'].append(residual_tumor)
                #     clinical_data['primary_therapy_outcome'].append(primary_therapy_outcome)
                # except Exception:
                #     print(f"Not all clinical data in file {xml_file}")

                if project_name == 'TCGA-BRCA':
                    clinical_data['er_status_by_ihc'].append(er_status_by_ihc)
                    clinical_data['pr_status_by_ihc'].append(pr_status_by_ihc)
                    clinical_data['her2_status_by_ihc'].append(her2_status_by_ihc)

        except etree.XMLSyntaxError as e:
            print(f"Syntax error in file {xml_file}: {e}")
        except Exception as e:
            print(f"Error in file {xml_file}: {e}")

    if len(clinical_data['patient_id']) == 0:
        raise ValueError("No clinical data extracted. See the XML files.")

    # --- Continue with the unchanged part of the function. ---
    df_clinical_data = pd.DataFrame(clinical_data)
    df_clinical_data['days_to_death'] = pd.to_numeric(df_clinical_data['days_to_death'], errors='coerce')
    df_clinical_data['days_to_last_follow_up'] = pd.to_numeric(df_clinical_data['days_to_last_follow_up'], errors='coerce')
    df_clinical_data['vital_status'] = df_clinical_data['vital_status'].apply(lambda x: 1 if x == "Dead" else 0)
    df_clinical_data['survival_time'] = df_clinical_data[['days_to_death', 'days_to_last_follow_up']].max(axis=1)

    df_clinical_data_ = df_[['sample_id','patient_id','sample_type']].merge(df_clinical_data, on='patient_id', how='left')
    df_clinical_data_ = df_clinical_data_.replace('', 'unknown').replace(np.nan, 'unknown')

    df_clinical_data_.index = df_clinical_data_['sample_id']
    df_clinical_data_.index.name = None

    df_clinical_data_['COHORT'] = project_name
    df_clinical_data_ = df_clinical_data_.transpose()
    print(df_clinical_data_.T['sample_type'].value_counts())

    if not os.path.exists(path_to_load):
        os.makedirs(path_to_load)

    df_clinical_data_.columns = df_clinical_data_.columns.str.lstrip()
    df_clinical_data_.to_csv(path_to_load+'metadata.txt', sep='\t')
    remove_first_character(path_to_load+'metadata.txt')

    return df_clinical_data_


def remove_first_character(file_path):
    try:
        # Open the file and read its content
        with open(file_path, 'r', encoding='utf-8') as file:
            content = file.read()
        
        if content:
            # Remove the first character
            new_content = content[1:]
            
            # Write the updated content back to the file
            with open(file_path, 'w', encoding='utf-8') as file:
                file.write(new_content)
            
            print("First character successfully removed.")
        else:
            print("The file is empty.")
    except FileNotFoundError:
        print("File not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

def load_counts_data(combined_counts_2, df_clinical_data_, path_to_load):
    """
    Load 

    """

    # keep samples with clinical data
    common_columns = list(set(df_clinical_data_.columns).intersection(set(combined_counts_2.columns)))
    df_counts_ = combined_counts_2[common_columns]
    # Remove genes with very low read counts.
    df_counts_ = df_counts_.loc[df_counts_.sum(axis=1)> 10]
    # Replace missing values with 0.
    df_counts_ = df_counts_.fillna(0)
    # Average genes with identical names.
    df_counts_ = df_counts_.groupby(df_counts_.index).mean()

    df_counts_.to_csv(path_to_load+'counts.txt', sep = '\t')
    # Drop the first caracter of the file
    remove_first_character(path_to_load+'counts.txt')

    df_counts_.index.to_series().to_csv(path_to_load+'genes.txt', header = False, index = False)

def parse_list(string):
    # Remove spaces and split comma-separated values into a list.
    return [s.strip() for s in string.split(',')]

# main
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("project_list_name", 
                    type=parse_list,
                    help="List of Names of the TCGA project in GDC, example 'TCGA-BRCA,TCGA-LUAD'")
    parser.add_argument("--path_to_load", type=str, help="Path to load data")
    parser.add_argument("--path_to_download", type=str, help="Path to download data")

    args = parser.parse_args()
    project_list_name = args.project_list_name

    for project_name in project_list_name:
        print(project_name)
        project_name_2 = project_name.split('-')[1].lower()
        print(f"Load data for {project_name} ({project_name_2})")

        if args.path_to_load:
            path_to_load = _path_with_trailing_slash(args.path_to_load)
        else:
            path_to_load = _path_with_trailing_slash(PROJECT_ROOT / "data" / project_name / "treated_files")
        if args.path_to_download:
            path_to_download = _path_with_trailing_slash(args.path_to_download)
        else:
            path_to_download = _path_with_trailing_slash(PROJECT_ROOT / "data" / project_name / "download_data")

        os.makedirs(path_to_load, exist_ok=True)
        os.makedirs(path_to_download, exist_ok=True)
        
        # Load RNA-seq data
        rna_seq_files = get_rnaseq_data(project_name)
        print(f"Number of  RNA-seq files found: {len(rna_seq_files)}")

        if not os.path.exists(path_to_download+'RNA_seq'):
            os.makedirs(path_to_download+'RNA_seq')
        
        for file in rna_seq_files:
            if file['file_name'] not in os.listdir(path_to_download+'RNA_seq'):
                download_file(file['file_id'], path_to_download+ f'RNA_seq/{file["file_name"]}')
        print(f"RNA-seq files downloaded to {path_to_download}RNA_seq/")

        combined_counts_2, df_ = get_counts_data(rna_seq_files, path_to_download = path_to_download)
        
        # clinical data
        download_clinical_data(project_name, path_to_download)
        df_clinical_data_ = load_clinical_data(project_name, df_, path_to_load,  path_to_download)
        
        load_counts_data(combined_counts_2, df_clinical_data_, path_to_load)

        # Remove the downloaded files (not a good idea because sometimes you have HTTP error)
        if 0:
            try:
                shutil.rmtree(path_to_download)
                print(f"Success for removing the downloaded files")
            except FileNotFoundError:
                print(f"{path_to_download} does not exist")
            except Exception as e:
                print(f"Error: {e}")
