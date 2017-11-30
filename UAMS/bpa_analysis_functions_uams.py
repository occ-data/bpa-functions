# UAMS group: Dr. Donald Johann Jr, Dr. Erich Peterson, Jason Liem

import subprocess
import os
import sys
import pysam
from pysam import VariantFile
import gzip
import shutil
import csv
import bpa_analysis_functions_v2 as bp

pysam.set_verbosity(0)

bp.add_keys('/home/ubuntu/.secrets')

def download_file(project, profile, path, file_name):
    ''' Automates downloading files and compressed files (and extracting) if a local copy does not exist '''
    
    # file is not compressed
    if not file_name.endswith('.gz'):        

        # download file if it does not exist locally
        if os.path.isfile(path + file_name):
            print("Local file copy found: " + file_name)
        else:
            print("Local file copy not found, downloading: " + file_name)
            bp.get_files_from_bucket(project, profile, path, file_name)

    # file is compressed
    else:
        unzipped_file_name = file_name[:-len('.gz')]

        # extracted local copy found
        if os.path.isfile(path + unzipped_file_name):
            print("Local extracted file copy found: " + unzipped_file_name)
        else:
            # compressed local copy found
            if os.path.isfile(path + file_name):
                print("Local compressed file copy found, extracting: " + file_name)
                extract_gz_file(path + file_name)
            # nothing exists locally, download and extract
            else:
                print("Downloading file: " + file_name)
                bp.get_files_from_bucket(project, profile, path, file_name)
                print("Extracting file: " + file_name)
                extract_gz_file(path + file_name)

    print('-'*80)

    
def download_gz_file_old(project, profile, path, file_name_gz):
    ''' Checks if specified .gz has an extracted version locally (ready to use).  
        If .gz file is found but not extracted yet, file is extracted. 
        If .gz file is not found, download file and extract it.'''
    
    if file_name_gz.endswith('.gz'):
        
        file_name = file_name_gz[:-len('.gz')]
        
    else:
        file_name = file_name_gz
        
    path_and_file_name = path + file_name
    path_and_file_name_gz = path + file_name_gz
        
    # check if extracted file exists
    if os.path.isfile(path_and_file_name):
        print("Extracted file found: " + file_name)
    else: 
        # compressed file found but not extracted; extract file
        if os.path.isfile(path_and_file_name_gz):
            print("Compressed file found, extracting: " + file_name_gz)
            extract_gz_file(path_and_file_name_gz)
        else:
            # compressed and extracted file not found; download and extract file
            print("File not found, downloading: " + file_name_gz)
            bp.get_files_from_bucket(project, profile, path, file_name_gz)
            print("Extracting file: " + file_name_gz)
            extract_gz_file(path_and_file_name_gz)
            
    print('-'*80)
            
        
def extract_gz_file(file_path):
    ''' Extract a *.gz file '''
    
    output_path = file_path[0:(file_path.rfind('.'))]

    with gzip.open(file_path, 'rb') as f_in, open(output_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
        
    print(file_path + " extracted to " + output_path)
    
    os.remove(file_path)
    
    print(file_path + " removed")

    
def query_somatic_mutations_by_case(project_id, case_id):
    ''' Query all VCFs and it's experimental strategy for a specific case'''
    
    query_txt = """ { project(project_id: "%s") {
                        studies {
                            cases(submitter_id: "%s") { 
                                biospecimens {
                                  samples {
                                    aliquots {
                                        analytes {
                                        read_groups {
                                          submitted_somatic_mutations {
                                            experimental_strategy,
                                            file_name
                                          }
                                        }
                                      } 
                                    }
                                  }
                                }
                            }
                        }
                      }
                    }
                """ % (project_id, case_id)
    
    data = bp.query_api(query_txt)
    
    return data

    
def query_cases_by_project(project_id):
    ''' Query all patients in a specific project '''
    
    query_txt = """ { project(project_id: "%s") {
                        studies {
                             cases {
                                submitter_id
                                
                            }
                        }
                      }
                    }
                    """ %(project_id)
    
    data = bp.query_api(query_txt)
    
    return data


def print_cases_by_project(project_id, path, file_name):
    ''' Print to screen and print to file a list of all patients in a project '''
    
    cases = query_cases_by_project(project_id)

    cases_file = open(path + file_name, 'w')
    
    for c in cases["data"]["project"][0]["studies"][0]["cases"]:
        print(c["submitter_id"])
        cases_file.write("%s\n" %c["submitter_id"])
        
    cases_file.close()
    
    
def assign_case_id_dict(project_id):
    
    cases = query_cases_by_project(project_id)
    
    cases_dict = {}
    
    case_id = 1
    
    for c in cases["data"]["project"][0]["studies"][0]["cases"]:
            cases_dict[c["submitter_id"]] = case_id
            case_id += 1
            
    return cases_dict
    
    
def dict_VCF_files_by_case(project_id, case_id):
    ''' Return dictionary of all VCF files (key) and experimental strategy (value) for a specific case'''
    
    data = query_somatic_mutations_by_case(project_id, case_id)

    vcf_dict = {}
    
    for project in data["data"]["project"]:
        for studies in project["studies"]:
            for cases in studies["cases"]:
                for biospecimens in cases["biospecimens"]:
                    for samples in biospecimens["samples"]:
                        for aliquots in samples["aliquots"]:
                            for analytes in aliquots["analytes"]:
                                for read_groups in analytes["read_groups"]:
                                    for ssm in read_groups["submitted_somatic_mutations"]:
                                        vcf_dict[ssm["file_name"]] = ssm["experimental_strategy"]

    return vcf_dict


def print_VCF_files_by_case(project, case, path, file_name):
    ''' Print to screen and print to file a list of all VCFs and experimental strategies for specific case '''

    vcf_dict = dict_VCF_files_by_case(project, case)
    
    if vcf_dict:
        case_vcfs_file = open(path + file_name, 'w')

        for vcf,strategy in vcf_dict.items():
            print(vcf + "\t" + strategy)
            case_vcfs_file.write("%s\t%s\n" %(vcf, strategy))

        case_vcfs_file.close()
        
    else:
        print("No files found")
    
    
def somatic_mutation_by_gene_and_mutation(project, profile, path, vcf_gz, gene, mutation, details = False):
    ''' Return all records in a specified VCF that contains a specified gene and mutation and the total found'''
    
    download_file(project, profile, path, vcf_gz)
    
    if vcf_gz.endswith('.gz'):
        vcf = vcf_gz[:-len('.gz')]
    else:
        vcf = vcf_gz
    
    vcf_reader = VariantFile(path + vcf, drop_samples = True)
    
    counter = 0

    print("Processing file: " + vcf + "\n")

    # iterate through each record
    for rec in vcf_reader.fetch():
        if "ANN=" in str(rec):

            # get list of Annotations
            ann_list = str(rec.info["ANN"]).replace(" ", "").split(",")

            found = False
            
            # parse each Annotation
            for a in ann_list:
                ann = a.split("|")

                if len(ann) >= 10:
                    if ann[3] == gene and ann[10] == mutation:
                        
                        # first occurrence of a match in list of Annotations; print entire record
                        if found is False:
                            found = True
                            print(str(rec.info["ANN"]) + "\n")
                            counter += 1
                        
                        # print significant fields from record for easier visibility
                        if details is True:
                            print("* gene: %s, chromosome: %s, mutation offset: %s, nucleotide change: %s, AA change: %s\n" 
                                  % (ann[3], rec.chrom, rec.pos, ann[9], ann[10]))
    
    print('-'*80)
    print("Records found: " + str(counter))
    print("Finished processing file: " + vcf)
    print('-'*80)
    
    return counter


def somatic_mutation_by_gene_and_mutation_list(project, profile, path, vcf_gz_csv, gene, mutation, details = False):
    ''' For a CSV string list of VCFs return all records per VCF that contains a specified 
        gene and mutation and aggregate total found '''
    
    # convert vcf__gz_csv to unique list (object) of vcf_gz file names
    vcf_gz_list = vcf_gz_csv.replace(" ", "").split(",")
    vcf_gz_list = list(set(vcf_gz_list))

    counter = 0
    
    for vcf_gz in vcf_gz_list:
        counter += somatic_mutation_by_gene_and_mutation(project, profile, path, vcf_gz, gene, mutation, details)

    return counter


def all_somatic_mutation_by_gene_and_mutation_by_case(project, profile, path, case, gene, mutation, details = False):
    
    ''' For a specificed patient, return all records for all VCFs that 
        contains a specified gene and mutation and the aggregate total found '''
    
    vcf_dict = dict_VCF_files_by_case(project, case)
    
    vcf_gz_list = list(vcf_dict.keys())
    
    vcf_gz_csv = ",".join(vcf_gz_list)    
    
    print("VCF files found for case %s: \n" %case)
    
    for vcf,strategy in vcf_dict.items():
        print(vcf + "\t" + strategy)
        
    print('-'*80)
    
    return somatic_mutation_by_gene_and_mutation_list(project, profile, path, vcf_gz_csv, gene, mutation, details)


