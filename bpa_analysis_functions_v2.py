from cdispyutils.hmac4 import get_auth
import subprocess
import glob
import os
import sys
import requests
import json
import pysam
import numpy as np
import matplotlib.pyplot as plt

pysam.set_verbosity(0)

auth = ''

main_header_order = [
    'VCF File',
    'Expectations',
    'True-Positive',
    'False-Positive',
    'Sensitivity',
    'Specificity'
]

data_types = {
    'VCF': 'submitted_somatic_mutations',
    'FASTQ': 'submitted_unaligned_reads_files',
    'BAM': 'submitted_aligned_reads_files',
    'CNV': 'submitted_copy_number_files' 
}

metadata_types = {
    'METADATA': 'experiment_metadata_files'
}


class MetricsTable(list):
    ''' Represent result tables in HTML format for visualization '''
    
    def _repr_html_(self):
        html = []
        html.append("<table style>")
        html.append("<thead>")
        for key in main_header_order:  
            html.append("<th>%s</th>" % key)
        html.append("</thead>")       
        for line in self:
            html.append("<tr>") 
            for key in main_header_order:  
               html.append("<td>%s</td>" % line[key])
            html.append("<tr>") 
        html.append("</table>")        
        
        return ''.join(html)


def add_keys(filename):
    ''' Get auth from our secret keys '''

    global auth 
    with open(filename,'r') as f:
        secrets = json.load(f)
    auth = get_auth(secrets['access_key'], secrets['secret_key'], 'submission')


def get_files_from_bucket(project, profile, files_path, include_type=None):
    ''' Transfer data from object storage to the VM in the private subnet '''

    # Create folder
    if not os.path.exists(files_path):
       os.makedirs(files_path)

    # Get bucket name and path
    bucket_name = project.replace('bpa-', 'BPA_')
    s3_path = 's3://bpa-data/' + bucket_name

    # Getting files
    print "Getting files..."
    cmd = ['aws', 's3', 'sync', s3_path, files_path, '--profile', profile]
    if include_type:
       cmd = cmd + ['--exclude', '*', '--include' , include_type]
    
    try:
        output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)       
    except Exception as e:
        output = e.output
        print "ERROR:" + output
    print "Finished"


def query_api(query_txt):
    ''' Request results for a specific query '''

    query = {'query': query_txt}
    output = requests.post('http://api.internal.io/v0/submission/graphql/', auth=auth, json=query).text
    data = json.loads(output) 

    if 'data' not in data and 'errors' in data:
        print query    
    
    return data


def query_project_samples(project_id):
    ''' Query samples for a specific project'''

    query_txt = """query Test { sample (first:1000, project_id: "%s") {
                               submitter_id}} """ % (project_id)

    data = query_api(query_txt)   

    return data


def query_sample(project_id, sample_id):
    ''' Query alignment files from one sample'''

    query_txt = """query Test { sample (project_id: "%s", submitter_id: "%s") {
                               submitter_id _aliquots_count aliquots {
                               aliquot_concentration  _read_groups_count read_groups {
                               _submitted_somatic_mutations_count submitted_somatic_mutations { file_name} 
                               _submitted_unaligned_reads_files_count submitted_unaligned_reads_files { file_name} 
                               _submitted_aligned_reads_files_count submitted_aligned_reads_files { file_name}
                               _submitted_copy_number_files_count submitted_copy_number_files { file_name}} } }} """ % (project_id, sample_id)

    data = query_api(query_txt)   

    return data


def query_summary_field(node, field, project_id=None):
    ''' Query summary counts for each data type '''
    
    if project_id == None:
        query_txt = """ query {project{project_id}}"""
        projects = query_api(query_txt)['data']['project']    
    else:
        projects = [{"project_id": project_id}]
    
    summary = {}
    for p in projects:
               
        query_txt = """query { %s(first:0, project_id: "%s") {%s}} """ % (node, p["project_id"], field) 
        
        data = query_api(query_txt)

        name = p["project_id"]
        summary[name] = {}
        for d in data['data'][node]:
            
            summary[name].setdefault(str(d[field]), 0)        
            summary[name][str(d[field])] += 1
    
    # plot_summary(summary, field, project_id)
    
    return summary


def plot_summary(summary_counts, field, project_id):
    ''' Plot summary results in a barplot ''' 
    
    N = len(summary_counts[project_id])

    values = []
    types = []
    for n in sorted(summary_counts[project_id]):
        values.append(summary_counts[project_id][n])
        types.append(n)
        
    positions = np.arange(N)        
    plt.figure(figsize=(2*N, N))   
    
    size_prop = (N/10) + 1
    plt.barh(positions, values, 0.2, align='center', alpha=0.5, color='b')
    plt.title('Summary counts by (' + field + ')', fontsize=10*size_prop)
    plt.xlabel('COUNTS', fontsize=10*size_prop)    
    plt.ylabel(field.upper(), fontsize=10*size_prop)  
    plt.yticks(positions, types, fontsize=10*size_prop)    
    
    for i, v in enumerate(values):
        plt.text(v, i, str(v), color='red', fontweight='bold', fontsize=10*size_prop)
    
    plt.show()
    

def list_samples(project_id):
    ''' Retrieve samples included in one specific project'''

    sample_data = query_project_samples(project_id)

    samples = []
    for s in sample_data["data"]["sample"]:
      samples.append(s['submitter_id'].encode('ascii'))

    return samples       


def query_experimental_metadata(project_id):
    ''' Query experimental metadata files from a specific project '''


    query_txt = """query Test { experiment (project_id: "%s") {   
                               experiment_metadata_files{file_name}}} """ % (project_id)

    data = query_api(query_txt)   

    return data


def query_project(project_id):
    ''' Retrieve all sample data for one specific project '''


    data = query_project_samples(project_id)
    for s in data['data']['sample']:
      sample = query_sample(project_id, s['submitter_id'])
      s.update(sample['data']['sample'][0])
 
    return data


def query_expectations(project_id, sample_id=None):
    ''' Retrieve all expected mutations for one sample or one project'''

    if sample_id:
      query_txt = """query Test { 
                          sample (project_id: "%s", submitter_id: "%s") { 
                              submitter_id 
                              _sample_expectations_count 
                              sample_expectations(first:100) { 
                                  expected_mutation_chromosome 
                                  expected_mutation_position 
                              }
                          }
                    }""" % (project_id, sample_id)
    else:
      query_txt = """query Test { 
                          sample (project_id: "%s") { 
                              submitter_id  
                              _sample_expectations_count
                              sample_expectations(first:100) { 
                                  expected_mutation_chromosome 
                                  expected_mutation_position 
                              }
                          }
                    }""" % (project_id)    

    data = query_api(query_txt) 

    return data


def search_files(query_data, file_type):
    ''' Retrieve file names from a sample query result'''

    node = data_types[file_type.upper()]

    files = {}

    for s in query_data["data"]["sample"]:
        sample_id = s['submitter_id'].encode('ascii')
        if not sample_id in files:
           files[sample_id] = []
        for a in s['aliquots']:
            for rg in a['read_groups']:
                for f in rg[node]:
                    if 'file_name' in f:
                       files[sample_id].append(f['file_name'].encode('ascii'))

    return files


def search_metadata(query_data, file_type):
    ''' Retrieve file names from an experimental metadata query result'''

    node = metadata_types[file_type.upper()]

    files = []
    for e in query_data["data"]["experiment"]:
        for em in e['experiment_metadata_files']:
           if 'file_name' in em:
              files.append(em['file_name'].encode('ascii'))


    return files


def list_files_by_type(project_id, file_type, sample_id=None):
    ''' Retrieve file names in the project/sample of a specific type'''

    file_type = file_type.upper()

    if file_type == 'METADATA':
      data = query_experimental_metadata(project_id)
      files = search_metadata(data, file_type)
    elif sample_id:
      data = query_sample(project_id, sample_id)
      files = search_files(data, file_type)  
    else:
      data = query_project(project_id) 
      files = search_files(data, file_type)
    
    return files   


def list_files(project_id, sample_id=None):
    ''' Retrieve all file names associated to a project/sample''' 

    if sample_id:
      data = query_sample(project_id, sample_id)  
      files = []
      for key in data_types.keys():
          type_files = search_files(data, key) 
          if type_files:
             files += type_files[sample_id]
    else:
      data = query_project(project_id)  
      files = []
      for key in data_types.keys():
          type_files = search_files(data, key) 
          if type_files:
             for sample_id in type_files:
                files += type_files[sample_id]    

      metadata = query_experimental_metadata(project_id)
      for key in metadata_types:
          type_files = search_metadata(metadata, key) 
          if type_files:
             files += type_files    

    return files


def count_file_types(project_id, sample_id=None):
    ''' Count file types associated to a project/sample ''' 


    data_files = list_files(project_id, sample_id)
    count_files = dict()
    for f in data_files:
       file_type = f.split('.')[-1]
       if file_type in count_files:
          count_files[file_type] += 1
       else:
          count_files[file_type] = 1        

    return count_files


def get_expected_mutations(project_id, sample_id=None):
    ''' Retrieve expected mutation from an expectation query ''' 


    data = query_expectations(project_id, sample_id)

    expectations = []
    for s in data["data"]["sample"]:   
        sample_id = s['submitter_id'].encode('ascii')
        for a in s['sample_expectations']:
            expectation = {'sample_id': sample_id,
                           'expected_mutation_chromosome': a['expected_mutation_chromosome'].encode('ascii').replace('chr', ''),
                           'expected_mutation_position': a['expected_mutation_position'].encode('ascii')}
            expectations.append(expectation)

    return expectations


def find_germlines(expectations, baseline):
    ''' Find potential germline variants from a baseline vcf (unexpected somatic variants) ''' 

    vcf_back = pysam.VariantFile(baseline, 'rb') 

    for rec in vcf_back.fetch():
      if 'PASS' in rec.filter:  
        chrom = rec.chrom.replace('chr', '')
        pos   = str(rec.pos)
        ref   = rec.alleles[0]
        alt   = rec.alleles[1]

        for var in expectations:
           if chrom == var['expected_mutation_chromosome'] and \
              pos == var['expected_mutation_position']:
                 expectations.remove(var)

    return expectations


def calculate_metrics_vcf(project, path, vcf_name, sample, baseline_vcf=None):
    ''' Calculate sensitivity/specificity for one VCF file and its corresponding expectations ''' 

    data = {'VCF File': '', 'Expectations': 0, 'True-Positive': 0, 'False-Positive': 0, 'Sensitivity': 0.0 , 'Specificity': 0.0}
    vcf_path = path + vcf_name
    vcf_in = pysam.VariantFile(vcf_path, 'rb') 

    expectations = get_expected_mutations(project, sample)
    if not expectations:
      print "ERROR: There are no expected mutations for this project or sample" 

    sample_expectations = [e for e in expectations if e['sample_id'] == sample]
    if baseline_vcf:
      sample_expectations = find_germlines(sample_expectations, path + baseline_vcf)  

    TP = 0
    FP = 0
    for rec in vcf_in.fetch():
      if 'PASS' in rec.filter and float(rec.info['MAF'][0]) < 0.1:  
        chrom = rec.chrom.replace('chr', '')
        pos   = str(rec.pos)
        ref   = rec.alleles[0]
        alt   = rec.alleles[1]

        # Not used yet
        if len(ref)>1:
           ref = ref[1:]
           alt = '-'
        if len(alt)>1:
           alt = alt[1:]
           ref = '-'             
  
        if any(var['expected_mutation_chromosome'] == chrom \
           and var['expected_mutation_position'] == pos for var in sample_expectations):
              TP += 1
        else:
              FP += 1

    P  = len(sample_expectations)
    TN = 169 - (TP + FP)
    data['VCF File'] = vcf_name
    data['Expectations'] = P
    data['True-Positive'] = TP
    data['False-Positive'] = FP
    data['Sensitivity'] = round(float(TP)/float(P), 3)
    data['Specificity'] = round(float(TN)/float(TN+FP),3)         
  
    return data

def calculate_metrics_all_vcf(project, path, vcfs_files, baseline_vcf=None):
    ''' Calculate sensitivity/specificity for a set of VCF files and create a table ''' 

    data_results = []

    for sample in vcfs_files: 
       for vcf in vcfs_files[sample]:       
          data = calculate_metrics_vcf(project, path, vcf, sample, baseline_vcf)
          data_results.append(data)

    table = MetricsTable(data_results)

    return table, data_results   


def plot_metrics(data_metrics, data_filter_metrics=None):
    ''' Plot table results in a barplot ''' 
    
    N = len(data_metrics)
    files         = [m['VCF File'] for m in data_metrics]
    sens_values   = [m['Sensitivity'] for m in data_metrics] 
    spec_values   = [m['Specificity'] for m in data_metrics]
    ind = np.arange(N) # the x locations for the groups
    width = 0.35       # the width of the bars
    
    if data_filter_metrics != None:
        f_sens_values = [m['Sensitivity'] for m in data_filter_metrics] 
        f_spec_values = [m['Specificity'] for m in data_filter_metrics] 
        ind = 2*ind


    fig, ax = plt.subplots(figsize=(15, 10))
    rects1 = ax.bar(ind, sens_values, width, color='#b0e0e6')
    rects2 = ax.bar(ind + width, spec_values, width, color='#87cefa')
    rects = (rects1, rects2)
    labels = ('Sensitivity', 'Specificity')
    if data_filter_metrics != None: 
        rects3 = ax.bar(ind + 2*width, f_sens_values, width, color='#4682b4')
        rects4 = ax.bar(ind + 3*width, f_spec_values, width, color='#0000cd')
        rects = (rects1, rects2, rects3, rects4)
        labels = ('Sensitivity', 'Specificity', 'Germline-filtering Sensitivity', 'Germline-filtering Specificity')
        
    # add some text for labels, title and axes ticks
    ax.set_ylabel('Metrics')
    ax.set_title('Specificity/Sensitivity Analysis')
    ax.set_xticks(ind + width)
    ax.set_xticklabels(files, rotation=90)

    ax.legend(rects, labels, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()    
