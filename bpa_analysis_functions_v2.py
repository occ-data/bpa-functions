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
from operator import add

pysam.set_verbosity(0)

auth = ''

main_header_order = [
    'Sample',
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


class arrayTable(list):
    ''' Represent result arrays in HTML format for visualization '''

    def _repr_html_(self):
        html = []
        html.append("<table style>")
        for value in self:
            html.append("<tr>")
            html.append("<td>%s</td>" % value)
            html.append("<tr>")
        html.append("</table>")

        return ''.join(html)

class SummaryTable(dict):
    ''' Represent result tables in HTML format for visualization '''

    def _repr_html_(self):
        html = []
        html.append("<table style>")
        html.append("<thead>")
        headers = []
        for key in self:
            for field in self[key]:
                if type(field) is dict:
                    if field not in headers:
                        headers.append(field)
                        html.append("<th>%s</th>" % (field))
        if not headers:
            html.append("<th>%s</th>" % (main_header_order[0]))
            html.append("<th>%s</th>" % (main_header_order[1]))
        html.append("</thead>")
        for key in self:
            if headers:
                html.append("<tr>")
                html.append("<td>%s</td>" % key)
                for h in headers:
                    if h in self[key]:
                        html.append("<td>%s</td>" % str(self[key][h]))
                    else:
                        html.append("<td>0</td>")
                html.append("</tr>")
            else:
                for value in self[key]:
                    html.append("<tr>")
                    html.append("<td>%s</td>" % key)
                    html.append("<td>%s</td>" % str(value))
                    html.append("</tr>")
        html.append("</table>")

        return ''.join(html)


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

    json_data=open(filename).read()
    keys = json.loads(json_data)
    auth = requests.post('https://data.bloodpac.org/user/credentials/cdis/access_token', json=keys)


def get_files_from_bucket(project, profile, files_path, files=None):
    ''' Transfer data from object storage to the VM in the private subnet '''

    # Create folder
    if not os.path.exists(files_path):
       os.makedirs(files_path)

    # Get bucket name and path
    bucket_name = project.replace('bpa-', 'BPA_')
    s3_path = 's3://bpa-data/' + bucket_name

    # If only one file or a pattern, create array
    if isinstance(files, str):
        files=[files]

    # Getting files
    print("Getting files...")
    if files:
       for f in files:
          s3_path = s3_path + '/'
          cmd = ['aws', 's3', 'cp', s3_path, files_path, '--recursive', '--profile', profile, '--exclude', '*', '--include', f]
          try:
              output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
          except Exception as e:
              output = e.output
              print("ERROR:", output)
    else:
       cmd = ['aws', 's3', 'sync', s3_path, files_path, '--profile', profile]
       try:
          output = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
       except Exception as e:
          output = e.output
          print("ERROR:", output)
    print("Finished")


def query_api(query_txt, variables = None):
    ''' Request results for a specific query '''

    if variables == None:
        query = {'query': query_txt}
    else:
        query = {'query': query_txt, 'variables': variables}

    output = requests.post('https://data.bloodpac.org/api/v0/submission/graphql/' ,
                           headers={'Authorization': 'bearer '+ auth.json()['access_token']},
                           json=query).text
    data = json.loads(output)

    if 'errors' in data:
        print(data)

    return data


def query_project_samples(project_id):
    ''' Query samples for a specific project'''

    query_txt = """query Test { sample (first:0, project_id: "%s") {
                               submitter_id}} """ % (project_id)

    data = query_api(query_txt)

    return data


def query_sample(project_id, sample_id):
    ''' Query alignment files from one sample'''

    query_txt = """{ sample (project_id: "%s", submitter_id: "%s") {
                               submitter_id
                               _aliquots_count
                               aliquots {
                               analytes {
                                  _read_groups_count
                                  read_groups {
                                     _submitted_somatic_mutations_count submitted_somatic_mutations { file_name}
                                     _submitted_unaligned_reads_files_count submitted_unaligned_reads_files { file_name}
                                     _submitted_aligned_reads_files_count submitted_aligned_reads_files { file_name}
                                     _submitted_copy_number_files_count submitted_copy_number_files { file_name}
                                  }
                               }
                               }
                         }
                    } """ % (project_id, sample_id)

    data = query_api(query_txt)

    return data

def query_field_counts(node, fields):
    ''' Query summary counts for each data type '''

    query_txt = """{ case(first:0, with_path_to: {type: "%s" """ % (node)
    for f in fields:
        query_txt += """,%s: "%s" """ % (f, fields[f])
    query_txt += """}){ project_id }}"""

    data = query_api(query_txt)

    summary = {}
    if 'data' in data and 'case' in data['data']:
        for d in data['data']['case']:
            project = d["project_id"]
            if project not in summary:
                summary[project] = {}

            summary[project].setdefault('COUNT', 0)
            summary[project]['COUNT'] += 1

    return SummaryTable(summary)


def query_summary_field(node, field, project_id=None):
    ''' Query summary counts for each data type '''

    if project_id != None:
        query_txt = """query { %s(first:0, project_id: "%s") {%s}} """ % (node, project_id, field)
    else:
        query_txt = """query { %s(first:0) {%s project_id}} """ % (node, field)

    data = query_api(query_txt)

    summary = {}
    total = []
    for d in data['data'][node]:

        if isinstance(d[field], float):
            d[field] = str(d[field])[:-2]

        if 'project_id' in d:
            summary.setdefault(d['project_id'], {})
            summary[d['project_id']].setdefault(d[field], 0)
            summary[d['project_id']][d[field]] += 1
            if d[field] not in total:
                total.append(d[field])
        else:
            summary.setdefault(d[field], 0)
            summary[d[field]] += 1

    #plot_summary(summary, field)

    if project_id != None:
        plot_field_metrics(summary, field)
    else:
        plot_overall_metrics(summary, field, total)

    return summary


def plot_field_metrics(summary_counts, field):
    ''' Plot summary results in a barplot '''

    N = len(summary_counts)

    values = []
    types = []
    for n in sorted(summary_counts):
        value = 0
        for p in summary_counts[n]:
            value += summary_counts[n][p]
        values.append(value)
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

def plot_overall_metrics(summary_counts, field, totals):
    ''' Visualize summary results across projects in a barplot '''

    results = {}
    projects = {}
    for project in summary_counts:

        results[project] = []
        projects.setdefault(project, 0)

        for value in totals:
            if value in summary_counts[project]:
                results[project].append(summary_counts[project][value])
                projects[project] += summary_counts[project][value]
            else:
                results[project].append(0)

    N = len(totals)
    positions = np.arange(N)
    sorted_projects = sorted(projects, key=projects.get, reverse=True)
    bar_size = 0.4
    size_prop = (N/10) + 1

    plots = []
    plt.figure(figsize=(2*N, N))
    left = [0]*N
    for pr in sorted_projects:
        p = plt.barh(positions, results[pr], bar_size, left, align='center', alpha=1)
        plots.append(p[0])
        left = map(add, left, results[pr])

    plt.title('Summary counts by (' + field + ')', fontsize=10*size_prop)
    plt.xlabel('COUNTS', fontsize=10*size_prop)
    plt.ylabel(field.upper(), fontsize=10*size_prop)
    plt.yticks(positions, totals, fontsize=10*size_prop)
    plt.legend(plots, sorted_projects, fontsize=10*size_prop)

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


def query_expectations(project_id, vcf_name):
    ''' Retrieve all expected mutations associated to one VCF in one project'''

    query_txt = """{
                   aliquot(project_id: "%s", with_path_to: {type: "submitted_somatic_mutation", file_name: "%s"}) {
                       _contrived_expectations_count
                       samples{
                           submitter_id
                       }
                       contrived_expectations(first:0) {
                          expected_mutation_chromosome
                          expected_mutation_position
                       }
                   }
                }""" % (project_id, vcf_name)

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
            for an in a['analytes']:
                for rg in an['read_groups']:
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


def get_expected_mutations(project_id, vcf_name):
    ''' Retrieve expected mutation from an expectation query '''


    data = query_expectations(project_id, vcf_name)

    expectations = []
    for a in data["data"]["aliquot"]:
        sample_id = a["samples"][0]["submitter_id"]
        for se in a["contrived_expectations"]:
            expectation = {'sample_id': sample_id, 'vcf': vcf_name,
                           'expected_mutation_chromosome': se['expected_mutation_chromosome'].encode('ascii').replace('chr', ''),
                           'expected_mutation_position': se['expected_mutation_position'].encode('ascii')}
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


def calculate_metrics_vcf(project, path, vcf_name, baseline_vcf=None):
    ''' Calculate sensitivity/specificity for one VCF file and its corresponding expectations '''

    data = {'Sample': '', 'VCF File': '', 'Expectations': 0, 'True-Positive': 0, 'False-Positive': 0, 'Sensitivity': 0.0 , 'Specificity': 0.0}
    vcf_path = path + vcf_name
    vcf_in = pysam.VariantFile(vcf_path, 'rb')

    expectations = get_expected_mutations(project, vcf_name)
    if not expectations:
       print("Warning: There are no expected mutations for %s VCF file" % vcf_name)
       return {}

    if baseline_vcf:
       expectations = find_germlines(expectations, path + baseline_vcf)

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
           and var['expected_mutation_position'] == pos for var in expectations):
              TP += 1
        else:
              FP += 1

    sample_id = [var['sample_id'] for var in expectations][0]
    P  = len(expectations)
    TN = 169 - (TP + FP)
    data['Sample'] = sample_id
    data['VCF File'] = vcf_name
    data['Expectations'] = P
    data['True-Positive'] = TP
    data['False-Positive'] = FP
    data['Sensitivity'] = round(float(TP)/float(P), 3)
    data['Specificity'] = round(float(TN)/float(TN+FP),3)

    return MetricsTable([data])

def calculate_metrics_all_vcf(project, path, vcfs_files, samples=None, baseline_vcf=None):
    ''' Calculate sensitivity/specificity for a set of VCF files and create a table '''

    data_results = []

    for sample in vcfs_files:
       if sample in samples:
           for vcf in vcfs_files[sample]:
              data = calculate_metrics_vcf(project, path, vcf, baseline_vcf)
              if data:
                 data_results = data_results + data

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
