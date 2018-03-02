# UAMS group: Dr. Donald Johann Jr, Dr. Erich Peterson, Jason Liem

import sqlite3
import csv
import os
import sys
import bpa_analysis_functions_uams as ar 

def create_db_tables(path, drop_if_exists=True):
    """Creates mandatory sqlite3 DB tables and returns the location of the sqlite3 DB."""

    # set filename for sqlite3 DB
    db_full_location = os.path.join(os.path.dirname(os.path.realpath(__file__)), path,
                                           'vcf.sqlite')
    
    if drop_if_exists or not os.path.exists(db_full_location):
        # connect to DB
        conn = sqlite3.connect(db_full_location)
        cursor = conn.cursor()

        # drop table if exists
        drop_tables_str = '''
            DROP TABLE IF EXISTS annotations;
            DROP TABLE IF EXISTS vcf;
            DROP TABLE IF EXISTS metadata;
        '''

        print('Dropping tables if they exist.')
        cursor.executescript(drop_tables_str)
        conn.commit()

        # create table in DB called vcf to store VCF and snpEff annotations data
        create_table_str = '''
           CREATE TABLE vcf (
                vcf_id INTEGER PRIMARY KEY,
                chrom TEXT,
                pos TEXT,
                id TEXT,
                ref TEXT,
                alt TEXT,
                qual TEXT,
                filter TEXT,
                info TEXT,
                format TEXT,
                normal TEXT,
                normal_dp INTEGER,
                normal_af REAL,
                tumor TEXT,
                tumor_dp INTEGER,
                tumor_af REAL,
                metadata_id INTEGER
            );
            CREATE TABLE annotations (
                ensembl_gene_id TEXT,
                effect TEXT,
                effect_impact TEXT,
                amino_acid_change TEXT,
                gene_symbol TEXT,
                dna TEXT,
                vcf_id INTEGER,
                FOREIGN KEY(vcf_id) REFERENCES vcf(vcf_id)
            );
            CREATE TABLE metadata (
                metadata_id INTEGER PRIMARY KEY,
                case_name TEXT,
                case_id INTEGER,
                experimental_strategy TEXT,
                file_name TEXT
            );
            '''

        print('Creating tables.')
        cursor.executescript(create_table_str)
        conn.commit()

        conn.close()

    return db_full_location


def import_vcf_data(db_full_location, path, vcf_file, metadata_id, import_only_pass=True):
    """Imports data for a given VCF file into the sqlite3 DB."""

    # connect to DB
    conn = sqlite3.connect(db_full_location)

    full_file_name = path + vcf_file
    
    print('\nOpening: ' + full_file_name)
    with open(full_file_name, 'rb') as vcf_fp:  # open vcf_file
        vcf_reader = csv.reader(vcf_fp, delimiter='\t')  # create csv reader

        print('Importing: ' + vcf_file)
        file_row_count = 0
        insert_row_count = 0
        # advance to first record (skipping header info)
        # and insert rows of data into vcf table
        for row in vcf_reader:
            if '#' == row[0][0]:
                continue

            if row[6] != 'PASS' and import_only_pass:
                continue

            file_row_count += 1
            if file_row_count % 50 == 0:
                print('Total File Records Read: ' + str(file_row_count) + '\r'),
                sys.stdout.flush()
            insert_row_count += insert_vcf_row(conn, row, metadata_id)

    print('Imported: ' + vcf_file + ' | Total File Records in File : ' + str(file_row_count) +
          ' | Total Rows Inserted: ' + str(insert_row_count))

    if file_row_count != insert_row_count:
        print('ERROR: Rows in file did not equal rows inserted.')

    conn.close()


def insert_vcf_row(conn, row, metadata_id):
    """Inserts one row of VCF data into the vcf table and snpEff data into the annotations table (if present)."""

    cursor = conn.cursor()

    params = ()

    for item in row:
        if item.strip() != '':
            params = params + (item, )

    params = params + (metadata_id, )

    if len(params) == 12:  # has two sample normal and tumor (IN THAT ORDER)
        # parse out tag values for DP and AF
        normal_dp = get_tag_value("DP", row[8], row[9])
        normal_af = get_tag_value("AF", row[8], row[9])
        tumor_dp = get_tag_value("DP", row[8], row[10])
        tumor_af = get_tag_value("AF", row[8], row[10])

        # add tag values to params
        params = params + (normal_dp, normal_af, tumor_dp, tumor_af, )

        cursor.execute('INSERT INTO vcf '
                       '(chrom, pos, id, ref, alt, qual, filter, info, format, '
                       'normal, tumor, metadata_id, '
                       'normal_dp, normal_af, tumor_dp, tumor_af) '
                       'VALUES '
                       '(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', params)
    elif len(params) == 11:  # has only one sample (the tumor only)
        # parse out tag values for DP and AF
        tumor_dp = get_tag_value("DP", row[8], row[9])
        tumor_af = get_tag_value("AF", row[8], row[9])

        # add tag values to params
        params = params + (tumor_dp, tumor_af,)
        cursor.execute('INSERT INTO vcf '
                       '(chrom, pos, id, ref, alt, qual, filter, info, format, tumor, metadata_id, '
                       'tumor_dp, tumor_af) '
                       'VALUES '
                       '(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)', params)
    else:
        print('ERROR: An unexpected number of columns were found in the vcf file. Exiting.')
        sys.exit(1)

    conn.commit()

    last_rowid = cursor.lastrowid

    # if an snpEff annotation is found insert it into the annotations table
    for tag in row[7].split(';'):
        if tag.startswith("ANN="):  # if the ANN tag is found in the INFO column
            predicted_effects = tag[4:].split(',')
            for effect in predicted_effects:
                annotation_row = effect.split('|')
                if len(annotation_row) >= 10:  # corner case where there was not enough columns
                    insert_annotation_row(conn, annotation_row, last_rowid)
        else:
            continue

    return 1


def insert_annotation_row(conn, annotation_row, last_rowid):
    """Inserts one row of snpEff annotation data into the annotations table."""
    cursor = conn.cursor()

    cursor.execute('INSERT INTO annotations '
                   '(ensembl_gene_id, effect, effect_impact, amino_acid_change, '
                   'gene_symbol, dna, vcf_id) '
                   'VALUES '
                   '(?, ?, ?, ?, ?, ?, ?)', (annotation_row[4], annotation_row[1], annotation_row[2],
                                             annotation_row[10], annotation_row[3], annotation_row[9],
                                             last_rowid))
    conn.commit()


def insert_metadata_row(db_full_location, case_name, case_id, experimental_strategy, file_name):
    """Inserts one row of metadata into the metadata table."""
    
    # connect to DB
    conn = sqlite3.connect(db_full_location)

    cursor = conn.cursor()

    # verify that this data has not already been imported by checking for matches in metadata table
    cursor.execute('SELECT metadata_id, case_name, case_id, experimental_strategy, file_name '
                   'FROM metadata '
                   'WHERE case_name = ? '
                   'AND case_id = ? '
                   'AND experimental_strategy = ? ', 
                   (case_name, case_id, experimental_strategy))
    
    rows = cursor.fetchall()
    
    # no match on case_name, case_id, experimental_strategy
    if len(rows) == 0:
        
        cursor.execute('INSERT INTO metadata '
                       '(case_name, case_id, experimental_strategy, file_name) '
                       'VALUES (?, ?, ?, ?)', (case_name, case_id, experimental_strategy, file_name))
        conn.commit()

        print("Row inserted into metadata table.")
        
        return cursor.lastrowid
        
    # BID in the file name indicates files belong to the same study.
    # Parse the file_name field (csv list of file names) to see if any file's BID matches the
    # BID of the file the user is attempting to import; if so we will append this file to the csv
    else:
        # get BID from file_name
        #file_name_substring = file_name.split("BID-", 1)[1]
        #file_name_bid = file_name_substring.split("-", 1)[0]
        
        file_name_bid = (file_name.split("BID-", 1)[1]).split("-", 1)[0]
        
        for row in rows:

            file_list = str(row[4]).split(",")
            
            # file already exists, abort
            if file_name in file_list:
                print("File already imported.  Abort insert to metadata table.")
                return -1            
          
            else:
                for f in file_list:
                    
                    #f_substring = f.split("BID-", 1)[1]
                    #f_bid = f_substring.split("-", 1)[0]
                    
                    f_bid = (f.split("BID-", 1)[1]).split("-", 1)[0]

                    # file found with a matching BID , append this file_name 
                    # to the file_name column as csv 
                    if file_name_bid == f_bid:    
                        updated_file_name = row[4] + "," + file_name

                        cursor.execute('UPDATE metadata '
                                       'SET file_name = ? '
                                       'WHERE metadata_id = ? ', (updated_file_name, row[0]))

                        conn.commit()

                        print("Metadata table row " + str(row[0]) + " updated.")

                        return row[0]
          
            # no matches found from rows with the same case and experimental strategy with the same BID
            # insert a new record
            cursor.execute('INSERT INTO metadata '
                       '(case_name, case_id, experimental_strategy, file_name) '
                       'VALUES (?, ?, ?, ?)', (case_name, case_id, experimental_strategy, file_name))
                               
            conn.commit()

            print("Row inserted into metadata table.")
            
            return cursor.lastrowid
        
        
def get_tag_value(tag_string, format_string, sample_string):
    """Returns the value of the tag in a given sample."""
    split_format_string = format_string.split(':')

    # find the position of the tag specified
    tag_index = 0
    found_tag = False
    for tag in split_format_string:
        if tag == tag_string:
            found_tag = True
            break
        tag_index += 1

    if found_tag:
        split_sample_string = sample_string.split(':')

        return split_sample_string[tag_index]
    else:
        return None
        
        
def import_case_sqlite(db_full_location, project, profile, path, case_dict, case):
    case_id = case_dict.get(case)
    case_vcfs = ar.dict_VCF_files_by_case(project, case)

    for vcf, strategy in case_vcfs.items():
        if strategy == 'Panel': # or strategy == 'Total RNA':
            metadata_id = insert_metadata_row(db_full_location, case, case_id, strategy, vcf)

            if metadata_id >= 0:
                ar.download_file(project, profile, path, vcf)
                extracted_file_name = vcf[:-len('.gz')]
                import_vcf_data(db_full_location, path, extracted_file_name, metadata_id)
            else:
                  print("File already imported.  Abort insert to vcf & annotations table.")
            