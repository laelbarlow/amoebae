#!/usr/bin/env python3
# Copyright 2018 Lael D. Barlow
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
"""Script to run a quick phylogenetic analysis with IQtree on a phylip
alignment on CIPRES without having to go through the graphical interface.

Instructions:
    1. Register and log into a CIPRES REST API account here: https://www.phylo.org/restusers/login.action

    2. From the Developer drop down menu choose "Application Management" and press the "Create New Application" button. 

    3. Complete the form:
        - Choose the default DIRECT Authentication Type
        - Make a note of the Application ID that is assigned.

    4. Define values for the following variables near the bottom of this script:
        - cra_user (the user name that you registered with in step 1)
        - application_name (the Application ID from step 3)

    5. Check that this script is executable:

        chmod a+x run_iqtree_on_cipres.py

    6. Run this script, providing a path to a phylip
    alignment file in any directory as follows:
        
        ./run_iqtree_on_cipres.py /path/to/your/phylip/file/here


Documentation for the CIPRES REST API:
    https://www.phylo.org/restusers/docs/qs.html
    https://www.phylo.org/restusers/docs/guide.html

IQ-tree-specific documentation:
    http://www.phylo.org/index.php/rest/iqtree_xsede.html

"""
import subprocess
import os
import sys
import getpass
from time import sleep
import re
import time
from string import Template



def get_job_status_url(submission_info):
    """Take submission info, and return the job status url.
    """
    line_num = 0
    for i in submission_info.split('\n'):
        line_num += 1
        if line_num == 4:
            #print(i.split('>', 1)[1].rsplit('<')[0])
            return i.split('>', 1)[1].rsplit('<')[0]


def get_job_status(job_status_url, cra_user, password, key):
    """Take submission info, and return the job status as True or False."""
    # Get info from url.
    status_info_string = str(subprocess.check_output(['curl', '-u', cra_user +\
        ':' + password, '-H', 'cipres-appkey:' + key, job_status_url]), 'utf-8')

    # Parse info.
    if ' Task is finished.' in status_info_string:
        return True
    else:
        return False


def get_output_dir_url(job_status_url, cra_user, password, key):
    """Take submission info and return job output directory url as a string."""
    # Get info from status url.
    status_info_string = str(subprocess.check_output(['curl', '-u', cra_user +\
        ':' + password, '-H', 'cipres-appkey:' + key, job_status_url]), 'utf-8')

    # Parse the xml.
    x = re.compile(r'resultsUri')
    working_dir_section = x.split(status_info_string)[1]
    y = re.compile(r'(<url>|</url>)')
    url = y.split(working_dir_section)[2]

    # Check that url was found.
    assert url is not None, """Could not find output directory url."""
    assert url.startswith('https'), """Could not find output directory url."""

    # Return the url for the output directory.
    return url


def get_output_file_urls(output_dir_url, cra_user, password, key):
    """Take an output directory url, and return the urls for all the files in
    that directory.
    """
    # Get info from directory url.
    dir_info_string = str(subprocess.check_output(['curl', '-u', cra_user +\
        ':' + password, '-H', 'cipres-appkey:' + key, output_dir_url]), 'utf-8')

    # Parse the xml.
    x = re.compile(r'(<jobfile>|</jobfile>)')
    jobfile_sections = x.split(dir_info_string)
    output_file_urls = []
    for i in jobfile_sections:
        if i.strip().startswith('<downloadUri>'):
            y = re.compile(r'(<url>|</url>)')
            url = y.split(i)[2]
            output_file_urls.append(url)

    # Return the list of output file urls.
    return output_file_urls


def download_output_files_to_dir(output_file_urls, outdir, cra_user, password,
        key):
    """Given a list of output file urls and a full path to an output directory,
    download all the output files to the directory.
    * Requires password, cra_user, and key variables to be set.
    """
    # Make the output directory if it doesn't exist.
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Loop over output file urls.
    for url in output_file_urls:
        # Download each file to the directory.
        subprocess.call(['curl', '-u', cra_user + ':' + password, '-H',\
            'cipres-appkey:' + key, '-O', '-J', url], cwd=outdir)


def manage_any_cipres_job(infile_path, submission_info, timestamp, cra_user,
        password, key):
    """Take any job submitted to CIPRES, check status regularly, and download
    all output files when complete.
    """

    # Get job status URL.
    print('\nGetting job status URL.')
    job_status_url = get_job_status_url(submission_info)
    
    # Check job status, and wait until complete.
    job_status = False
    while not job_status:
        sleep(30)
        print('\nChecking job status.\n')
        job_status = get_job_status(job_status_url, cra_user, password,
            key)
    print('\nJob complete\n')
    
    # Get output directory url.
    print('\nGetting output directory url.\n')
    output_dir_url = get_output_dir_url(job_status_url, cra_user, password, key)
    
    # Get output file urls.
    print('\nGetting output file urls.\n')
    output_file_urls = get_output_file_urls(output_dir_url, cra_user, password,
            key)
    
    # Define output directory.
    output_dir = infp + '_CIPRES_output_' + timestamp 
    
    # Download files.
    print('\nDownloading output files.\n')
    download_output_files_to_dir(output_file_urls, output_dir, cra_user,
            password, key)

    # Report location of downloaded files.
    print('\nOutput files downloaded to %s\n' % output_dir)


def submit_raxml_job_to_cipres_1(infile_path, timestamp, cra_user, password,
        key):
    """Take a phylip file and run a raxml analysis as per the quickstart
    example on the CIPRES website.
    """
    # Upload file. 
    print('\nUploading file, and initiating job\n')
    submission_info = str(subprocess.check_output(['curl',
                                                   '-u',
                                                   cra_user + ':' + password,
                                                   '-H',
                                                   'cipres-appkey:' + key,
                                                   url + '/job/' + cra_user,
                                                   '-F',
                                                   'tool=RAXMLHPC8_REST_XSEDE',
                                                   '-F',
                                                   'input.infile_=@./' + infile_path,
                                                   '-F',
                                                   'metadata.statusEmail=true'
                                                   ]), 'utf-8')

    # Manage job submission.
    manage_any_cipres_job(infile_path, submission_info, cra_user, password, key)


def submit_iqtree_job_to_cipres_1(url, cra_user, password, key, timestamp, infile_path):
    """Take a phylip file and run a basic phylogenetic analysis with IQtree on
    CIPRES.
    """
    # Determine alignment length of input alignment.
    alignment_length_from_file = None
    with open(infile_path) as f:
        alignment_length_from_file = f.readline().strip().split(' ')[1]
    assert alignment_length_from_file is not None

    ## Define path for partition file.
    ##partitionfp = infile_path.rsplit('.', 1)[0] + '_partitions.nex'
    #partitionfp = infile_path.rsplit('.', 1)[0] + '.partitions'

    #    # Write partition file, specifying only one partition for the input
    #    # alignment.
    #    template_string = Template(r"""#nexus
    #begin sets;
    #charset part1=1-$alignment_length;
    #end;
    #""")

    ## Write partition file, specifying only one partition for the input
    ## alignment.
    #template_string = Template(r"""AA, part1 = 1-$alignment_length""")

    #with open(partitionfp, 'w') as o:
    #    o.write(template_string.substitute(alignment_length=alignment_length_from_file))

    # Construct command to pass through the CIPRES REST API.
    command_list = ['curl',
                    '-u',
                    cra_user + ':' + password,
                    '-H',
                    'cipres-appkey:' + key,
                    url + '/job/' + cra_user,
                    '-F',
                    'metadata.statusEmail=true', # Send email updates.
                    '-F',
                    'tool=IQTREE_XSEDE', # Specify which tool to use.
                    '-F',
                    'input.infile_=@./' + os.path.basename(infile_path), # Specify path to main input file (phylip alignment).
                    '-F',
                    'vparam.specify_numparts_=' + '1', # Specify number of partitions.
                    #'-F',
                    #'vparam.partition_type_=' + '-spp', # Specify partition type; -q; -spp; -sp; no default.
                    #'-F',
                    #'input.partition_file_=@./' + partitionfp, # Specify path to partition file. 
                    '-F',
                    'vparam.runtime_=' + '1.0', # Specify maximum hours to run (float).
                    '-F',
                    'vparam.specify_mem_=' + '120', # Specify maximum hours to run (float).
                    #'-F',
                    #'vparam.specify_autocores_=' + '1', # Automatically select number of threads.
                    '-F',
                    'vparam.specify_runtype_=' + '2', # Specify the nrun type; 1 (Model Selection); 2 (Tree Inference); no default.
                    '-F',
                    'vparam.bootstrap_type_=' + 'bb', # (Excl) - Bootstrap Type; bb, b, bc; (default=bb)
                    '-F',
                    'vparam.num_bootreps_=' + '1000', # Specify number of bootstraps.
                    '-F',
                    'vparam.sh_test_=' + '1', # Do aLRT tests.
                    '-F',
                    'vparam.num_replicates_=' + '1000', # (Integer) - Specify number of replicates (1000 or less); no default.
                    '-F',
                    'vparam.sequence_type_=' + 'AA' # Specify sequence type (-st); DNA; AA; BIN; MORPH; CODON; NT2AA (default = DNA).
                    ]

    # Upload file. 
    print('\nUploading file, and initiating job\n')
    submission_info = str(subprocess.check_output(command_list, cwd=os.path.dirname(infile_path)), 'utf-8')

    # Check that it was submitted without error.
    assert '<error>' not in submission_info, """An error occured when
    submitting the job."""
        
    # Manage job submission.
    manage_any_cipres_job(infile_path, submission_info, timestamp, cra_user,
            password, key)



if __name__ == '__main__':

    # Define timestamp for output.
    timestamp = timestamp = time.strftime("%Y%m%d%H%M%S")
    
    # Get input alignment file path from command line.
    command_line_list = sys.argv
    infp = str(command_line_list[1])
    
    # Get login info.
    url = 'https://cipresrest.sdsc.edu/cipresrest/v1'
    password = getpass.getpass('Password for CIPRES: ')

    #################################################
    # Define these variables manually:
    cra_user = 'laelbarlow'
    application_name = 'archaeplastid_trees-1D575AEF608B49A1A8070479BADD492F'
    #################################################

    # Call function to run the tree on CIPRES. 
    submit_iqtree_job_to_cipres_1(url, cra_user, password, application_name, timestamp, infp)



