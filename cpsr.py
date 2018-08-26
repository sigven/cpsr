#!/usr/bin/env python

import csv
import re
import argparse
import os
import subprocess
import logging
import sys
import getpass
import platform
import toml

pcgr_version = '0.6.3'
cpsr_version = '0.1.0'
db_version = 'PCGR_DB_VERSION = 20180825'
vep_version = '93'
global vep_assembly


def __main__():
   
   parser = argparse.ArgumentParser(description='Cancer Predisposition Sequencing Report (CPSR) - report of cancer-predisposing germline variants',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('--input_vcf', dest = "input_vcf", help='VCF input file with somatic query variants (SNVs/InDels).')
   parser.add_argument('--force_overwrite', action = "store_true", help='By default, the script will fail with an error if any output file already exists. You can force the overwrite of existing result files by using this flag')
   parser.add_argument('--version', action='version', version='%(prog)s ' + str(cpsr_version))
   parser.add_argument('--basic',action="store_true",help="Run functional variant annotation on VCF through VEP/vcfanno, omit report generation (STEP 4)")
   parser.add_argument('--docker-uid', dest='docker_user_id', help='Docker user ID. Default is the host system user ID. If you are experiencing permission errors, try setting this up to root (`--docker-uid root`)')
   parser.add_argument('--no-docker', action='store_true', dest='no_docker', default=False, help='Run the CPSR workflow in a non-Docker mode (see install_no_docker/ folder for instructions')
   parser.add_argument('pcgr_base_dir',help='Directory that contains the PCGR data bundle directory, e.g. ~/pcgr-0.6.3')
   parser.add_argument('output_dir',help='Output directory')
   parser.add_argument('genome_assembly',choices = ['grch37','grch38'], help='Genome assembly build: grch37 or grch38')
   parser.add_argument('configuration_file',help='Configuration file (TOML format)')
   parser.add_argument('sample_id',help="Sample identifier - prefix for output files")
   
   docker_image_version = 'sigven/pcgr:' + str(pcgr_version)
   args = parser.parse_args()
   
   overwrite = 0
   if args.force_overwrite is True:
      overwrite = 1
   
   if args.no_docker:
      docker_image_version = None
   else:
      # check that script and Docker image version correspond
      check_docker_command = 'docker images -q ' + str(docker_image_version)
      output = subprocess.check_output(str(check_docker_command), stderr=subprocess.STDOUT, shell=True)
      logger = getlogger('cpsr-validate-config')
      if(len(output) == 0):
          err_msg = 'Docker image ' + str(docker_image_version) + ' does not exist, pull image from Dockerhub (docker pull ' + str(docker_image_version) + ')'
          error_message(err_msg,logger)
   
   config_options = {}
   if os.path.exists(args.configuration_file):
      config_options = read_config_options(args.configuration_file, args.pcgr_base_dir, args.genome_assembly, logger)
   else:
      err_msg = "PCGR configuration file " + str(args.configuration_file) + " does not exist - exiting"
      error_message(err_msg,logger)

   logger = getlogger('pcgr-check-files')
   host_directories = verify_input_files(args.input_vcf, args.configuration_file, config_options, args.pcgr_base_dir, args.output_dir, args.sample_id, args.genome_assembly, overwrite, logger)

   run_cpsr(host_directories, docker_image_version, config_options, args.sample_id, args.genome_assembly, cpsr_version, args.basic, docker_user_id=args.docker_user_id)


def read_config_options(configuration_file, pcgr_dir, genome_assembly, logger):
   
   ## read default options
   cpsr_config_options = {}
   cpsr_config_file_default = os.path.join(pcgr_dir,'data',str(genome_assembly),'cpsr_configuration_default.toml')
   if not os.path.exists(cpsr_config_file_default):
      err_msg = "Default CPSR configuration file " + str(cpsr_config_file_default) + " does not exist - exiting"
      error_message(err_msg,logger)
   try:
      cpsr_config_options = toml.load(cpsr_config_file_default)
   except (IndexError,TypeError):
      err_msg = 'Configuration file ' + str(configuration_file) + ' is not formatted correctly'
      error_message(err_msg, logger)

   ## override with options set by the users
   try:
      user_options = toml.load(configuration_file)
   except (IndexError,TypeError):
      err_msg = 'Configuration file ' + str(configuration_file) + ' is not formatted correctly'
      error_message(err_msg, logger)

   for section in cpsr_config_options:
      if section in user_options:
         for var in cpsr_config_options[section]:
            if not var in user_options[section]:
               continue
            if isinstance(cpsr_config_options[section][var],bool) and not isinstance(user_options[section][var],bool):
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + ' cannot be parsed properly (expecting boolean)'
               error_message(err_msg, logger)
            if isinstance(cpsr_config_options[section][var],int) and not isinstance(user_options[section][var],int):
               err_msg = 'Configuration value \"' + str(user_options[section][var]) + '\" for ' + str(var) + ' cannot be parsed properly (expecting integer)'
               error_message(err_msg, logger)
            if isinstance(cpsr_config_options[section][var],float) and (not isinstance(user_options[section][var],float) and not isinstance(user_options[section][var],int)):
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + ' cannot be parsed properly (expecting float)'
               error_message(err_msg, logger)
            if isinstance(cpsr_config_options[section][var],str) and not isinstance(user_options[section][var],str):
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + ' cannot be parsed properly (expecting string)'
               error_message(err_msg, logger)
            theme_options = ['default', 'cerulean', 'journal', 'flatly', 'readable', 'spacelab', 'united', 'cosmo', 'lumen', 'paper', 'sandstone', 'simplex','yeti']
            if var == 'report_theme' and not str(user_options[section][var]) in theme_options:
               err_msg = 'Configuration value ' + str(user_options[section][var]) + ' for ' + str(var) + ' cannot be parsed properly (expecting \'default\', \'cerulean\', \'journal\', \'flatly\', \'readable\', \'spacelab\', \'united\', \'cosmo\', \'lumen\', \'paper\', \'sandstone\', \'simplex\',or \'yeti\')'
               error_message(err_msg, logger)
           
   return cpsr_config_options


def error_message(message, logger):
   logger.error('')
   logger.error(message)
   logger.error('')
   exit(0)

def pcgr_warn_message(message, logger):
   logger.warning('')
   logger.warning(message)
   logger.warning('')

def verify_input_files(input_vcf, configuration_file, cpsr_config_options, base_pcgr_dir, output_dir, sample_id, genome_assembly, overwrite, logger):
   """
   Function that checks the input files and directories provided by the user and checks for their existence
   """
 
   input_vcf_dir = "NA"
   input_conf_dir = "NA"
   db_dir = "NA"
   base_dir = "NA"
   output_dir_full = "NA"
   input_vcf_basename = "NA"
   input_conf_basename = "NA"
   
   ## check that either input vcf or cna segments exist
   if input_vcf is None:
      err_msg = "Please specifiy a VCF input file (--input_vcf)"
      error_message(err_msg,logger)
   
   ## check the existence of given output folder
   output_dir_full = os.path.abspath(output_dir)
   if not os.path.isdir(output_dir_full):
      err_msg = "Output directory (" + str(output_dir_full) + ") does not exist"
      error_message(err_msg,logger)
   
   ## check if input vcf exist
   if not input_vcf is None:
      if not os.path.exists(os.path.abspath(input_vcf)):
         err_msg = "Input file (" + str(input_vcf) + ") does not exist"
         error_message(err_msg,logger)

      if not (os.path.abspath(input_vcf).endswith('.vcf') or os.path.abspath(input_vcf).endswith('.vcf.gz')):
         err_msg = "VCF input file (" + os.path.abspath(input_vcf) + ") does not have the correct file extension (.vcf or .vcf.gz)"
         error_message(err_msg,logger)

      ## check that tabix file exist if bgzipped files is given
      if os.path.abspath(input_vcf).endswith('.vcf.gz'):
         tabix_file = input_vcf + '.tbi'
         if not os.path.exists(os.path.abspath(tabix_file)):
            err_msg = "Tabix file (i.e. '.gz.tbi') is not present for the bgzipped VCF input file (" + os.path.abspath(input_vcf) + "). Please make sure your input VCF is properly compressed and indexed (bgzip + tabix)"
            error_message(err_msg,logger)

      input_vcf_basename = os.path.basename(str(input_vcf))
      input_vcf_dir = os.path.dirname(os.path.abspath(input_vcf))

      ## if output vcf exist and overwrite not set
      output_vcf = os.path.join(str(output_dir_full),str(sample_id)) + '.cpsr.' + str(genome_assembly) + '.vcf.gz'
      if os.path.exists(output_vcf) and overwrite == 0:
         err_msg = "Output files (e.g. " + str(output_vcf) + ") already exist - please specify different sample_id or add option --force_overwrite"
         error_message(err_msg,logger)
   
   if not configuration_file is None:
      if not os.path.exists(os.path.abspath(configuration_file)):
         err_msg = "Input file (" + str(configuration_file) + ") does not exist"
         error_message(err_msg,logger)

      if not os.path.abspath(configuration_file).endswith('.toml'):
         err_msg = "Configuration file (" + os.path.abspath(configuration_file) + ") does not have the correct file extension (.toml)"
         error_message(err_msg,logger)

      input_conf_basename = os.path.basename(str(configuration_file))
      input_conf_dir = os.path.dirname(os.path.abspath(configuration_file))
   
   ## check the existence of base folder
   base_dir = os.path.abspath(base_pcgr_dir)
   if not os.path.isdir(base_dir):
      err_msg = "Base directory (" + str(base_dir) + ") does not exist"
      error_message(err_msg,logger)
   
   ## check the existence of data folder within the base folder
   db_dir = os.path.join(os.path.abspath(base_pcgr_dir),'data')
   if not os.path.isdir(db_dir):
      err_msg = "Data directory (" + str(db_dir) + ") does not exist"
      error_message(err_msg,logger)
   
   ## check the existence of specified assembly data folder within the base folder
   db_assembly_dir = os.path.join(os.path.abspath(base_pcgr_dir),'data',genome_assembly)
   if not os.path.isdir(db_assembly_dir):
      err_msg = "Data directory for the specified genome assembly (" + str(db_assembly_dir) + ") does not exist"
      error_message(err_msg,logger)
   
   ## check the existence of RELEASE_NOTES
   rel_notes_file = os.path.join(os.path.abspath(base_pcgr_dir),'data',genome_assembly,'RELEASE_NOTES')
   if not os.path.exists(rel_notes_file):
      err_msg = 'The PCGR data bundle is outdated - please download the latest data bundle (see github.com/sigven/cpsr for instructions)'
      error_message(err_msg,logger)
      
   f_rel_not = open(rel_notes_file,'r')
   compliant_data_bundle = 0
   for line in f_rel_not:
      if db_version in line:
         compliant_data_bundle = 1
   
   f_rel_not.close()
    
   if compliant_data_bundle == 0:
      err_msg = 'The PCGR data bundle is not compliant with the software version - please download the latest software and data bundle (see https://github.com/sigven/cpsr for instructions)'
      error_message(err_msg,logger)
   
   host_directories = {}
   host_directories['input_vcf_dir_host'] = input_vcf_dir
   host_directories['input_conf_dir_host'] = input_conf_dir
   host_directories['db_dir_host'] = db_assembly_dir
   host_directories['base_dir_host'] = base_dir
   host_directories['output_dir_host'] = output_dir_full
   host_directories['input_vcf_basename_host'] = input_vcf_basename
   host_directories['input_conf_basename_host'] = input_conf_basename

   return host_directories
   

def check_subprocess(command):
   try:
      output = subprocess.check_output(str(command), stderr=subprocess.STDOUT, shell=True)
      if len(output) > 0:
         print (str(output.decode()).rstrip())
   except subprocess.CalledProcessError as e:
      print (e.output.decode())
      exit(0)

def getlogger(logger_name):
   logger = logging.getLogger(logger_name)
   logger.setLevel(logging.DEBUG)

   # create console handler and set level to debug
   ch = logging.StreamHandler(sys.stdout)
   ch.setLevel(logging.DEBUG)

   # add ch to logger
   logger.addHandler(ch)
   
   # create formatter
   formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s", "20%y-%m-%d %H:%M:%S")
   
   #add formatter to ch
   ch.setFormatter(formatter)
   
   return logger

def run_cpsr(host_directories, docker_image_version, config_options, sample_id, genome_assembly, cpsr_version, basic, docker_user_id=None):
   """
   Main function to run the cpsr workflow using Docker
   """
   
   ## set basic Docker run commands
   output_vcf = 'None'
   output_pass_vcf = 'None'
   output_pass_tsv = 'None'
   uid = ''
   gencode_version = 'release 27'
   vep_assembly = 'GRCh38'
   if genome_assembly == 'grch37':
      gencode_version = 'release 19'
      vep_assembly = 'GRCh37'

   logger = getlogger('cpsr-get-OS')

   if docker_user_id:
      uid = docker_user_id
   elif platform.system() == 'Linux' or platform.system() == 'Darwin' or sys.platform == 'darwin' or sys.platform == 'linux2' or sys.platform == 'linux':
      uid = os.getuid()
   else:
      if platform.system() == 'Windows' or sys.platform == 'win32' or sys.platform == 'cygwin':
         uid = getpass.getuser()

   if uid == '':
      logger.warning('Was not able to get user id/username for logged-in user on the underlying platform (platform.system(): ' + str(platform.system()) + ', sys.platform: ' + str(sys.platform) + '), now running CPSR as root')
      uid = 'root'
   
   vepdb_dir_host = os.path.join(str(host_directories['db_dir_host']),'.vep')

   input_vcf_docker = 'None'
   input_conf_docker = 'None'
   
   if docker_image_version:
      if host_directories['input_vcf_basename_host'] != 'NA':
         input_vcf_docker = '/workdir/input_vcf/' + str(host_directories['input_vcf_basename_host'])
      if host_directories['input_conf_basename_host'] != 'NA':
         input_conf_docker = '/workdir/input_conf/' + str(host_directories['input_conf_basename_host'])

      docker_command_run1 = 'NA'
      if host_directories['input_vcf_dir_host'] != 'NA':
         docker_command_run1 = "docker run --rm -t -u " + str(uid) + " -v=" + str(host_directories['base_dir_host']) + ":/data -v=" + str(vepdb_dir_host) + ":/usr/local/share/vep/data -v=" + str(host_directories['input_vcf_dir_host']) + ":/workdir/input_vcf -v=" + str(host_directories['output_dir_host']) + ":/workdir/output -w=/workdir/output " + str(docker_image_version) + " sh -c \""      
      if host_directories['input_conf_dir_host'] != 'NA':
         docker_command_run1 = "docker run --rm -t -u " + str(uid) + " -v=" + str(host_directories['base_dir_host']) + ":/data -v=" + str(vepdb_dir_host) + ":/usr/local/share/vep/data -v=" + str(host_directories['input_vcf_dir_host']) + ":/workdir/input_vcf -v=" + str(host_directories['input_conf_dir_host']) + ":/workdir/input_conf -v=" + str(host_directories['output_dir_host']) + ":/workdir/output -w=/workdir/output " + str(docker_image_version) + " sh -c \""
      docker_command_run2 = "docker run --rm -t -u " + str(uid) + " -v=" + str(host_directories['base_dir_host']) + ":/data -v=" + str(host_directories['output_dir_host']) + ":/workdir/output -w=/workdir " + str(docker_image_version) + " sh -c \""
      docker_command_run_end = '\"'

      data_dir = '/data'
      output_dir = '/workdir/output'
      vep_dir = '/usr/local/share/vep/data'
      r_scripts_dir = '/'
      python_scripts_dir = ''

   else:
      if host_directories['input_vcf_basename_host'] != 'NA':
         input_vcf_docker = os.path.join(host_directories['input_vcf_dir_host'], host_directories['input_vcf_basename_host'])
      if host_directories['input_conf_basename_host'] != 'NA':
         input_conf_docker = os.path.join(host_directories['input_conf_dir_host'], host_directories['input_conf_basename_host'])

      docker_command_run1 = ''
      docker_command_run2 = ''
      docker_command_run_end = ''

      data_dir = host_directories['base_dir_host']
      output_dir = host_directories['output_dir_host']
      vep_dir = vepdb_dir_host

      src_dir = os.path.join(host_directories['base_dir_host'], 'src')
      python_scripts_dir = os.path.join(src_dir, 'pcgr')
      os.environ['PYTHONPATH'] = os.path.join(python_scripts_dir, 'lib') + ((':' + os.environ['PYTHONPATH']) if 'PYTHONPATH' in os.environ else '')
      r_scripts_dir = src_dir

   check_subprocess(docker_command_run1.replace("-u " + str(uid), "") + 'mkdir -p ' + output_dir + docker_command_run_end)
   
   ## verify VCF and CNA segment file
   logger = getlogger('cpsr-validate-input')
   logger.info("STEP 0: Validate input data")
   vcf_validate_command = str(docker_command_run1) + "cpsr_validate_input.py " + str(data_dir) + " " + str(input_vcf_docker) + " " + str(input_conf_docker) + " " + str(genome_assembly) + docker_command_run_end
   check_subprocess(vcf_validate_command)
   ## Log tumor type of query genome
   logger.info('Finished')
   
   if not input_vcf_docker == 'None':
      
      ## Define input, output and temporary file names
      pcgr_model = 'cpsr'
      output_vcf = os.path.join(output_dir, str(sample_id) + '.' + str(pcgr_model) + '.' + str(genome_assembly) + '.vcf.gz')
      output_pass_vcf = os.path.join(output_dir, str(sample_id) + '.' + str(pcgr_model) + '.' + str(genome_assembly) + '.pass.vcf.gz')
      output_pass_tsv = os.path.join(output_dir, str(sample_id) + '.' + str(pcgr_model) + '.' + str(genome_assembly) + '.pass.tsv')
      input_vcf_cpsr_ready = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_ready.vcf.gz',host_directories['input_vcf_basename_host']))
      input_vcf_cpsr_ready_uncompressed = os.path.join(output_dir, re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_ready.vcf',host_directories['input_vcf_basename_host']))
      vep_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_vep.vcf',input_vcf_cpsr_ready)
      vep_vcfanno_vcf = re.sub(r'(\.vcf$|\.vcf\.gz$)','.cpsr_vep.vcfanno.vcf',input_vcf_cpsr_ready)
      vep_tmp_vcf = vep_vcf + '.tmp'
      vep_vcfanno_annotated_vcf = re.sub(r'\.vcfanno','.vcfanno.annotated',vep_vcfanno_vcf) + '.gz'
      vep_vcfanno_annotated_pass_vcf = re.sub(r'\.vcfanno','.vcfanno.annotated.pass',vep_vcfanno_vcf) + '.gz'

      fasta_assembly = os.path.join(vep_dir, "homo_sapiens", str(vep_version) + "_" + str(vep_assembly), "Homo_sapiens." + str(vep_assembly) + ".dna.primary_assembly.fa.gz")
      vep_options = "--vcf --check_ref --flag_pick_allele --force_overwrite --species homo_sapiens --assembly " + str(vep_assembly) + " --offline --fork " + str(config_options['other']['n_vep_forks']) + " --hgvs --dont_skip --failed 1 --af --af_1kg --af_gnomad --variant_class --regulatory --domains --symbol --protein --ccds --uniprot --appris --biotype --canonical --gencode_basic --cache --numbers --total_length --no_stats --allele_number --no_escape --xref_refseq --plugin LoF --dir " + str(vep_dir)
      if config_options['other']['vep_skip_intergenic'] == 1:
         vep_options = vep_options + " --no_intergenic"
      vep_main_command = str(docker_command_run1) + "vep --input_file " + str(input_vcf_cpsr_ready) + " --output_file " + str(vep_tmp_vcf) + " " + str(vep_options) + " --fasta " + str(fasta_assembly) + docker_command_run_end
      vep_sed_command =  str(docker_command_run1) + "sed -r 's/:p\.[A-Z]{1}[a-z]{2}[0-9]+=//g' " + str(vep_tmp_vcf) + " > " + str(vep_vcf) + docker_command_run_end
      vep_bgzip_command = str(docker_command_run1) + "bgzip -f " + str(vep_vcf) + docker_command_run_end
      vep_tabix_command = str(docker_command_run1) + "tabix -f -p vcf " + str(vep_vcf) + ".gz" + docker_command_run_end
      logger = getlogger('cpsr-vep')

      print()
      logger.info("STEP 1: Basic variant annotation with Variant Effect Predictor (" + str(vep_version) + ", GENCODE " + str(gencode_version) + ", " + str(genome_assembly) + ")")
      #print(str(vep_main_command))
      #return

      check_subprocess(vep_main_command)
      check_subprocess(vep_sed_command)
      check_subprocess(vep_bgzip_command)
      check_subprocess(vep_tabix_command)
      logger.info("Finished")
   
      ## vcfanno command
      print()
      logger = getlogger('cpsr-vcfanno')
      logger.info("STEP 2: Annotation for cancer predisposition with pcgr-vcfanno (ClinVar, dbNSFP, UniProtKB, CiVIC)")
      pcgr_vcfanno_command = str(docker_command_run2) + os.path.join(python_scripts_dir, "pcgr_vcfanno.py") + " --num_processes "  + str(config_options['other']['n_vcfanno_proc']) + " --dbnsfp --clinvar --civic --uniprot --pcgr_onco_xref " + str(vep_vcf) + ".gz " + str(vep_vcfanno_vcf) + " " + os.path.join(data_dir, "data", str(genome_assembly)) + docker_command_run_end
      check_subprocess(pcgr_vcfanno_command)
      logger.info("Finished")
   
      ## summarise command
      print()
      logger = getlogger("cpsr-summarise")
      pcgr_summarise_command = str(docker_command_run2) + os.path.join(python_scripts_dir, "pcgr_summarise.py") + " " + str(vep_vcfanno_vcf) + ".gz " + os.path.join(data_dir, "data", str(genome_assembly)) + " --cpsr" + docker_command_run_end
      logger.info("STEP 3: Cancer gene annotations with cpsr-summarise")
      check_subprocess(pcgr_summarise_command)
      
      create_output_vcf_command1 = str(docker_command_run2) + 'mv ' + str(vep_vcfanno_annotated_vcf) + ' ' + str(output_vcf) + "\""
      create_output_vcf_command2 = str(docker_command_run2) + 'mv ' + str(vep_vcfanno_annotated_vcf) + '.tbi ' + str(output_vcf) + '.tbi' + "\""
      create_output_vcf_command3 = str(docker_command_run2) + 'mv ' + str(vep_vcfanno_annotated_pass_vcf) + ' ' + str(output_pass_vcf) + "\""
      create_output_vcf_command4 = str(docker_command_run2) + 'mv ' + str(vep_vcfanno_annotated_pass_vcf) + '.tbi ' + str(output_pass_vcf) + '.tbi' + "\""
      clean_command = str(docker_command_run2) + 'rm -f ' + str(vep_vcf) + '* ' + str(vep_vcfanno_annotated_vcf) + ' ' + str(vep_vcfanno_annotated_pass_vcf) + '* ' + str(vep_vcfanno_vcf) + '* ' +  str(input_vcf_cpsr_ready_uncompressed) + "* "  + docker_command_run_end
      check_subprocess(create_output_vcf_command1)
      check_subprocess(create_output_vcf_command2)
      check_subprocess(create_output_vcf_command3)
      check_subprocess(create_output_vcf_command4)
      cpsr_vcf2tsv_command = str(docker_command_run2) + os.path.join(python_scripts_dir, "vcf2tsv.py") + " " + str(output_pass_vcf) + " --compress " + str(output_pass_tsv) + docker_command_run_end
      logger.info("Converting VCF to TSV with https://github.com/sigven/vcf2tsv")
      check_subprocess(cpsr_vcf2tsv_command)
      check_subprocess(clean_command)
      logger.info("Finished")

      #return
  
   print()
   
   ## Generation of HTML reports for VEP/vcfanno-annotated VCF file
   if not basic: 
      logger = getlogger('cpsr-writer')
      logger.info("STEP 4: Generation of output files - Cancer predisposition sequencing report")
      cpsr_report_command = (docker_command_run1 + os.path.join(r_scripts_dir, "cpsr.R") + " " + output_dir + " " + str(output_pass_tsv) + ".gz " +  str(sample_id)  + " " + str(input_conf_docker) + " " + str(cpsr_version) + " " + str(genome_assembly) + " " + os.path.join(data_dir, docker_command_run_end))
      #print(pcgr_report_command)
      check_subprocess(cpsr_report_command)
      logger.info("Finished")
   
   
   
if __name__=="__main__": __main__()

