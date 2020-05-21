#!/usr/bin/env python

# Script to run GSEA from command line to run on large dataset or for multiple comparison: 
# Author Yogesh Dhungana
# Load java version java/openjdk-11 in HPCF

import subprocess
import os
import re
import logging
import sys
import math
import argparse
from itertools import permutations
from datetime import date
logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',
    stream=sys.stdout)
log = logging.getLogger(__name__)
def parse_arguments():
    parser = argparse.ArgumentParser(prog='Wrapper for GSEA.', description='')
    parser.add_argument('--gct',nargs=1,type=str , help='Gene expression matrix in .gct format.')
    parser.add_argument('--cls', nargs=1, type=str, help='Path to class file in .cls format.') ## GSEA for all possible comparisons will be run to save time and effort. 
    parser.add_argument('--gmt', nargs=1, type=str, help='Path to geneset file in .gmt format.')
    parser.add_argument('--chip', nargs=1, type=str, help='Path to annotation file in .chip format.')
    parser.add_argument('--projectname', type=str, help='Prefix for project with absolute path.')
    parser.add_argument('--nplots', default=50, type=int, help='Number of plots to generate.')
    parser.add_argument('--nperms', default=1000, type=int, help='Number of permutations.')
    parser.add_argument('--metric', type=str,default='Signal2Noise',choices=['Signal2Noise','tTest', 'log2_Ratio_of_Classes'], 
                        help='Gene ranking method to choose (If number of samples in each group is more than 3, choose "Signal2Noise") (If number of samples in each group is less than 3, choose either "tTest" or "log2_Ratio_of_Classes")')
    
    parser.add_argument('--ispreranked',action='store_false', help='Analysis is preranked GSEA.')
    parser.add_argument('--rnk',nargs='*',type=str,required=False, help='Preranked genelist in .rnk format.') ## this is to take into account if one has to run pre rank gsea for multiple comparision. 
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')

    args = parser.parse_args()
    log.setLevel(args.log_level)
    log.info(sys.argv)
    
    if not args.ispreranked and args.rnk==None: 
        raise ValueError("Missing .rnk file(s), please provide rnk file(s) and rerun the program.")
        sys.exit()
    
    return args
    

def mkdir_p(dirname):
    if dirname and not os.path.exists(dirname):
        os.makedirs(dirname)


def gsea_cls_parser(cls):
    """Extract class(phenotype) name from .cls file.

    :param cls: the a class list instance or .cls file which is identical to GSEA input .
    :return: phenotype name and a list of class vector.
    """
    all_possible_comp = []
    if isinstance(cls, list) :
        classes = cls
        sample_name= unique(classes)
    elif isinstance(cls, str) :
        with open(cls) as c:
            file = c.readlines()
        classes = file[2].strip('\n').split(" ")
        sample_name = file[1].lstrip("# ").strip('\n').strip('\r').split(" ")
    else:
        raise Exception('Error parsing sample name!')
        
    perm = permutations(sample_name, 2)
    for i in list(perm): all_possible_comp.append(i)
    
    return all_possible_comp




def run_shell_cmd(cmd):
    p = subprocess.Popen(
        ['/bin/bash', '-o', 'pipefail'],  # to catch error in pipe
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        preexec_fn=os.setsid)  # to make a new process with a new PGID
    pid = p.pid
    pgid = os.getpgid(pid)
    log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd))
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    err_str = 'PID={}, PGID={}, RC={}\nSTDERR={}\nSTDOUT={}'.format(
        pid, pgid, rc, stderr.strip(), stdout.strip())
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(err_str)
    else:
        log.info(err_str)
    return stdout.strip('\n')




def run_gsea(expression, class_, geneset, annotation, analysisName ,metric,nplots, nperm):
    """ run GSEA for microarray dataa. 
    
    :parm expression: gene expression matrix in .gct format.
    :parm class_: the class file in .cls format identical to GSEA input.
    :geneset: geneset to run GSEA against supplied as .gmt file. 
    :annotation: chip platform used to perfrom microarray in .chip format.
    :analysisName: Name of the analysis directory. 
    :metric: ranking method for gene it can be either signal2noise, tTest, log2FC (look up --help option)
    :nplots: number of plots to generate for top pathways by default 20 top pathways will be used to generate the plot.
    :nperm: number of permutations to perfrom by default 1000 permutations will be done. 
    """
    comparisonlist = gsea_cls_parser(class_)
    today = date.today()
    d1 = today.strftime("%Y_%m_%d") ## time stamp the analysis directory
    result_dir = analysisName+'_'+d1 
    comp = []
    prefix = []
    for i in range(len(comparisonlist)):
        comp.append("#"+comparisonlist[i][0]+"_versus_"+comparisonlist[i][1])
        prefix.append(comparisonlist[i][0]+"_versus_"+comparisonlist[i][1])
    for val in range(len(comp)):
        cmd1 = "/home/ydhungan/GSEA_4.0.3/gsea-cli.sh GSEA -res {} -cls {} -gmx {} -chip {} -collapse true -mode Max_probe -norm meandiv -nperm {} "
        cmd1 += "-permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label {} "
        cmd1 += "-metric {} -sort real -order descending -include_only_symbols true -make_sets true "
        cmd1 += "-median false -num 100 -plot_top_x {} -rnd_seed timestamp -save_rnd_lists false "
        cmd1 += "-set_max 500 -set_min 15 -zip_report false -out {}"
        cmd1 = cmd1.format(expression, class_+comp[val], geneset, annotation, nperm, prefix[val], metric, nplots, result_dir)
        run_shell_cmd(cmd1)


def strip_ext_rnk(rankfile):
    return re.sub(r'\.(rnk)$', '',
                  str(rankfile))

def run_prerank_gsea(ranklist, geneset, analysisName, nplots, nperm):
    """ run GSEA for preranked file or files 
    
    :parm ranklist: gene expression ranking in .rnk format can be a single file or list of files.
    :geneset: geneset to run GSEA against supplied as .gmt file. 
    :analysisName: Name of the project include the geneset name used. 
    :nplots: number of plots to generate for top pathways by default 20 top pathways will be used to generate the plot.
    :nperm: number of permutations to perfrom by default 1000 permutations will be done. 
    """
    today = date.today()
    d1 = today.strftime("%Y_%m_%d")
    if len(ranklist) == 1:
        head, tail = os.path.split(ranklist[0])
        prefix = strip_ext_rnk(tail)
        result_dir = head + '/' + prefix + '_' + analysisName+ '_' + d1 
        cmd1 = "/home/ydhungan/GSEA_4.0.3/gsea-cli.sh GSEAPreranked -rnk {} -gmx {} -nperm {} -scoring_scheme weighted " 
        cmd1 += "-set_max 500 -set_min 15 -plot_top_x {} -rnd_seed timestamp -zip_report false -rpt_label {} "
        cmd1 += "-out {} -norm meandiv "
        cmd1 = cmd1.format(ranklist[0], geneset, nperm, nplots, prefix, result_dir)
        run_shell_cmd(cmd1)
    else:
        for file in range(len(ranklist)):
            head, tail = os.path.split(ranklist[file])
            prefix = strip_ext_rnk(tail)
            result_dir = head + '/' + prefix + '_' + analysisName + '_' +  d1
            cmd1 = "/home/ydhungan/GSEA_4.0.3/gsea-cli.sh GSEAPreranked -rnk {} -gmx {} -nperm {} -scoring_scheme weighted " 
            cmd1 += "-set_max 500 -set_min 15 -plot_top_x {} -rnd_seed timestamp -zip_report false -rpt_label {} "
            cmd1 += "-out {} -norm meandiv "
            cmd1 = cmd1.format(ranklist[file], geneset, nperm, nplots, prefix, result_dir)
            run_shell_cmd(cmd1)
    




def main():
    log.info('Parsing arguments and making Project directories...')
    
    args = parse_arguments()
    print(args)
    #cmd = 'module load java/openjdk-11'   
    #os.system(cmd)
    if args.ispreranked:
        log.info('Running GSEA in normal mode...')
        run_gsea(args.gct[0], args.cls[0], args.gmt[0], args.chip[0], args.projectname, args.metric, args.nplots, args.nperms)
        
    else:
        log.info('Running GSEA in preranked mode...')
        run_prerank_gsea(args.rnk, args.gmt[0], args.projectname, args.nplots, args.nperms)
             
    log.info('All Done...')
if __name__ == '__main__':
    main()

