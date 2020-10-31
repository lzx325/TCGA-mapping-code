#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import subprocess as sp
import pandas as pd
import os
import shlex
import pprint
import time
import glob
import sys
WORK_DIR="/ibex/scratch/projects/c2066/"
file_pattern={
    'step2':"02.downTCGA/*",
    'step3':"03.STAR_1stAlign/TCGA-*.SJ.out.tab",
    'step5':"05.STAR_2ndAlign/TCGA-*.rmdup.bam",
    'step6':"06.buildWhippletIndex/*.jls",
    'step7':"07.WhippletQuant/*.psi.gz"
}
class slurmjob(object):
    def __init__(self,cmd,job_name="TCGA",log_dir="slurm",memory=100,time=1,n_cpus=8):
        self.cmd=cmd
        self.job_name=job_name
        self.log_dir=log_dir
        self.memory=memory
        if not os.path.isdir(self.log_dir):
            os.makedirs(self.log_dir)
        self.__jobid=None
        self.__exitcode=None
        self.__jobstate=None
        self.time=time
        self.n_cpus=n_cpus
    def run(self):
        if self.__jobid != None:
            raise RuntimeError("Job is already run")
        slurm_cmd='''\
sbatch <<- EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=batch
#SBATCH -J {job_name}
#SBATCH -o slurm/{job_name}.%J.out
#SBATCH -e slurm/{job_name}.%J.err
#SBATCH --time={time}-00:00:00
#SBATCH --mem={memory}G
#SBATCH --cpus-per-task={n_cpus}
#run the application:
{cmd}
EOF
'''.format(job_name=self.job_name,cmd=self.cmd,memory=self.memory,time=self.time,n_cpus=self.n_cpus)
        print("command:")
        print(slurm_cmd)
        print("="*20)
        proc=sp.Popen(slurm_cmd,shell=True,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)
        stdout,stderr=proc.communicate()
        if proc.returncode==0:
            stdout=stdout.rstrip()
            slurm_jobid=int(stdout.split(' ')[-1])
            self.__jobid=slurm_jobid
            print("submitted job %d"%(slurm_jobid))
        else:
            print(stderr)
            raise RuntimeError("sbatch exited with status %d"%(proc.returncode))
    def update(self):
        if self.__jobid ==None:
            raise RuntimeError("Job is not run")
        # must sleep, will fail otherwise
        time.sleep(10)
        sacct_cmd="sacct -X -o state,exitcode -n -j {jobid}".format(jobid=self.__jobid)
        proc=sp.Popen(shlex.split(sacct_cmd),stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)
        stdout,stderr=proc.communicate()
        if proc.returncode==0:
            stdout=stdout.strip()
            stdout_split=stdout.split()
            if len(stdout_split)!=2:
                print("sacct did not return correct format. sacct gives result %r"%(stdout_split),file=sys.stderr)
                return
            state,exitcode=stdout.split()
            exitcode=tuple(int(c) for c in exitcode.split(':'))
            self.__jobstate=state
            self.__exitcode=exitcode
        else:
            print(stderr)
            raise RuntimeError("sacct exited with status %d"%(proc.returncode))

    def get_status(self):
        self.update()
        return (self.__jobstate,self.__exitcode)

def direct_run(cmd,shell=False):
    if not shell:
        cmd=shlex.split(cmd)
    proc=sp.Popen(cmd,shell=shell,stdout=sp.PIPE,stderr=sp.PIPE,universal_newlines=True)
    stdout,stderr=proc.communicate()
    print("="*20)
    print("%s finished")
    print("="*20)
    print("stdout:")
    print(stdout)
    print("="*20)
    print("stderr:")
    print(stderr)
    print("="*20)

def check_sample_list(cancer,step,check_download=True,skip=True,verbose=True):
    assert step in ["step3","step5","step6","step7"]
    manifest_fn="TCGA_all_sample_info.txt"
    manifest=pd.read_csv(os.path.join(WORK_DIR,"anz/data/manifest",manifest_fn),sep="\t")
    manifest_cancer=manifest.groupby("cases").get_group(cancer)
    manifest_file_id_list=manifest_cancer["file_id"].values.tolist()
    if check_download:
        download_dir=os.path.join(WORK_DIR,"anz/result",cancer,"02.downTCGA")
        downloaded_set=set(os.listdir(download_dir))
        for f in manifest_file_id_list:
            if f not in downloaded_set:
                raise FileNotFoundError("%s not found in download"%(f))
    manifest_barcodes=manifest_cancer["associated_entities"]
    if skip:
        processed_barcodes=get_processed_barcodes(cancer,step)
        indices = [i+1 for i,barcode in enumerate(manifest_barcodes) if barcode not in processed_barcodes]
        processed_indices=[i+1 for i,barcode in enumerate(manifest_barcodes) if barcode in processed_barcodes]
        if verbose:
            print("unprocessed indices:",indices)
            print("processed indices:",processed_indices)
        return indices,processed_indices
    else:
        indices = list(range(1,manifest_cancer.shape[0]+1))
        return indices


def main2():
    DEVNULL=open(os.devnull,"w")
    proc=sp.Popen(["bash","test.sh"],env={"myvar":"3"},stdout=sp.PIPE,stderr=DEVNULL,universal_newlines=True)
    stdout,stderr=proc.communicate()
    print("from main.py",type(stdout))
    print("from main.py",stderr)
    print("retcode",proc.returncode)

def step2(cancer_type):
    script_dir="02.downTCGA"
    script_file="02.downTCGA.sh"
    script_path=os.path.join(script_dir,script_file)
    os.system("bash {script_path} {work_dir} {cancer_type}".format(
        script_path=script_path,
        work_dir=WORK_DIR,
        cancer_type=cancer_type
    ))
def step3(cancer_type):
    script_dir="03.STAR_1stAlign"
    script_file="03.STAR_1stAlign.sh"
    script_path=os.path.join(script_dir,script_file)
    unprocessed_indices,processed_indices=check_sample_list(cancer_type,"step3")
    thread=30
    print("#samples in %s:"%(cancer_type),len(unprocessed_indices))
    for index in unprocessed_indices:
        cmd="bash {script_path} {work_dir} {cancer_type} {index} {thread}".format(
            script_path=script_path,
            work_dir=WORK_DIR,
            cancer_type=cancer_type,
            index=index,
            thread=thread
        )
        job=slurmjob(cmd,cancer_type,"slurm")
        job.run()
        try:
            print("%d: %r"%(index,job.get_status()))
        except RuntimeError as e:
            print(e)

def step4(cancer_type):
    script_dir="04.STAR_2ndIndex"
    script_file="04.STAR_2ndIndex.sh"
    script_path=os.path.join(script_dir,script_file)
    thread=30
    memory=300
    cmd="bash {script_path} {work_dir} {cancer_type} {thread}".format(
            script_path=script_path,
            work_dir=WORK_DIR,
            cancer_type=cancer_type,
            thread=thread

        )
    job=slurmjob(cmd,cancer_type,"slurm",memory)
    job.run()
    print(job.get_status())
def step5(cancer_type,time,n_cpus):
    script_dir="05.STAR_2ndAlign"
    script_file="05.STAR_2ndAlign.sh"
    script_path=os.path.join(script_dir,script_file)
    thread=n_cpus
    unprocessed_indices,processed_indices=check_sample_list(cancer_type,"step5")
    for index in unprocessed_indices:
        cmd="bash {script_path} {work_dir} {cancer_type} {index} {thread}".format(
                script_path=script_path,
                work_dir=WORK_DIR,
                cancer_type=cancer_type,
                index=index,
                thread=thread
            )
        job=slurmjob(cmd,cancer_type,"slurm",time=time,n_cpus=n_cpus)
        job.run()
        print(job.get_status())
def step6(cancer_type,time,n_cpus):
    script_dir="06.buildWhippletIndex"
    script_file="06.buildWhippletIndex.sh"
    script_path=os.path.join(script_dir,script_file)
    whipplet_src="/home/liz0f/.julia/v0.6/Whippet/bin/whippet-index.jl"
    memory=300
    unprocessed_indices,_=check_sample_list(cancer_type,"step6")
    for index in unprocessed_indices:
        cmd="bash {script_path} {work_dir} {cancer_type} {whipplet_src} {index}".format(
                script_path=script_path,
                work_dir=WORK_DIR,
                cancer_type=cancer_type,
                whipplet_src=whipplet_src,
                index=index
            )
        job=slurmjob(cmd,cancer_type,"slurm",memory=memory,time=time,n_cpus=n_cpus)
        job.run()
        print(job.get_status())
def step7(cancer_type,time,n_cpus):
    script_dir="07.WhippletQuant"
    script_file="07.WhippletQuant.sh"
    script_path=os.path.join(script_dir,script_file)
    whipplet_src="/home/liz0f/.julia/v0.6/Whippet/bin/whippet-quant.jl"
    memory=100
    unprocessed_indices,processed_indices=check_sample_list(cancer_type,"step7")
    for index in unprocessed_indices:
        cmd="bash {script_path} {work_dir} {cancer_type} {whipplet_src} {index}".format(
                script_path=script_path,
                work_dir=WORK_DIR,
                cancer_type=cancer_type,
                whipplet_src=whipplet_src,
                index=index
            )
        job=slurmjob(cmd,cancer_type,"slurm",memory=memory,time=time,n_cpus=n_cpus)
        job.run()
        print(job.get_status())
def sambaba():
    input_dir="/scratch/dragon/intel/liz0f/KEEPME/anz/result/TCGA-CHOL/05.STAR_2ndAlign"
    infn="star2pass_merged.sorted.bam"
    outfn="star2pass_merged.sorted.rmdup.bam"
    infp=os.path.join(input_dir,infn)
    outfp=os.path.join(input_dir,outfn)
    cmd = "sambamba markdup --hash-table-size 67108864 -r -t 50 -p %s %s"%(infp,outfp)
    memory=300
    job=slurmjob(cmd,"TCGA-CHOL","slurm",memory=memory)
    job.run()
    print(job.get_status())
def glob_name(cancer_types):
    base_dir=os.path.join(WORK_DIR,"anz/result")
    for i,cancer_type in enumerate(cancer_types):
        n_step2=len(glob.glob(os.path.join(base_dir,cancer_type,file_pattern["step2"])))
        n_step3=len(glob.glob(os.path.join(base_dir,cancer_type,file_pattern["step3"])))
        n_step5=len(glob.glob(os.path.join(base_dir,cancer_type,file_pattern["step5"])))
        n_step6=len(glob.glob(os.path.join(base_dir,cancer_type,file_pattern["step6"])))
        n_step7=len(glob.glob(os.path.join(base_dir,cancer_type,file_pattern["step7"])))
        print("Cancer type:",cancer_type,"Index:", i)
        print("n_step2:",n_step2)
        print("n_step3:",n_step3)
        print("n_step5:",n_step5)
        print("n_step6:",n_step6)
        print("n_step7:",n_step7)

def get_processed_barcodes(cancer_type,step):
    base_dir=os.path.join(WORK_DIR,"anz/result")
    assert step in file_pattern
    fps=glob.glob(os.path.join(base_dir,cancer_type,file_pattern[step]))
    barcodes=[os.path.basename(fp).split('.')[0] for fp in fps]
    return barcodes

def get_barcodes(cancer):
    manifest_fn="TCGA_all_sample_info.txt"
    manifest=pd.read_csv(os.path.join(WORK_DIR,"anz/data/manifest",manifest_fn),sep="\t")
    manifest_cancer=manifest.groupby("cases").get_group(cancer)
    return manifest_cancer["associated_entities"]

def ls_dir(cancer_types):
    base_dir=os.path.join(WORK_DIR,"anz/result")
    for cancer_type in cancer_types:
        print(os.listdir(os.path.join(base_dir,cancer_type)))

def write_failed_list(cancer_types,step="step7",out_file="failed_list.csv"):
    data_dict={
        'cancer_type':list(),
        'barcode':list(),
    }
    for cancer_type in cancer_types:
        unprocessed,_=check_sample_list(cancer_type,step,verbose=False)
        ae=get_barcodes(cancer_type)
        data_dict["cancer_type"]+=[cancer_type]*len(unprocessed)
        data_dict["barcode"]+=ae.iloc[[i-1 for i in unprocessed]].tolist()
    df=pd.DataFrame(data_dict)
    df.to_csv(out_file,header=True,index=False)
if __name__=="__main__":

    batch1=[
        "TCGA-CHOL",
        "TCGA-DLBC",
        "TCGA-UCS",
        "TCGA-ACC",
        "TCGA-UVM",
        "TCGA-MESO",
        "TCGA-KICH",
        "TCGA-THYM",
        "TCGA-READ",
        "TCGA-TGCT",
        "TCGA-LAML",
    ]

    batch2=["TCGA-PAAD","TCGA-PCPG","TCGA-GBM","TCGA-SARC"]
    batch2_remain=["TCGA-SARC"]
    batch3=["TCGA-CESC","TCGA-UCEC","TCGA-KIRP","TCGA-COAD"]
    batch4=["TCGA-BLCA","TCGA-LIHC","TCGA-SKCM","TCGA-ESCA","TCGA-PRAD","TCGA-LUAD"]
    if len(sys.argv)>=2 and sys.argv[1]=="run":
        for cancer in batch3[0:2]:
            # step2(cancer)
            # step6(cancer,time=3,n_cpus=16)
            # step5(cancer,time=1,n_cpus=16)
            # step3(cancer)
            step4(cancer)
            # step7(cancer, time=5, n_cpus=8)

    else:
        glob_name(batch3[0:2])
        

    # glob_name(batch3)


    # unprocessed,processed=check_sample_list("TCGA-READ","step7",verbose=False)
    # print(processed)
    # print(unprocessed)

    #step7("TCGA-READ",time=5,n_cpus=8)

    # write_failed_list(step_05[1:])



