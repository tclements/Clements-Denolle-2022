import os 
import glob
import shlex
import subprocess
import numpy as np
import pandas as pd
import re

def call_cvsmi(lon, lat, min_depth, max_depth, model, OUT, step=100):
    """Depth profile for California Community Velocity Model"""

    if not os.path.isdir(OUT):
        os.makedirs(OUT)
    locs = os.path.join(OUT,'locs')
    disp = os.path.join(OUT,'disp')

    # create points
    with open(locs, 'w') as f:
        for ii in np.arange(min_depth, max_depth, step):
            point = '{} {} {}\n'.format(lon, lat, ii)
            f.write(point)

    # run cvmsi_query
    # write bash script
    bash = os.path.join(OUT,'run_cvmsi.sh')
    to_call = 'cvmsi_query -m {} < {} > {}'.format(model, locs, disp)
    with open(bash, 'w') as f:
        f.write(to_call + '\n')
    pwd = os.getcwd()
    os.chdir(OUT)
    subprocess.call('chmod +x {}'.format(os.path.basename(bash)),shell=True)
    subprocess.call('./' + os.path.basename(bash),shell=True)
    os.chdir(pwd)

    df = pd.read_csv(disp, delim_whitespace=True,header=None)

    return df


def make_model(model,OUT):
    out = """MODEL.01
TEST MODEL
ISOTROPIC
KGS
FLAT EARTH
1-D
CONSTANT VELOCITY
LINE08
LINE09
LINE10
LINE11
  H(KM) VP(KM/S) VS(KM/S) RHO(GM/CC)   QP   QS  ETAP  ETAS  FREFP  FREFS\n"""
    H = np.ones(model.shape[0]) * np.diff(model[:, 0])[0]
    vp = model[:, 1] 
    vs = model[:, 2] 
    rho = model[:,3] 

    for ii in range(len(H)):
        layer = '{:04.2f}   {:02.1f} {:02.1f} {:02.1f}'.format(H[ii],vp[ii],vs[ii],rho[ii])
        layer += ' 0.0 0.0 0.0 0.0 1.0 1.0\n'
        out += layer
    last_layer = '{:03.1f} {:02.1f} {:02.1f} {:02.1f}'.format(0.,vp[ii],vs[ii],rho[ii])
    last_layer += ' 0.0 0.0 0.0 0.0 1.0 1.0\n'

    # out to model
    model96 = os.path.join(OUT,'CUS.mod')
    with open(model96, 'w') as f:
        f.write(out)

    return None


def make_sobs(OUT):
    sobs = os.path.join(OUT,'sobs.d')
    text = "  4.99999989E-03  4.99999989E-03   0.0000000      4.99999989E-03   0.0000000\n" + \
           "    1    2    2    2    4    4    4    0    1    0\n" + \
           "CUS.mod\n" + \
           "tdisp.d\n"

    with open(sobs, 'w') as f:
        f.write(text)
    return None


def make_tdisp(OUT,wave,vel_type,period,dispersion,flag='X',mode=0,err=0.001):
    tdisp = os.path.join(OUT,'tdisp.d')

    if not os.path.isfile(tdisp):
        open(tdisp,'a').close()
    with open(tdisp,'a') as f:
        text = 'SURF96 {} {} {}   {}     {:06.4f}     {:05.4f}     {:05.4f}\n'.format(wave,vel_type,flag,mode,period,dispersion,err)
        f.write(text)
    return None


def run_surf96(OUT):
    pwd = os.getcwd()
    os.chdir(OUT)
    subprocess.call('surf96 39',shell=True)
    subprocess.call('surf96 1',shell=True)
    subprocess.call('srfker96 > srfker96.txt',shell=True)
    subprocess.call('surf96 39',shell=True)

    # need to split output 
    files = open(os.path.join(OUT,'srfker96.txt'),'r').read().split('___________________________________________________________________________________________')
    files = [f for f in files if f != ''] # remove blanks
    files = [f for f in files if f != '\n'] # remove end char
    files = [f[1:] for f in files] # remove starting end character
    for file in files:
        file_name = file.split('\n')[0]
        wave_type,period = file_name.split(':')
        wave_type = list(filter(None,wave_type.split(' ')))
        period = list(filter(None,period.split(' ')))
        out_file = '_'.join([wave_type[0],wave_type[1],period[1]]).lower()
        out_file = os.path.join(OUT,out_file + '.txt')
        with open(out_file, 'w') as f:
            f.write(file)

    os.chdir(pwd)
    return None   

if __name__ == "__main__":
    stadf = pd.read_csv("/home/timclements/CALI/CIstations.csv")
    step = 50
    min_depth, max_depth = 0, 4050
    wave = 'R'
    vel_type = 'U' 
    dispersion = 3.0

    # output directory 
    OUT = '/media/FOUR/data/DISPERSION'

    # community velocity model 
    # change this to where the cvm model of your choice is installed on your machine
    # if having trouble installing the cvm on your machine
    # we recommend using http://moho.scec.org/UCVM_web/web/viewer.php
    cvm='/Users/thclements/CVM-SI/model/cvmsi/'

    # get dispersion for any station 
    station = "LJR"
    ind = stadf.loc[stadf['Station'] == station].index.values[0]
    df = call_cvsmi(stadf.loc[ind,"Longitude"], stadf.loc[ind, "Latitude"], min_depth, max_depth, cvm, OUT, step=step)
    model = df.values.T[[2,6,7,8], :].T / 1000
    make_model(model,OUT)

    # make dispersion curves at different periods 
    # 0.5, 1.0, 2.0, 4.0, 8.0 Hz
    for period in [2.0 , 1.0, 0.5, 0.25, 0.125]:
        make_tdisp(OUT,wave,vel_type,period,dispersion,flag='X',mode=0,err=0.001)

    # make control file
    make_sobs(OUT)

    # run model 
    run_surf96(OUT)