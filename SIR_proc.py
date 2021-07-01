'''
    Author: Nicholas J. Sisco, Ph.D. 
    Email: X@barrowneuro.org where X = nicholas.sisco
    Affiliation: Barrow Neurological Institue in Phoenix, AZ

    Title: Selective inversion recovery fitting using least squares and Levenberg-Marquardt 
    - This is a basic script written in Python for preprocessing signal acquired from SIR-qMT and calling Julia to fit to a double exponential 
    equation described in (1,2,3,4). 


    Requirements:
    - An OS with a unix-like environment, i.e. Linux, MacOS, or Windows Subsystem Linux
    - A working ANTs install with $ANTSPATH set
    - Python3.8
    - Julia 1.5

    Input and Output Parameters for Code

    -Currently this has only been tested with PAR files. It can easily be used with Dicom/nii files, but will need to be coded.

    User supplied data:
        - Command line usage example
        $ python ./SIR_proc_20210617.py -i <FULL_PATH>/qMT_.PAR

        - In the Python processed directory, there should be nt number of preprocessed nifti files, corresponding 
        to the number of dynamic time points in the data, e.g. nt = 4 means there are ti = [15,15,...,...], etc.

    Output:
        - The output from the Python portion of this file will be motion corrected SIR-qMT. One can change this to denoise and N4bias correct, which has been coded below but is turned off.


    Future Updates TO Do:
        1) user defined ti and td
        2) user defined kmf
        3) user defined preprocessing niftis, i.e. flexible input names
        4) Need to add dicom and nii functionality.
        5) module creation
        6) unit testing
        7) depolyable docker

    References:
    1. R. D. Dortch, J. Moore, K. Li, M. Jankiewicz, D. F. Gochberg, J. A. Hirtle, J. C. Gore, S. A. Smith, Quantitative magnetization transfer imaging of human brain at 7T. Neuroimage. 64, 640–649 (2013).
    2. F. Bagnato, G. Franco, F. Ye, R. Fan, P. Commiskey, S. A. Smith, J. Xu, R. Dortch, Selective inversion recovery quantitative magnetization transfer imaging: Toward a 3 T clinical application in multiple sclerosis. Mult. Scler. J. 26, 457–467 (2020).
    3. R. D. Dortch, F. Bagnato, D. F. Gochberg, J. C. Gore, S. A. Smith, Optimization of selective inversion recovery magnetization transfer imaging for macromolecular content mapping in the human brain. Magn. Reson. Med. 80, 1824–1835 (2018).
    4. R. D. Dortch, K. Li, D. F. Gochberg, E. B. Welch, A. N. Dula, A. A. Tamhane, J. C. Gore, S. A. Smith, Quantitative magnetization transfer imaging in human brain at 3 T via selective inversion recovery. Magn. Reson. Med. 66, 1346–1352 (2011).

    Change log:
    "history (of nifti library changes):\n"
    "\n",
    "0.0  June 18, 2021 [nsisco]\n"
    "     (Nicholas J. Sisco, Ph.D. of Barrow Neurological Institue)\n"
    "   - initial version \n"

'''
import os
import nibabel as nib
import numpy as np
from datetime import date
import argparse


class VerboseStore(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None:
            raise ValueError('nargs not allowed')
        super(VerboseStore, self).__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        print('The Path to PAR/REC is: %r using %r option' % (values,option_string))
        setattr(namespace, self.dest, values)

my_parser = argparse.ArgumentParser()
my_parser.add_argument('-i', '--input', action=VerboseStore, type=str,
                    help='Full Path of the SIR PAR/REC. Will eventually accept dicoms and niftis.'
                    )

args = my_parser.parse_args()

USR_INPUT = args.input
path = os.path.split(USR_INPUT)[0]
fname = USR_INPUT

today = date.today()
label=today.strftime("%Y%m%d")

out_dir = os.path.join(path,'output_{}/'.format(label))
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
img_fname = os.path.join(out_dir,'SIR_{}.nii.gz'.format(label))
print("Base Path for file => {}".format(path))
print("New Output Directory => {}".format(out_dir))
print("Image Name => {}".format(fname))
print("New Image Name => {}".format(img_fname))

base_dir = out_dir
rawData = nib.load(fname)
nx,ny,nz,ntmp = rawData.shape
reshaped_data = rawData.get_fdata().reshape(nx,ny,nz,3,4) # verify with matlab
new_sir_data = np.squeeze(reshaped_data[:,:,:,0,:])
header = rawData.header
affine = rawData.affine
new_img = nib.Nifti1Image(new_sir_data,affine,header=header)
nib.save(new_img,img_fname)
sir_fname = img_fname
sir_fname = base_dir+'SIR_{}.nii.gz'.format(label)

mo_corr_path = os.path.join(out_dir,'SIR_mo_corr')


# Denoising is not working right now for SIR
# I think it is actually the N4 or a combination of both. The N4 is definitely making a mess of the fitting. It's only used for registration now.
N4 = False
DENOISE = False
if N4:
    if DENOISE:
        print('Doing Rician Model Denoise with ANTs')
        denoised = base_dir+'SIR_{}_noise_corrected.nii.gz'.format(label)
        noise_img = base_dir+'SIR_{}_noise.nii.gz'.format(label)

        if not os.path.isfile(denoised):
            cmd = 'DenoiseImage -d 4 -i {} -n Rician -o [{},{}] -s 3 -r 8 -p 4 -v'.format(sir_fname,denoised,noise_img,label)
            print(cmd)
            os.system(cmd)
        
        def minmax(path):
            vol=nib.load(path).get_fdata().astype('float32')
            min=vol.min()
            max=vol.max()
            return min,max

        min,max = minmax(denoised)
        rescaled_dn_img = base_dir+'SIR_{}_noise_corrected_rescaled.nii.gz'.format(label)
        cmd = 'ImageMath 4 {} RescaleImage {} {} {}'.format(rescaled_dn_img,denoised,min,max)
        os.system(cmd)

        n4corrected = base_dir+'SIR_{}_noise_corrected_n4Corr.nii.gz'.format(label)
        n4corrected_f = base_dir+'SIR_{}_noise_corrected_n4Corr_field.nii.gz'.format(label)
        cmd = 'N4BiasFieldCorrection -d 4 --input-image {} --convergence [ 100x100x100x50, 1e-06 ] --output [{},{}] --shrink-factor 3'.format(denoised,n4corrected,n4corrected_f)
        print(cmd)
        os.system(cmd)
    else:        
        def minmax(path):
            vol=nib.load(path).get_fdata().astype('float32')
            min=vol.min()
            max=vol.max()
            return min,max

        min,max = minmax(sir_fname)
        rescaled_img = base_dir+'SIR_{}_rescaled.nii.gz'.format(label)
        cmd = 'ImageMath 4 {} RescaleImage {} {} {}'.format(rescaled_img,sir_fname,10,100)
        os.system(cmd)

        n4corrected = base_dir+'SIR_{}_n4Corr.nii.gz'.format(label)
        n4corrected_f = base_dir+'SIR_{}_n4Corr_field.nii.gz'.format(label)
        cmd = 'N4BiasFieldCorrection -d 4 --input-image {} --convergence [ 100x100x100x50, 1e-06 ] --output [{},{}] --shrink-factor 4'.format(rescaled_img,n4corrected,n4corrected_f)
        print(cmd)
        os.system(cmd)

        rescaled_img = base_dir+'SIR_{}_rescaled.nii.gz'.format(label)
        cmd = 'ImageMath 4 {} RescaleImage {} {} {}'.format(rescaled_img,sir_fname,min,max)
        os.system(cmd)

# Motion Correction
def mo_corr(in_fname,o_fname):
    cmd = 'antsMotionCorr -d 3 -n 10 -a {} -o {}_avg.nii.gz -v'.format(in_fname,o_fname)
    os.system(cmd)

    cmd = 'antsMotionCorr  -d 3 -o [{},{}.nii.gz,{}_avg.nii.gz] -m gc[ {}_avg.nii.gz , {} , 1 , 1 , Random, 0.05  ] -t Rigid[ 0.005 ] -i 20 -u 1 -e 1 -s 0 -f 1 -n 10'.format(o_fname,o_fname,o_fname,o_fname,in_fname)
    print(cmd)
    os.system(cmd)


if not os.path.isfile(mo_corr_path + '.nii.gz'):
    mo_corr(sir_fname,mo_corr_path)
else:
    print('Skipping Motion Correction')

proc_dir = base_dir+'/proc_{}/'.format(label)
if not os.path.isdir(proc_dir):
    os.mkdir(proc_dir)


''' 
Splitting up the files
'''
info=nib.load(mo_corr_path + '.nii.gz')
DATA = info.get_fdata().astype('float64')
nx,ny,nz,nt = DATA.shape
header = info.header
affine = info.affine
for ii in range(nt):
    new_data = DATA[:,:,:,ii]
    new_img = nib.Nifti1Image(new_data,affine,header=header)
    nib.save(new_img,proc_dir+'preprocessed_{}.nii.gz'.format(ii))


os.system('cp {} {}'.format(base_dir+'preprocessed_{}.nii.gz'.format(0),proc_dir+'/preprocessed_{}.nii.gz'.format(0)))
for ii in [1,2,3]:
    os.system('cp {} {}'.format(base_dir+'preprocessed_{}_2_0.nii.gz'.format(0),proc_dir+'/orig_{}_2_0.nii.gz'.format(0)))
cmd = 'bet {} {}brain.nii.gz -m '.format(proc_dir+'preprocessed_{}.nii.gz'.format(0),proc_dir)
print(cmd)
os.system(cmd)

cmd = 'julia SIR_PING_Brain_20210617_v2.jl %s' % proc_dir
print(cmd)
os.system(cmd)