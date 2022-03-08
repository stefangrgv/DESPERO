###################################################################
#
#   Functions for the
#   Data reduction pipeline for ESPERO @ NAO Rozhen
#
#
#   Stefan Georgiev
#   stefan.grgv@gmail.com
#
#   last modification: 2 September 2020
#
#
###################################################################

import numpy as np
import datetime
import shutil
from astropy.io import fits
import os
import sys
import tkinter as tk
tk.Tk().withdraw()

from tkinter.filedialog import askdirectory
from tkinter.filedialog import askopenfilename

print('Loading PYRAF...')
from pyraf import iraf

fontname = 'Times New Roman'
fontsize = 18

def update_progressbar(current_task, **kwargs):
    global progress, progress_label, progress_file

    progress_label.configure(text=current_task, font=(fontname, fontsize))
    f = ''
    if 'f' in kwargs:
        f = kwargs['f']

    progress_file.configure(text=f, fg='white', bg='blue', width=progress_width, font=(fontname, fontsize))

    progress.update_idletasks()
    progress.update()


def create_progressbar():
    global progress, progress_label, progress_file, progress_width
    progress = tk.Tk()
    progress_width = progress.winfo_screenwidth()
    progress_height = fontsize*3
    progress.geometry('%sx%s' % (progress_width, progress_height))
    #progress.overrideredirect(True)
    progress.title('DESPERO')
    #progressbar = Progressbar(progress, orient=HORIZONTAL, length=300,mode='determinate')
    progress_label = tk.Label(progress)
    progress_label.pack()
    progress_file = tk.Label(progress, width=progress_width)
    progress_file.pack()

def set_bias_flat(b, f):
    global boolean_bias, boolean_flat

    boolean_bias = b
    boolean_flat = f

def set_telluric(correct):
    global telluric_correct
    telluric_correct = correct

def get_directory():
    update_progressbar(current_task='Select OBSERVATIONS folder')
    global working_directory
    working_directory = askdirectory(title='Select OBSERVATIONS folder')
    if not working_directory:
        print('Observations folder not selected, exiting...')
        sys.exit(1)

    return(working_directory)

def load_packages():
    update_progressbar(current_task='Loading IRAF packages')
    print('Loading IRAF packages...')
    iraf.rv()
    iraf.imred()
    iraf.ccdred()
    iraf.echelle()
    iraf.imutil()

def load_parameters(procedure, filename):
    #sorry for the non-pythonic approach, but numpy doesn't handle these files properly
    with open(os.path.join( os.getcwd(), os.path.join('settings', filename) ), 'r') as f:
        for line in f:
            parameter = str(line.split(' = ')[0]).strip(' ')
            value = str(line.split(' = ')[1].split('#')[0].strip(' '))
            procedure[parameter] = value

def pyraf_load_settings():
    update_progressbar(current_task='Loading IRAF settings')
    global parameters_ccdproc, parameters_zerocombine, parameters_flatcombine 
    parameters_ccdproc = {}
    parameters_zerocombine = {}
    parameters_flatcombine = {}

    load_parameters(parameters_ccdproc,  'ccdproc.ini')
    for par, value in parameters_ccdproc.items():
        iraf.ccdproc.setParam(par, value)

    load_parameters(parameters_zerocombine, 'zerocombine.ini')
    for par, value in parameters_zerocombine.items():
        iraf.zerocombine.setParam(par, value)

    load_parameters(parameters_flatcombine, 'flatcombine.ini')
    for par, value in parameters_flatcombine.items():
        iraf.flatcombine.setParam(par, value)

    iraf.ccdred.setParam('instrument','NULL')
    iraf.echelle.setParam('dispaxi', '1')

    os.chdir(working_directory)

def check_for_existing_file(fname):
    f = Path(fname)
    if f.is_file():
        input('\n\n/!\ I\'M TRYING TO CREATE FILE <' + fname + '>, BUT IT ALREADY EXISTS. REMOVE IT AND PRESS ENTER TO CONTINUE OR CTRL+C TO ABORT PROGRAM\n')

def clean_header(filename):
    '''
    Removes unnecessary keys from header
    '''
    iraf.hedi(images=filename, fields='datasec', value='', delete='yes', ver='no')
    iraf.hedi(images=filename, fields='ccdsec', value='', delete='yes', ver='no')

def load_journal():
    update_progressbar(current_task='Loading journal')
    global obs_file, obs_date, obs_exp, obs_type
    #read the journal file
    obs_file, dates, exp, obs_type = np.loadtxt('Journal.txt', dtype='str', unpack=True)
    #convert str to corresponding type
    obs_date = [datetime.datetime.strptime(date, '%Y-%m-%dT%H:%M:%S') for date in dates]
    obs_exp = [float(ex) for ex in exp]

def obtain_rdtime():
    '''
    Obtain the readout time
    '''
    global obs_rdspeed, all_readspeeds
    obs_rdspeed, all_readspeeds = [], []

    for current_file in obs_file:
        fits_file = fits.open(current_file + '.fits')
        current_readtime = float(fits_file[0].header['READTIME'])
        fits_file.close()

        obs_rdspeed.append(current_readtime)
        if current_readtime not in all_readspeeds:
            all_readspeeds.append(current_readtime)

def sort_fits_type():
    '''
    Sorts the .fits according to the type of exposure: bias, flat or comparison spectrum
    '''

    global list_bias, list_flat, list_comp
    list_bias = []
    list_flat = []
    list_comp = []

    for j in range(len(obs_file)):
        if obs_type[j] == 'zero':
            list_bias.append([obs_file[j], obs_rdspeed[j]])
        elif obs_type[j] == 'flat':
            list_flat.append([obs_file[j], obs_rdspeed[j]])
        elif obs_type[j] == 'comp':
            list_comp.append([obs_file[j], obs_rdspeed[j]])

        #clean unnecesarry fields in the header that interfere with PYRAF
        clean_header(obs_file[j])

def rv_correct():
    '''
    Heliocentric radial velocity correction
    '''

    update_progressbar(current_task='Correcting for heliocentric radial velocity')
    for j in range(len(obs_file)):
        update_progressbar(current_task='Correcting for heliocentric radial velocity', f=obs_file[j])
        iraf.rvcorrect(images=obs_file[j], imup='yes', observatory='rozhen')

def make_master_bias_flat():
    global master_flats
    master_flats = []

    update_progressbar(current_task='Building master bias & flat')
    if boolean_bias:
        for speed in all_readspeeds:
            current_bias_speed = '' #a list of all the bias files with the current readspeed; pyraf's zerocombine requires a list input, that's the reason for this
            for bias in list_bias:
                if bias[1] == speed:
                    current_bias_speed += bias[0] + '.fits, '
            current_bias_speed = current_bias_speed[0:-2] #cut the ', ' at the end

            if current_bias_speed:
                master_bias_name = 'master_bias-' + str(speed) + '.fits'

                print('Creating master bias for rdspeed ' + str(speed) + '...')
                iraf.zerocombine(input=current_bias_speed, output=master_bias_name)
                clean_header(master_bias_name)
                
            else:
                print('\nWARNING:\n\tSOME FILES HAVE RDSPEED = ' + str(speed) + '\n\tNO BIAS FILE IS AVAILABLE FOR THIS RDSPEED')
                input('(press any key to continue)')

    if boolean_flat:
        for speed in all_readspeeds:
            current_flat_speed = '' #same reason
            for flat in list_flat:
                if flat[1] == speed:
                    current_flat_speed += flat[0] + '.fits, '
            current_flat_speed = current_flat_speed[0:-2]

            if current_flat_speed:
                master_flat_name = 'master_flat-' + str(speed) + '.fits'
                master_flats.append(master_flat_name)
                print('Creating master flat for rdspeed ' + str(speed) + '...')
                iraf.flatcombine(input=current_flat_speed, output=master_flat_name)
                clean_header(master_flat_name)
            else:
                print('\nWARNING:\n\tSOME FILES HAVE RDSPEED = ' + str(speed) + '\n\tNO FLAT FILE IS AVAILABLE WITH THIS RDSPEED')
                input('(press any key to continue)')
                
def correct_for_bias(filenum):
    update_progressbar(current_task='Correcting for bias')
    masterbias = 'master_bias-' + str(obs_rdspeed[filenum]) + '.fits'

    if os.path.exists(masterbias):
        update_progressbar(current_task='Correcting for bias', f=obs_file[filenum])
        print(str(obs_file[filenum]) + '.fits corrected for bias... ')
        iraf.ccdproc(images=obs_file[filenum] + '.fits', zerocor='yes', zero=masterbias)
    else:
        print('>\t' + str(obs_file[filenum]) + ' not corrected for bias.')

def correct_master_flat_with_bias():
    for speed in all_readspeeds:
        masterbias = 'master_bias-' + str(speed) + '.fits'
        masterflat = 'master_flat-' + str(speed) + '.fits'
        if os.path.exists(masterbias) and os.path.exists(masterflat):
            print('Master flat ' + str(speed) + ' corrected for bias... ')
            iraf.ccdproc(images=masterflat, zerocor='yes', zero=masterbias)
        else:
            print('Master flat ' + str(speed) + ' not corrected for bias or does not exist.')
            

def correct_for_flat(filenum):
    update_progressbar(current_task='Correcting for flat')
    masterflat = 'master_flat-' + str(obs_rdspeed[filenum]) + '.ec.fits'

    clean_header(masterflat)

    if os.path.exists(masterflat):
        update_progressbar(current_task='Correcting for flat', f=obs_file[filenum])
        clean_header(obs_file[filenum] + '.ec.fits')

        print(str(obs_file[filenum]) + '.ec.fits corrected for flat... ')
        iraf.ccdproc(images=obs_file[filenum] + '.ec.fits', flatcor='yes', flat=masterflat)
    else:
        print('>\t' + str(obs_file[filenum]) + ' not corrected for flat.')

def select_thar(mode):
    '''
    For each observation, select the comparison spectrum which is most suited for it
    '''

    update_progressbar(current_task='Selecting ThAr spectra')
    global thar
    thar = []

    if mode == 'closest': # picks the comparsion spectrum which is closest in time to the observation
        for j in range(len(obs_file)):
            if obs_type[j] == 'object':
                first_thar_index, second_thar_index = -1, -1
                for i in range(j, -1, -1):
                    if obs_type[i] == 'comp':
                        first_thar_index = i
                        break
                for i in range(j, len(obs_file)):
                    if obs_type[i] == 'comp':
                        second_thar_index = i
                        break
                
                if first_thar_index != -1 and second_thar_index != -1:
                    first_thar_timedifference = (obs_date[j] - obs_date[first_thar_index]).total_seconds()
                    second_thar_timedifference = (obs_date[j] - obs_date[second_thar_index]).total_seconds() + obs_exp[j]

                    if first_thar_timedifference < second_thar_timedifference:
                        final_thar = obs_file[first_thar_index]
                    else:
                        final_thar = obs_file[second_thar_index] + '.ec.fits'
                else:
                    if first_thar_index != -1:
                        final_thar = obs_file[first_thar_index] + '.ec.fits'
                    elif second_thar_index != -1:
                        final_thar = obs_file[second_thar_index] + '.ec.fits'
                    else:
                        print('ERROR: FOUND NO COMP SPECTRUM FOR FILE ' + str(obs_file[j]) + '!')

                print(obs_file[j] + '\t' + final_thar)
                thar.append(final_thar.split('.')[0])
                iraf.hedi(images=obs_file[j] + '.fits', fields='refspec1', value=str(final_thar.split('.fits')[0]), add='yes', addonly='no', delete='no', verify='no', show='yes', update='yes')

def apall_reference():
    # apall is short for aperture allocate; this iraf procedure finds the pixels in the .fits where the actual spectrum is found
    update_progressbar(current_task='Select APALL reference file')
    global apall_ref
    apall_path = askopenfilename(title='Select APALL reference file')
    if not apall_path:
        print('APALL reference file not selected, exiting...')
        sys.exit(1)

    apall_ref = apall_path.split('/')[-1]

    iraf.imred.echelle.apall(apall_ref, extract='yes', recenter='no', resize='no', trace='yes', fittrac='yes', review='no', nfind=68, t_funct='spline3', t_order=1)

def apall_spectra(apref):
    if apref:
        global apall_ref
        apall_ref = apref

    for j in range(len(thar)):
        clean_header(thar[j])

    for j in range(len(obs_file)):
        if obs_type[j] == 'object':
            update_progressbar(current_task='Extracting echelle aperture', f=obs_file[j])
            print('Extracting echelle spectrum from file ' + obs_file[j] + '...')
            iraf.imred.echelle.apall(input=obs_file[j] + '.fits', referen=apall_ref, edit='no', trace='no', extract='yes', resize='no', recente='no', review='no', interac='no')

    solved = []
    for j in range(len(thar)):
        if str(thar[j]) not in solved:
            update_progressbar(current_task='Extracting echelle aperture', f=thar[j])
            print('Extracting echelle spectrum from file ' + thar[j] + '...')
            iraf.imred.echelle.apall(input=thar[j], referen=apall_ref, edit='no', trace='no', extract='yes', resize='no', recente='no', review='no', interac='no')
            solved.append(str(thar[j]))

def apall_master_flat():
    for mf in master_flats:
        update_progressbar(current_task='Extracting echelle aperture', f=mf)
        print('Extracting echelle spectrum from file ' + mf + '...')
        iraf.imred.echelle.apall(input=mf, referen=apall_ref, edit='no', trace='no', extract='yes', resize='no', recente='no', review='no', interac='no')


def solve_comp(default):
    '''
    Prompts for the solved spectrum that will later be used for the wavelength calibration of the observations
    '''
    if not default:
        update_progressbar(current_task='Select solved spectrum')
        wl_solution_path = askopenfilename(title='Select solved spectrum', filetypes=[('.EC files', '*.ec')])
        if not wl_solution_path:
            print('Solved spectrum not selected, exiting...')
            sys.exit(1)

        wl_solution = wl_solution_path.split('/')[-1][2:]

    else:
        wl_solution = 'refspectrum.ec'

    solved = []
    for j in range(len(thar)):
        if str(thar[j]) not in solved:
            update_progressbar(current_task='Solving ThAr for wavelength', f=thar[j])
            iraf.ecreidentify(images=(thar[j] + '.ec'), referenc=wl_solution, shift=0., cradius=5., threshold=10., refit='yes', database='database')
            solved.append(str(thar[j]))

def telluric():
    '''
    Clean telluric lines
    '''
    update_progressbar(current_task='Cleaning telluric lines')

    telluric_standard = askopenfilename(title='Select telluric standard file', filetypes=[('.EC.FITS files', '*.ec.fits')])
    if not telluric_standard:
        print('Telluric file not selected, exiting...')
        sys.exit(1)

    telluric_standard = telluric_standard.split('/')[-1]

    for j in range(len(obs_file)):
        if obs_type[j] == 'object':
            print('Cleaning telluric lines from file ' + obs_file[j] + '...')
            update_progressbar(current_task='Cleaning telluric lines', f=obs_file[j])
            
            corrected_file = obs_file[j] + '.tel'
            iraf.imarith(obs_file[j] + '.ec.fits', '/', telluric_standard, corrected_file)

def calibrate_wl():
    global wl_fits
    wl_fits = []
    update_progressbar(current_task='Wavelength calibration')

    #move all spectra to a seperate folder
    if not os.path.isdir( os.path.join(working_directory, 'spectra') ):
        os.mkdir( os.path.join(working_directory, 'spectra') )

    for j in range(len(obs_file)):
        if obs_type[j] == 'object':
            update_progressbar(current_task='Wavelength calibration', f =obs_file[j])
            obs_file_wl = obs_file[j] + '.wl.fits'

            if telluric_correct:
                extention = '.tel.fits'
            else:
                extention = '.ec.fits'

            iraf.dispcor(input=obs_file[j] + extention, output=obs_file_wl, lineari='yes', databas='database', table='', w1='INDEF', w2='INDEF', dw='INDEF', nw='INDEF', flux='yes', blank=0., samedisp='no', ignoreaps='no', confirm='no', listonly='no', verbose='yes')
            iraf.dopcor(input=obs_file_wl, output=obs_file_wl, redshift='-vhelio', isvelocity='yes', add='no', dispersion='yes', flux='no', factor=3., apertures='', verbose='yes')

#            shutil.move(os.path.join(working_directory, obs_file_wl), os.path.join(working_directory, 'spectra', obs_file_wl))

#            wl_fits.append( os.path.join(working_directory, 'spectra', obs_file_wl) )

def fits_to_ascii():
    '''
    Converts the .wl.fits files into ASCII spectra
    '''

    n_apertures = 69 # number of apertures in ESpeRo

    update_progressbar(current_task='Converting to ASCII')
    for obs_file_wl in wl_fits:
        update_progressbar(current_task='Converting to ASCII', f=obs_file_wl)
        print('Converting ' + obs_file_wl.split('/')[-1].split('.')[0] + ' to ASCII...')

        aps = []
        for ap in range(1,n_apertures):
            output_fits = obs_file_wl + '-' + str(ap)
            iraf.scopy(obs_file_wl, output_fits, apertures=str(ap))

            output = obs_file_wl.split('.wl.fits')[0] + '-' + str(ap) + '.%04i.txt'%ap
            iraf.wspectext(obs_file_wl + '-' + str(ap) + '.%04i.fits'%ap, output, header='no')
            aps.append(output)

            os.remove(obs_file_wl + '-' + str(ap) + '.%04i.fits'%ap)

        obs_file_ascii = obs_file_wl.replace('.wl.fits', '.es')

        #to ensure the file is empty
        f = open(obs_file_ascii, 'w')
        f.close()

        with open(obs_file_ascii, 'a') as f:
            for ap in reversed(aps):
                with open(ap, 'r') as inf:
                    f.write(inf.read())

        #remove the ASCII files of the individual apertures
        for ap in aps:
            os.remove(ap)

def exit():
    print('All done.')
    sys.exit(0)

def dobias():
    for j in range(len(obs_file)):
        if obs_type[j] != 'zero' and obs_type[j] != 'flat':
            correct_for_bias(j)

def doflat():
    for j in range(len(obs_file)):
            if obs_type[j] != 'zero' and obs_type[j] != 'flat':
                correct_for_flat(j)

