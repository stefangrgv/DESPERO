###################################################################
#
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

import functions as despero
import tkinter as tk
import shutil
import os

root = tk.Tk()
bias_correct = False
flat_correct = False
auto_apall = False
default_ec = False
telluric_correct = False
ascii = False

def go():
    '''
    This function runs the pipeline
    '''
    root.destroy()

    despero.create_progressbar()
    despero.load_packages()

    # copy the reference spectrum solution
    script_directory = os.getcwd()
    script_database  = os.path.join(script_directory, 'database')

    working_directory = despero.get_directory()
    working_database  = os.path.join(working_directory, 'database')

    if script_directory == working_directory:
        print('DESPERO is in the same directory as your observations. Please move the pipeline outside of there and re-run the script.')
        return(1)
    if not os.path.isdir(working_database):
        os.mkdir(working_database)
    shutil.copy( os.path.join(script_database, 'ecrefspectrum.ec'), os.path.join(working_database, 'ecrefspectrum.ec') )
    apref = ''
    if auto_apall:
        apref = 'apallrefspectrum.fits'
        shutil.copy( os.path.join(script_database, 'apapallrefspectrum'), os.path.join(working_database, 'apapallrefspectrum') )
        shutil.copy( os.path.join(script_directory, 'apallrefspectrum.fits'), os.path.join(working_directory, 'apallrefspectrum.fits') )

    despero.pyraf_load_settings()
    despero.set_bias_flat(bias_correct, flat_correct)
    despero.load_journal()
    despero.obtain_rdtime()
    despero.sort_fits_type()
    despero.rv_correct()
    despero.make_master_bias_flat()
    despero.select_thar('closest')

    if bias_correct:
        despero.dobias()

    if bias_correct and flat_correct:
        despero.correct_master_flat_with_bias()

    if not auto_apall:
        despero.apall_reference()
    despero.apall_spectra(apref)

    if flat_correct:
        if auto_apall:
            despero.apall_master_flat()
        despero.doflat()
        
    despero.solve_comp(default_ec)

    despero.set_telluric(telluric_correct)
    if telluric_correct:
        despero.telluric()

    despero.calibrate_wl()

    if ascii:
        despero.fits_to_ascii()

    despero.exit()

# get screen width and height
ws = root.winfo_screenwidth() # width of the screen
hs = root.winfo_screenheight() # height of the screen

w = 350 # width for the Tk root
h = 200 # height for the Tk root

# calculate x and y coordinates for the Tk root window
x = (ws/2) - (w/2)
y = (hs/2) - (h/2)

root.geometry('%dx%d+%d+%d' % (w, h, x, y)) 
root.title('DESPERO')

root.resizable(0,0)

button_color = 'orange'

def switch_bias():
    global bias_correct
    if bias_correct:
        bias_correct = False
    else:
        bias_correct = True

def switch_flat():
    global flat_correct
    if flat_correct:
        flat_correct = False
    else:
        flat_correct = True    

def switch_apall():
    global auto_apall
    if auto_apall:
        auto_apall = False
    else:
        auto_apall = True    

def switch_ec():
    global default_ec
    if default_ec:
        default_ec = False
    else:
        default_ec = True    

def switch_telluric():
    global telluric_correct
    if telluric_correct:
        telluric_correct = False
    else:
        telluric_correct = True    

def switch_ascii():
    global ascii
    if ascii:
        ascii = False
    else:
        ascii = True    

# populate the GUI with buttons
bias_checkbutton = tk.Checkbutton(root, command=switch_bias, highlightthickness=0, text='BIAS correction').pack(side=tk.TOP)
flat_checkbutton = tk.Checkbutton(root, command=switch_flat, highlightthickness=0, text='FLAT correction').pack(side=tk.TOP)
apall_checkbutton = tk.Checkbutton(root, command=switch_apall, highlightthickness=0, text='Automatic aperture extraction').pack(side=tk.TOP)
ec_checkbutton = tk.Checkbutton(root, command=switch_ec, highlightthickness=0, text='Default wavelength solution').pack(side=tk.TOP)
telluric_checkbutton = tk.Checkbutton(root, command=switch_telluric, highlightthickness=0, text='Telluric correction').pack(side=tk.TOP)
ascii_checkbutton = tk.Checkbutton(root, command=switch_ascii, highlightthickness=0, text='Convert to ASCII').pack(side=tk.TOP)
go_button = tk.Button(root, highlightthickness=0, text='GO', command=go)

go_button.pack(side=tk.BOTTOM)

root.mainloop() 
