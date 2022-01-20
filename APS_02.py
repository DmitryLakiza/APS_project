import random
import astropy

import matplotlib.pyplot as plt
import numpy as np

from astropy.io import fits
from formulas import gauss
from astropy.modeling import models, fitting
from math import sqrt, pi, log
from scipy.signal import savgol_filter as savgol
import scipy.integrate as integrate


my_file = fits.open(r'D:\projectAPS\objects\spec-1373-53063-0583.fits')

# finding dates from a file
data = my_file[1].data

# enumerating table columns
columns = my_file[1].columns

# assigning column data to variables
flux = data.field('flux')
loglam = data.field('loglam')
sky = data.field('sky')
model = data.field('model')

requests_list = ('lam graph',
                 'loglam graph',
                 'lam graph model',
                 'loglam graph model',
                 'gauss lam graph model',
                 'lam graph and model',
                 "analyzed line (requests completion of number 6 or 7) (gauss' info)",
                 'histogram',
                 'filtering flux',
                 'end')

loglam_axis = np.linspace(min(loglam), max(loglam), len(loglam))


def print_help(my_list):
    for i in range(1, len(my_list)+1):
        print(f'{i}. {my_list[i-1]}')
    return None


def log_graph(loglam_axis=loglam_axis, value=flux, model=0):
    plt.plot(loglam_axis, value-model)
    return None


def graph(loglam_axis=loglam_axis, value=flux, model=0):
    plt.plot(10**loglam_axis, value-model)
    return None


def gauss_graph(left, right, model=0):  # !available just for /lam graph model\
    while 1:
        gauss_request = input('Handmade or auto?: ')
        if gauss_request == 'handmade' or gauss_request == '1':

            coeffs = tuple(map(float, input('Input coeffs (a,b,c): ').split()))
            #print(coeffs)
            layout = np.linspace(left, right, 1500)
            #for item in layout:
            #    print(gauss(item, coeffs[0], coeffs[1], coeffs[2]))
            plt.plot(layout, gauss(layout, coeffs[0], coeffs[1], coeffs[2]))
            graph(loglam_axis, flux, model)
            plt.show()
            gauss_request = input('ok?: ')
            if gauss_request == 'yes':
                break
            elif gauss_request == 'no':
                continue
            else:
                continue

        elif gauss_request == 'auto' or gauss_request == '2':

            layout = np.linspace(left, right, 1500)

            lam_range = (left <= 10**loglam) & (10**loglam <= right)

            auto_gauss_flux = (flux-model)[lam_range]
            auto_gauss_loglam = 10**loglam[lam_range]

            some_sum = 0
            whole_flux = 0

            for idx in range(len(auto_gauss_flux)):
                some_sum += auto_gauss_flux[idx] * auto_gauss_loglam[idx]
                whole_flux += auto_gauss_flux[idx]

            average_flux = some_sum / whole_flux
            LSQ_sum = 0
            for idx in range(len(auto_gauss_loglam)):
                LSQ_sum += auto_gauss_flux[idx] * ((auto_gauss_loglam[idx] - average_flux) ** 2)
            sigma = sqrt(LSQ_sum / whole_flux)

            a = max(auto_gauss_flux) + (1 / (sigma * sqrt(2*pi)))
            b = average_flux
            c = sigma

            result = integrate.quad(lambda x: a*np.exp(-( ( (x-b) ** 2)/(2 * (c**2))) ), 6730, 6742)
            print(result)

            graph(loglam_axis, flux, model)
            plt.plot(layout, gauss(layout, a, b, c))
            plt.show()
            gauss_request = input('ok?: ')
            if gauss_request == 'yes':
                break
            elif gauss_request == 'no':
                continue
            else:
                continue
    return a, b, c


def savgol_filtering(flux):
    window = int(input('Input window length (it must be odd number): '))
    power = int(input("Input polynom's power(it must be less "
                      "or equal to 'window': "))
    filtered_flux = savgol(flux, window, power)
    return filtered_flux


def floating_average(flux):
    flux = list(flux)
    filtered_flux = []
    filtered_flux.append((flux[0] + flux[1]) / 2)

    for i in range(1, len(flux)-1):
        average_flux_by_i = (flux[i-1] + flux[i] + flux[i+1]) / 3
        filtered_flux.append(average_flux_by_i)
    filtered_flux.append((flux[-1] + flux[-2]) / 2)
    flux = np.array(flux)
    filtered_flux = np.array(filtered_flux)
    return filtered_flux


print_help(requests_list)
while 1:
    request = input('Input number of operation: ')
    if request == '1':
        graph(loglam_axis, flux)
        plt.show()
    elif request == '2':
        log_graph(loglam_axis, flux)
        plt.show()
    elif request == '3':
        graph(loglam_axis, flux, model)
        plt.show()
    elif request == '4':
        log_graph(loglam_axis, flux, model)
        plt.show()
    elif request == '5':
        graph(loglam_axis, flux, model)
        plt.show()
        left = float(input('The end of left side is at (Angstroms): '))
        right = float(input('The end of right side is at (Angstroms): '))
        coeffs = gauss_graph(left, right, model)
    elif request == '6':
        graph(loglam_axis, flux)
        graph(loglam_axis, model)
        plt.show()
    elif request == '7':
        print(f'a = {a}')
        print(f'b = {b}')
        print(f'c = {c}')
        print(f'left = {left}')
        print(f'right = {right}')
    elif request == '8':
        lam_range = (left <= 10**loglam) & (10**loglam <= right)
        shoots = int(input('Input number of shoots: '))
        x_coord_of_hist = []
        for i in range(shoots):
            x_coord_of_hist.append(random.gauss(coeffs[1], coeffs[2]))
        x_coord_of_hist.sort()
        x_coord_of_hist = np.array(x_coord_of_hist)
        graph(loglam_axis, flux, model)
        layout = np.linspace(left, right, 1500)
        plt.plot(layout, gauss(layout, coeffs[0], coeffs[1], coeffs[2]))
        plt.hist(x_coord_of_hist, bins=len(10**loglam[lam_range])+2, rwidth=2)
        plt.show()

    elif request == '9':
        ask = input('Savgol or floating average?: ')
        if ask == '1':
            filtered_flux = savgol_filtering(flux-model)
        elif ask == '2':
            filtered_flux = floating_average(flux-model)
        graph(loglam_axis, flux-model)
        graph(loglam_axis, filtered_flux)
        plt.show()
    elif request == '10':
        break
    else:
        continue
    print('-'*40)
