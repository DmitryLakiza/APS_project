import random
from math import sqrt, pi, log

import astropy
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from formulas import gauss
from astropy.modeling import models, fitting

# constants:
sqrt2pi = sqrt(2 * pi)
c = 299792.458


# functions for graphs:
def lam_graph():
    fig = plt.subplot()
    x = np.linspace(min(needed_loglam), max(needed_loglam), len(needed_flux))
    fig = plt.plot(10**x, needed_flux[range(len(needed_flux))])


def loglam_graph():
    fig = plt.subplot()
    x = np.linspace(min(needed_loglam), max(needed_loglam), len(needed_flux))
    fig = plt.plot(x, needed_flux[range(len(needed_flux))])


def lam_graph_model():
    fig = plt.subplot()
    x = np.linspace(min(needed_loglam), max(needed_loglam), len(needed_flux))
    fig = plt.plot(10 ** x, needed_flux[range(len(needed_flux))] - needed_model[range(len(needed_model))])


def loglam_graph_model():
    fig = plt.subplot()
    x = np.linspace(min(needed_loglam), max(needed_loglam), len(needed_flux))
    fig = plt.plot(x, needed_flux[range(len(needed_flux))] - needed_model[range(len(needed_model))])


# file opening
m_file = fits.open(r'D:\projectAPS\objects\spec-1373-53063-0583.fits')

# request module
requests_list = ('lam graph',
                 'loglam graph',
                 'element',
                 'lam graph model',
                 'loglam graph model',
                 'gauss lam graph',
                 'gauss lam graph model',
                 'lam graph and model',
                 "analyzed line (requests completion of number 6 or 7) (gauss' info)",
                 'histogram',
                 'end'
                 )

# finding dates from a file
needed_data = m_file[1].data

# enumerating table columns
columns = m_file[1].columns

# assigning column data to variables
needed_flux = needed_data.field('flux')
needed_loglam = needed_data.field('loglam')
needed_sky = needed_data.field('sky')
needed_model = needed_data.field('model')

# creating elements value according to their axes
# Halpha
ys_Halpha = []
xs_Halpha = []

# Oxygen at 5008 Angstrem
ys_O_Three = []
xs_O_Three = []

# right Nitrogen
ys_Nit_right = []
xs_Nit_right = []

# left Nitrogen
ys_Nit_left = []
xs_Nit_left = []

# left Sulfur
ys_S_left = []
xs_S_left = []

# right sulfur
ys_S_right = []
xs_S_right = []

# list of redshifts
cz_list = []

# making list of values of elements by OY axis
for i in range(int(min(needed_flux)) - 10, int(max(needed_flux)) + 10):
    # Halpha
    xs_Halpha.append(6564.6140)
    ys_Halpha.append(i)

    # Nitrogen right
    xs_Nit_right.append(6585.27)
    ys_Nit_right.append(i)

    # Oxygen_third
    xs_O_Three.append(5008.240)
    ys_O_Three.append(i)

    # Nitrogen left
    xs_Nit_left.append(6549.86)
    ys_Nit_left.append(i)

    # S right
    xs_S_right.append(6732.67)
    ys_S_right.append(i)

    # S left
    xs_S_left.append(6718.29)
    ys_S_left.append(i)

xs_Halpha = np.array(xs_Halpha)
ys_Halpha = np.array(ys_Halpha)

xs_Nit_right = np.array(xs_Nit_right)
ys_Nit_right = np.array(ys_Nit_right)

xs_O_Three = np.array(xs_O_Three)
ys_O_Three = np.array(ys_O_Three)

xs_Nit_left = np.array(xs_Nit_left)
ys_Nit_left = np.array(ys_Nit_left)

xs_S_left = np.array(xs_S_left)
ys_S_left = np.array(ys_S_left)

xs_S_right = np.array(xs_S_right)
ys_S_right = np.array(ys_S_right)

i = 1
print('Available requests are: ')
for item in requests_list:
    print(str(i) + '.', item, sep=' ')
    i += 1
print()
# deleting i
del i

while 1:
    request = input('What is needed?: ')
    if request == 'lam graph' or request == '1':
        lam_graph()
        plt.show()
        


    elif request == 'loglam graph' or request == '2':
        loglam_graph()
        plt.show()
        


    elif request == 'element' or request == '3':
        print("Next graph'll show the loglam!")

        loglam_graph()
        plt.show()
        # inputing the value of element's loglam
        loglam_element = float(input('Print the value of loglam of element: '))
        # choosing the element
        element = input('Choose the element: ')
        # searching for agreement
        if element == 'Halpha':
            Halpha_loglam_lab = 3.81720893
            cz = ((10 ** loglam_element - 10 ** Halpha_loglam_lab) / (10 ** Halpha_loglam_lab)) * c
            print(f'cz = {cz} km/s')
            delt_lam = loglam_element - Halpha_loglam_lab

        elif element == 'O':
            Oxygen_loglam_lab = 3.699685133
            cz = ((10 ** loglam_element - 10 ** Oxygen_loglam_lab) / (10 ** Oxygen_loglam_lab)) * c
            print(f'cz = {cz} km/s')
            delt_lam = loglam_element - Oxygen_loglam_lab

        elif element == 'N right':
            Nitrogen_right_loglam_lab = 3.818573586
            cz = ((10 ** loglam_element - 10 ** Nitrogen_right_loglam_lab) / (10 ** Nitrogen_right_loglam_lab)) * c
            print(f'cz = {cz} km/s')
            delt_lam = loglam_element - Nitrogen_right_loglam_lab

        elif element == 'N left':
            Nitrogen_left_loglam_lab = 3.816232017
            cz = ((10 ** loglam_element - 10 ** Nitrogen_left_loglam_lab) / (10 ** Nitrogen_left_loglam_lab)) * c
            print(f'cz = {cz} km/s')
            delt_lam = loglam_element - Nitrogen_left_loglam_lab

        elif element == 'S right':
            S_right_loglam_lab = 3.828187328
            cz = ((10 ** loglam_element - 10 ** S_right_loglam_lab) / (10 ** S_right_loglam_lab)) * c
            print(f'cz = {cz} km/s')
            delt_lam = loglam_element - S_right_loglam_lab

        elif element == 'S left':
            S_left_loglam_lab = 3.827258747
            cz = ((10 ** loglam_element - 10 ** S_left_loglam_lab) / (10 ** S_left_loglam_lab)) * c
            print(f'cz = {cz} km/s')
            delt_lam = loglam_element - S_left_loglam_lab

        """
        delt_lam = shows diffrence between real and seen lam
        """
        cz_list.append(cz)

        # creating of "normal" spectrum
        x_2 = np.linspace(min(needed_loglam) - delt_lam, max(needed_loglam) - delt_lam, len(needed_loglam))
        plt.plot(10 ** x_2, (needed_flux - needed_model)[range(len(needed_flux))])

        # putting in element's lines
        Halpha = plt.plot(xs_Halpha, ys_Halpha)
        Nit_right = plt.plot(xs_Nit_right, ys_Nit_right)
        O = plt.plot(xs_O_Three, ys_O_Three)
        Nit_left = plt.plot(xs_Nit_left, ys_Nit_left)
        S_left = plt.plot(xs_S_left, ys_S_left)
        S_right = plt.plot(xs_S_right, ys_S_right)
        plt.show()

        

    elif request == 'lam graph model' or request == '4':
        lam_graph_model()
        plt.show()
        
    elif request == 'loglam graph model' or request == '5':
        loglam_graph_model()
        plt.show()
        
    elif request == 'gauss lam graph' or request == '6':
        question_Gauss_1 = input('Handmade or auto?: ')
        # question_Gauss_1 - asks about way of building
        # question_Gauss_0 - asks for satisfying about the function itself
        if question_Gauss_1 == 'Handmade' or question_Gauss_1 == 'handmade':
            while True:

                print('According to graph, find ends and "h": ')
                lam_graph()
                plt.show()

                # creating new graph for building Gauss' function
                fig = plt.subplot()
                x = np.linspace(min(needed_loglam), max(needed_loglam), len(needed_flux))
                fig = plt.plot(10 ** x, needed_flux[range(len(needed_flux))])

                left_side = float(input('The end of left side is at (Angstroms): '))
                right_side = float(input('The end of right side is at (Angstroms): '))
                print('f(x) = a*e**(-((x-b)**2))/(2*(c**2)))')
                # coeffs
                x_a = float(input('Coef a equals: '))
                x_b = float(input('Coef b equals: '))
                x_c = float(input('Coef c equals: '))
                h = float(input('h equals: '))
                x_a = x_a - h
                
                # Gauss' formula
                y = lambda x_1: x_a * np.exp(-(((x_1 - x_b) ** 2) / (2 * (x_c ** 2)))) + h

                # creating Gauss' function
                x_1 = np.linspace(left_side, right_side, 1500)
                plt.plot(x_1, y(x_1))

                plt.show()

                # repeat if not satisfied
                question_Gauss_0 = input('Ok?: ')
                if question_Gauss_0 == 'no' or question_Gauss_0 == 'No':
                    pass
                elif question_Gauss_0 == 'yes' or question_Gauss_0 == 'Yes':
                    break
        elif question_Gauss_1 == 'Auto' or question_Gauss_1 == 'auto':
            while True:

                print('According to graph, find ends and "h": ')
                lam_graph()
                plt.show()
                fig = plt.subplot()
                x = np.linspace(min(needed_loglam), max(needed_loglam), len(needed_flux))
                fig = plt.plot(10 ** x, needed_flux[range(len(needed_flux))])

                left_side = float(input('Left end (Angstroms) equals: '))
                right_side = float(input('Right end (Angstroms) equals: '))
                h = float(input('h equals: '))

                p = ((10 ** needed_loglam) >= left_side) & ((10 ** needed_loglam) <= right_side)

                needed_flux_new = (needed_flux)[p]
                needed_loglam_new = 10 ** needed_loglam[p]

                sum_nx = 0
                sum_n = 0

                for i in range(len(needed_loglam_new)):
                    sum_nx += needed_flux_new[i] * needed_loglam_new[i]
                    sum_n += needed_flux_new[i]
                average_x = sum_nx / sum_n
                sum_xx_squared = 0
                for i in range(len(needed_loglam_new)):
                    sum_xx_squared += needed_flux_new[i] * ((needed_loglam_new[i] - average_x) ** 2)
                sigma = sqrt(sum_xx_squared / sum_n)

                x_a = max(needed_flux_new)
                x_c = sigma
                x_b = average_x
                
                y = lambda x_1: (x_a - h) * np.exp(-(((x_1 - x_b) ** 2) / (2 * (x_c ** 2)))) + h

                x_1 = np.linspace(left_side, right_side, 1500)

                plt.plot(x_1, y(x_1))

                plt.show()

                question_Gauss_0 = input('Ok?: ')
                if question_Gauss_0 == 'no' or question_Gauss_0 == 'No':
                    pass
                elif question_Gauss_0 == 'yes' or question_Gauss_0 == 'Yes':
                    break
        

    elif request == 'gauss lam graph model' or request == '7':
        question_Gauss_1 = input('Handmade or auto?: ')
        if question_Gauss_1 == 'Handmade' or question_Gauss_1 == 'handmade':
            while True:
                print('According to graph, find ends and "h": ')

                lam_graph_model()
                plt.show()
                fig = plt.subplot()
                x = np.linspace(min(needed_loglam), max(needed_loglam), len(needed_flux))
                fig = plt.plot(10 ** x, needed_flux[range(len(needed_flux))] - needed_model[range(len(needed_model))])

                left_side = float(input('The end of left side is at (Angstroms): '))
                right_side = float(input('The end of right side is at (Angstroms): '))
                print('f(x) = a*e**(-((x-b)**2))/(2*(c**2)))')
                x_a = float(input('Coef a equals: '))
                x_b = float(input('Coef b equals: '))
                x_c = float(input('Coef c equals: '))
                h = float(input('h equals: '))
                x_a = x_a - h

                y = lambda x_1: x_a * np.exp(-(((x_1 - x_b) ** 2) / (2 * (x_c ** 2)))) + h

                x_1 = np.linspace(left_side, right_side, 1500)
                plt.plot(x_1, y(x_1))

                plt.show()

                question_Gauss_0 = input('Ok?: ')
                if question_Gauss_0 == 'no' or question_Gauss_0 == 'No':
                    pass
                elif question_Gauss_0 == 'yes' or question_Gauss_0 == 'Yes':
                    break

        elif question_Gauss_1 == 'Auto' or question_Gauss_1 == 'auto':
            while True:
                print('According to graph, find ends and "h": ')
                lam_graph_model()
                plt.show()
                fig = plt.subplot()
                x = np.linspace(min(needed_loglam), max(needed_loglam), len(needed_flux))
                fig = plt.plot(10 ** x, needed_flux[range(len(needed_flux))] - needed_model[range(len(needed_model))])

                left_side = float(input('Left end (Angstroms) equals: '))
                right_side = float(input('Right end (Angstroms) equals: '))
                h = float(input('h equals: '))

                p = ((10 ** needed_loglam) >= left_side) & ((10 ** needed_loglam) <= right_side)

                # print(p)
                needed_flux_new = (needed_flux - needed_model)[p]
                needed_loglam_new = 10 ** needed_loglam[p]

                # print(needed_flux_new)
                # print(needed_loglam_new)

                sum_nx = 0
                sum_n = 0

                for i in range(len(needed_loglam_new)):
                    sum_nx += needed_flux_new[i] * needed_loglam_new[i]
                    sum_n += needed_flux_new[i]
                average_x = sum_nx / sum_n
                sum_xx_squared = 0
                for i in range(len(needed_loglam_new)):
                    sum_xx_squared += needed_flux_new[i] * ((needed_loglam_new[i] - average_x) ** 2)
                sigma = sqrt(sum_xx_squared / sum_n)

                # print(sigma)
                x_a = max(needed_flux_new)
                x_c = sigma
                x_b = average_x
                # print(x_a, x_b, x_c, sep=' ')

                y = lambda x_1: (x_a - h) * np.exp(-(((x_1 - x_b) ** 2) / (2 * (x_c ** 2)))) + h

                x_1 = np.linspace(left_side, right_side, 1500)
                # print(x_1)
                # print(y(x_1))
                plt.plot(x_1, y(x_1))

                plt.show()

                question_Gauss_0 = input('Ok?: ')
                if question_Gauss_0 == 'no' or question_Gauss_0 == 'No':
                    pass
                elif question_Gauss_0 == 'yes' or question_Gauss_0 == 'Yes':
                    break
        

    elif request == 'lam graph and model' or request == '8':
        fig = plt.subplot()
        x = np.linspace(min(needed_loglam), max(needed_loglam), len(needed_flux))
        fig = plt.plot(10 ** x, needed_flux[range(len(needed_flux))])
        fig_2 = plt.plot(10 ** x, needed_model[range(len(needed_model))])
        plt.show()

    elif request == '':
        pass
    elif request == 'analyzed line' or request == '9':
        print(f'left side = {left_side}')
        print(f'right side = {right_side}')
        print(f'h = {h}')
        print(f'Sigma = {sigma}')
        print(f'a = {x_a}')
        print(f'b = {x_b}')
        
    elif request == 'histogram' or request == '10':
        x_histogram = []
        shoots = int(input('Input number of "shoots": '))
        for i in range(shoots):
            x_histogram.append(random.gauss(x_b, sigma))
        x_histogram.sort()
        x_histogram = np.array(x_histogram)
        lam_graph_model()
        plt.plot(x_histogram, gauss(x_histogram, x_a, x_b, sigma))
        plt.hist(x_histogram)
        plt.show()
        


        

    elif request == 'end' or request == '12':
        break
    else:
        continue

m_file.close()
