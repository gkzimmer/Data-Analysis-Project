#This is the main python file for reading in, performing calculations on, and graphing the spectral data
#of red giant field stars.
#We are focused on three heavy, n-capture elements: Barium, Europium, and Lanthanum.
#The data is a mixture of individual data (analyzed through MOOG) and external data from published papers

import numpy as np
import math
import matplotlib.pyplot as plt

#First, read in the data.
spectral_data = np.genfromtxt('spectra.txt', dtype=None, skip_header=1)

#Now, we should convert the LogE units to x/Fe (Note: we could do the other way around, but since most of the
#literature is x/Fe, I will use that conversion)
la_fe_solar = -10.73
ba_fe_solar = -9.87
eu_fe_solar = -11.82
x_marker_metal = 0

#Get all the abundances with LogE units
for n in range(spectral_data.shape[0]):

    if spectral_data[n][1] == 'logE':

        if spectral_data[n][5] != 'x':
            #Now, find out which elements to convert- 'x' notates that we do not have the abundance for that element
            if spectral_data[n][2] != 'x':
                spectral_data[n][2] = format((float(spectral_data[n][2]) - 12.0 - float(spectral_data[n][5]) - ba_fe_solar), 'f')
            if spectral_data[n][3] != 'x':
                spectral_data[n][3] = format((float(spectral_data[n][3]) - 12.0 - float(spectral_data[n][5]) - eu_fe_solar), 'f')
            if spectral_data[n][4] != 'x':
                spectral_data[n][4] = format((float(spectral_data[n][4]) - 12.0 - float(spectral_data[n][5]) - la_fe_solar), 'f')

        if spectral_data[n][5] == 'x':
            x_marker_metal = x_marker_metal + 1
            
#First, get rid of all stars without metallicity abundances
for n in range(spectral_data.shape[0] - x_marker_metal):
    if spectral_data[n][5] == 'x':
        spectral_data = np.delete(spectral_data, n, 0)

#This is a function for deleting lines with 'x', which represents no abundance for the selected element
#Some stars have only one or two of the element abundances, so the number needed to graph elements against one another
#is different for every graph.
def delete_x(num, place1, place2, array, char):
    spectra_list_x = []
    spectra_list_y = []
    for n in range(num):
        if array[n][place1] != 'x' and array[n][place2] != 'x' and array[n][6] != char:
            spectra_list_x.append((array[n][place1]))
            spectra_list_y.append((array[n][place2]))
    spectra_array = np.array([spectra_list_x, spectra_list_y])
    spectra_array = spectra_array.astype(np.float)
    return spectra_array

#Function for calculating the best line fit for the graphs
#Based on exercise 3.8 in CP
def line_fit(num_points, x, y):
    e_x = (np.sum(x) / num_points)
    e_y = (np.sum(y) / num_points)
    e_xx = (np.sum(x**2) / num_points)
    e_xy = (np.sum(x * y) / num_points)
    m = (e_xy - (e_x * e_y)) / (e_xx - (e_x**2))
    c = ((e_xx * e_y) - (e_x * e_xy)) / (e_xx - (e_x**2))

    print(num_points)
    line_fit = np.empty([num_points])
    for n in range(num_points):
        line_fit[n] = m * x[n] + c

    return line_fit

#Get all our data organized, with no x values

#All our literature values
spectral_ba_metal_lit = delete_x(spectral_data.shape[0],2,5,spectral_data, 'I')
spectral_eu_metal_lit = delete_x(spectral_data.shape[0],3,5,spectral_data, 'I')
spectral_la_metal_lit = delete_x(spectral_data.shape[0],4,5,spectral_data, 'I')
spectral_ba_eu_lit = delete_x(spectral_data.shape[0],2,3,spectral_data, 'I')
spectral_ba_la_lit = delete_x(spectral_data.shape[0],2,4,spectral_data, 'I')
spectral_eu_la_lit = delete_x(spectral_data.shape[0],3,4,spectral_data, 'I')

#All the individually-sourced (through MOOG spectral analysis) data
spectral_ba_metal_ind = delete_x(spectral_data.shape[0],2,5,spectral_data, 'E')
spectral_eu_metal_ind = delete_x(spectral_data.shape[0],3,5,spectral_data, 'E')
spectral_la_metal_ind = delete_x(spectral_data.shape[0],4,5,spectral_data, 'E')
spectral_ba_eu_ind = delete_x(spectral_data.shape[0],2,3,spectral_data, 'E')
spectral_ba_la_ind = delete_x(spectral_data.shape[0],2,4,spectral_data, 'E')
spectral_eu_la_ind = delete_x(spectral_data.shape[0],3,4,spectral_data, 'E')

#All data altogether, for the line fit calculations
spectral_ba_metal = delete_x(spectral_data.shape[0],2,5,spectral_data, '')
spectral_eu_metal = delete_x(spectral_data.shape[0],3,5,spectral_data, '')
spectral_la_metal = delete_x(spectral_data.shape[0],4,5,spectral_data, '')
spectral_ba_eu = delete_x(spectral_data.shape[0],2,3,spectral_data, '')
spectral_ba_la = delete_x(spectral_data.shape[0],2,4,spectral_data, '')
spectral_eu_la = delete_x(spectral_data.shape[0],3,4,spectral_data, '')


ba_metal_linefit = line_fit(spectral_ba_metal.shape[1],spectral_ba_metal[1,:],spectral_ba_metal[0,:])
eu_metal_linefit = line_fit(spectral_eu_metal.shape[1],spectral_eu_metal[1,:],spectral_eu_metal[0,:])
la_metal_linefit = line_fit(spectral_la_metal.shape[1],spectral_la_metal[1,:],spectral_la_metal[0,:])
ba_eu_linefit = line_fit(spectral_ba_eu.shape[1],spectral_ba_eu[0,:],spectral_ba_eu[1,:])
ba_la_linefit = line_fit(spectral_ba_la.shape[1],spectral_ba_la[0,:],spectral_ba_la[1,:])
eu_la_linefit = line_fit(spectral_eu_la.shape[1],spectral_eu_la[0,:],spectral_eu_la[1,:])

#Seperate into invidiual (MOOG) data, and data we gathered from literature searches

#Stick to blue/orange/red for colorblind individuals
plt.scatter(spectral_ba_eu_ind[0,:],spectral_ba_eu_ind[1,:],marker='s',label='Individual Data Points',color='tab:orange')
plt.scatter(spectral_ba_eu_lit[0,:],spectral_ba_eu_lit[1,:],marker='o',label='Literature Data Points',color='b')
plt.plot(spectral_ba_eu[0,:],ba_eu_linefit,label='Line Fit',color='r')
plt.xlabel('Barium Abundances (Ba/Fe)',fontsize=15)
plt.ylabel('Europium Abundances (Eu/Fe)',fontsize=15)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.title('Europium v Barium Abundances of Red Giant Field Stars',fontsize=15)
plt.savefig('Ba_v_Eu.png')
plt.clf()

plt.scatter(spectral_ba_la_ind[0,:],spectral_ba_la_ind[1,:],marker='s',label='Individual Data Points',color='tab:orange')
plt.scatter(spectral_ba_la_lit[0,:],spectral_ba_la_lit[1,:],marker='o',label='Literature Data Points',color='b')
plt.plot(spectral_ba_la[0,:],ba_la_linefit,label='Line Fit',color='r')
plt.xlabel('Barium Abundances (Ba/Fe)',fontsize=15)
plt.ylabel('Lanthanum Abundaces (La/Fe)',fontsize=15)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.title('Barium v Lanthanum Abundances of Red Giant Field Stars',fontsize=15)
plt.savefig('Ba_v_La.png')
plt.clf()

plt.scatter(spectral_eu_la_ind[0,:],spectral_eu_la_ind[1,:],marker='s',label='Individual Data Points',color='tab:orange')
plt.scatter(spectral_eu_la_lit[0,:],spectral_eu_la_lit[1,:],marker='o',label='Literature Data Points',color='b')
plt.plot(spectral_eu_la[0,:],eu_la_linefit,label='Line Fit',color='r')
plt.xlabel('Europium Abundances (Eu/Fe)',fontsize=15)
plt.ylabel('Lanthanum Abundances (La/Fe)',fontsize=15)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.title('Europium v Lanthanum Abundances of Red Giant Field Stars',fontsize=15)
plt.savefig('Eu_v_La.png')
plt.clf()

plt.scatter(spectral_ba_metal_ind[1,:],spectral_ba_metal_ind[0,:],marker='s',label='Individual Data Points',color='tab:orange')
plt.scatter(spectral_ba_metal_lit[1,:],spectral_ba_metal_lit[0,:],marker='o',label='Literature Data Points',color='b')
plt.plot(spectral_ba_metal[1,:],ba_metal_linefit,label='Line Fit',color='r')
plt.ylabel('Barium Abundances (Ba/Fe)',fontsize=15)
plt.xlabel('Metallicity (Fe/H)',fontsize=15)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.title('Barium Abundances v Metallicity of Red Giant Field Stars',fontsize=15)
plt.savefig('Ba_v_Metal.png')
plt.clf()

plt.scatter(spectral_eu_metal_ind[1,:],spectral_eu_metal_ind[0,:],marker='s',label='Individual Data Points',color='tab:orange')
plt.scatter(spectral_eu_metal_lit[1,:],spectral_eu_metal_lit[0,:],marker='o',label='Literature Data Points',color='b')
plt.plot(spectral_eu_metal[1,:],eu_metal_linefit,label='Line Fit',color='r')
plt.ylabel('Europium Abundances (Eu/Fe)',fontsize=15)
plt.xlabel('Metallicity (Fe/H)',fontsize=15)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.title('Europium Abundances v Metallicity of Red Giant Field Stars',fontsize=15)
plt.savefig('Eu_v_Metal.png')
plt.clf()

plt.scatter(spectral_la_metal_ind[1,:],spectral_la_metal_ind[0,:],marker='s',label='Individual Data Points',color='tab:orange')
plt.scatter(spectral_la_metal_lit[1,:],spectral_la_metal_lit[0,:],marker='o',label='Literature Data Points',color='b')
plt.plot(spectral_la_metal[1,:],la_metal_linefit,label='Line Fit',color='r')
plt.ylabel('Lanthanum Abundances (La/Fe)',fontsize=15)
plt.xlabel('Metallicity (Fe/H)',fontsize=15)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.title('Lanthanum Abundances v Metallicity of Red Giant Field Stars',fontsize=15)
plt.savefig('La_v_Metal.png')
