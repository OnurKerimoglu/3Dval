# Routine to read in WFD shapefile area means, and plot them
# Requires a shapefile with WFD polygons
#
# example calls:
# USE PYTHON 3.?
# with default arguments - defined below
# python showShapefile.py 
# with explicit arguments
# arguments explained:
# dataroot: root path for shapefile
# yearlist: list of years to plot (give in list format)
# scenarios: list of scenarios (give as list of strings). Scenario names must match those in shapefile!
# mos: mean or standard deviation. Leave empty for mean (or put dummy), put 'std' or '_std' for standard deviation.
# python showShapefile.py dataroot yearlist scenarios mos
# Author: Daniel Thewes, 25.7.22

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import geoplot
import mapclassify
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import math
import sys

# main routine - define default arguments below
def main(dataroot,yearlist,scenarios,lonlatextent,varns,mos):
    # define mos2 for title string
    mos2=mos.replace('_',' ')
    
    # To show the figure in python, set this switch to true:
    show_plot_switch=False
    
    # To save the figure, set this switch to true:
    print_plot_switch=True
    if not show_plot_switch:
        print_plot_switch=True # if it's not shown, it's saved

    # adapt figure format to your liking
    figuresize=(20,11)
    dpi=120
    
    # expand this if you want to plot any other variable. Modify if necessary.
    # To put proper variable name above each panel:
    prettytitle = {'Chl': 'Summer Chl-a', 'DIP': 'Winter DIP', 'DIN': 'Winter DIN'}
    # Color limits:
    vlimsNull=[0,0]
    if not mos=='_std':
        vlims = {'Chl': [0,25], 'DIN': [0,100], 'DIP': [0,2]}
    else:
        vlims = {'Chl': [0,8], 'DIN': [0,30], 'DIP': [0,0.5]}
        
    # Color limits for difference plots:
    # define explicit color limits in vlimsdiff, or set to [0,0], so that no explicit limits are applied
    vlimsdiff=vlimsNull 
    
    # Units (TEX format):
    units = {'Chl': 'mg~CHL~m^{-3}', 'DIN': 'mmol~N~m^{-3}', 'DIP': 'mmol~P~m^{-3}'}

    # define number of rows by number of scenarios: 
    # If there are more than two scenarios, add comparison line
    if len(scenarios)>1:
        nrows = len(scenarios)+1
    else:
        nrows = len(scenarios)
    
    # Choose which scenarios you want to compare:
    if len(scenarios) == 2:
        comparescenarios = [0,1]    
    else:
        comparescenarios = [len(scenarios)-1,len(scenarios)]
        
    if len(scenarios)>1:
        print(f'scenarios {scenarios[comparescenarios[0]]} and {scenarios[comparescenarios[1]]} will be compared')
        
    print(f'Years to plot: {yearlist}')
    for y in yearlist:
        infile = f'{dataroot}WFD_assessment_areameans_{y}.shp'
        
        # Reading files
        print('reading shapefile: '+infile)
        amshp = gpd.read_file(infile)
        
        #filtering out negative values (and NaN)
        amshp = amshp[amshp.DIN_CS>0]
        
        # name column
        idcol = "EU_CD_CW"
        areaids = list(amshp[idcol])

        # define figure
        f = plt.figure(figsize=figuresize, dpi=dpi)
        
        # Loop over scenarios
        for si,scen in enumerate(scenarios):
            # subplot index will be variable index + subplot index = vari+ssi:
            ssi=1+si*3 # new scenario: new line
                        
            # loop over variable
            for vari,varn in enumerate(varns):
                # generate panel
                panelindex = vari + ssi
                #ax1= f.add_subplot(nrows,3, vari+ssi,projection=ccrs.PlateCarree())

                # define title string
                titlestr=prettytitle[varn]+mos2+' '+scen+' [$'+units[varn]+'$]'
                
                # call panel plotting function
                plot_panel(f,nrows,panelindex,amshp,f'{varn}_{scen}{mos}',vlims[varn],lonlatextent,titlestr)

        # Plot differences between two scenarios
        # Define which scenarios you want to compare in variable comparescenarios, above.
        # Default are the last two scenarios in variable scenarios.
        if len(scenarios)>1:
            for vari,varn in enumerate(varns):
                # generate panel
                panelindex = vari + len(scenarios)*3 + 1
                #ax1 = f.add_subplot(nrows, 3, vari + 7, projection=ccrs.PlateCarree())
                
                # add column with differences
                tmpdiff=amshp[f'{varn}_{scenarios[comparescenarios[1]]}']-amshp[f'{varn}_{scenarios[comparescenarios[0]]}']
                columndiff=f'{varn}{mos}_diff'
                amshp[columndiff]=tmpdiff
                
                # define title string
                titlestr = prettytitle[varn] + f'{mos2} difference {scenarios[comparescenarios[1]]}-{scenarios[comparescenarios[0]]} [$' + units[varn] + '$]'
                                                                                                                                
                # call panel plotting function
                plot_panel(f,nrows,panelindex,amshp,columndiff,vlimsdiff,lonlatextent,titlestr)
            
        # plot cosmetics
        f.subplots_adjust(wspace=0.1, hspace=0.5)
        
        # output
        if print_plot_switch:
            if mos=='_std':
                plt.savefig(f'{dataroot}WFD_area_std.jpg', dpi = dpi)
                print(f'{dataroot}WFD_area_std.jpg saved!')
            else:
                plt.savefig(f'{dataroot}WFD_area_means.jpg', dpi = dpi)
                print(f'{dataroot}WFD_area_means.jpg saved!')
        if show_plot_switch:
            plt.show()

# Function to plot panels
def plot_panel(f,nrows,panelindex,amshp,column,vlims,lonlatextent,titlestr):
    ax = f.add_subplot(nrows, 3, panelindex, projection=ccrs.PlateCarree())
    if vlims[0]==vlims[1]:
        amshp.plot(column = column, legend = True, ax = ax)
    else:
        amshp.plot(column = column, legend = True, ax = ax, vmin=vlims[0], vmax=vlims[1])
                   
    # add coastlines etc.
    ax.add_feature(cfeature.OCEAN)
    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(cfeature.BORDERS)
    
    # set x and y axes limits
    ax.set_extent(lonlatextent, ccrs.PlateCarree())
    # define aspect ratio
    ax.set_aspect('auto')
    # set axis ticks
    ax.set_xticks(range(math.ceil(lonlatextent[0]),math.floor(lonlatextent[1])+1), crs=ccrs.PlateCarree())
    ax.set_yticks(range(math.ceil(lonlatextent[2]),math.floor(lonlatextent[3])+1), crs=ccrs.PlateCarree())
    
    # plot title
    plt.title(titlestr,size=12.0)
    
#default values
if __name__=='__main__':
    if len(sys.argv) > 1:
        dataroot = sys.argv[1]
    else:
        dataroot = '/home/daniel/levante_work/IR/Harmonization/WFD/'
        
    if len(sys.argv) > 2:
        yearlist = sys.argv[2]
    else:
        yearlist = [2017]
        
    if len(sys.argv) > 3:
        scenarios = sys.argv[3]
    else:
        scenarios = ['CS','28']
        
    # mean or standard deviation?
    # mos needs to be either 'std' or '_std'. Underscore will be added, if absent.
    if len(sys.argv) > 4:
        mos = sys.argv[4]
    else:
        mos = '' 
        
    if mos == 'std':
        mos = '_std'
    elif not mos == '_std':
        mos = ''
        
    # lonlatextent = [-0.2, 9.5, 51, 55.8]  
    lonlatextent = [4.5, 9.5, 52.8, 55.8] # focus on Wadden Sea
    varns = ['Chl','DIN','DIP']
    
    # call main function
    main(dataroot,yearlist,scenarios,lonlatextent,varns,mos)
