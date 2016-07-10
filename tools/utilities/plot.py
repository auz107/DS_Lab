from __future__ import division
import sys, time
sys.path.append('../')
import os
import matplotlib
from imp import load_source
from matplotlib import colors 
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from tools.userError import userError
import numpy as np

#------------------------------------------
# Ali R. Zomorrodi - Segre Lab @ BU
# Last updated: 06-15-2016
#------------------------------------------

"""
TIPS:

1. The following line makes the math font regular (the same size the normal font)
Source: http://stackoverflow.com/questions/15836050/matplotlib-change-math-font-size 
matplotlib.rcParams['mathtext.default']='regular'
This is not easy to set back to default. See:
http://stackoverflow.com/questions/22569071/matplotlib-change-math-font-size-and-then-go-back-to-default 
Instead you can consider using \mathdefault{...} or its alias \mathregular{...}
Source: http://matplotlib.org/users/mathtext.html#fonts

3. To make part of mathtext bold use \bf{} or \mathbf
3. To remove italic in part of a mathtext use \rm{} 
"""

#--- Define some global default values ---
# Defaulty figure font name (set by rcParam)
figure_default_fontname = 'serif'

# Defaulty figure font size (set by rcParam)
figure_default_fontsize = 20

# Defaulty figure font weight  
figure_default_fontweight = 'regular'

# Default font size for title
title_default_fontname = 'Arial'
title_default_fontsize = 24

# Defaul font size for axes and colorbar labels
axes_label_default_fontname = 'Arial'
axes_label_default_fontsize = 23

# Default font name and size for tick labels 
ticklabels_default_fontname = None
ticklabels_default_fontsize = 22

# Default font name and size for legend 
legend_default_fontname = None
legend_default_fontsize = 23

# Defaulty figure size
figsize_default = None 

# default colormap
#colormap_default = matplotlib.cm.RdBu
colormap_default = None 

class axis(object):
    """
    A class holding the axes properties
    """
    def __init__(self,label = '', label_format = {'fontname': axes_label_default_fontname,'fontweight':figure_default_fontweight,'fontsize':axes_label_default_fontsize, 'distance_from_ticklabels':None, 'position': None}, limits = None, set_minorticks = False, minorticks_spacing = None, majorticks_spacing = None, custom_ticks = None, custom_ticklabels = None, ticklabels_format = {'fontname':ticklabels_default_fontname, 'fontweight':figure_default_fontweight,'fontsize':ticklabels_default_fontsize,'rotation':0, 'string_format':None, 'position':None,'horizontalalignment': None, 'verticalalignment':None}, plot_gridlines = False, gridlines_format = {'color':'black', 'linestyle':'dashed', 'linewidth':1}, spines_format = {}, invert = False, scale = None):
        """
        INPUTS:
        -------
        label: 
        Axis label (a string)

        label_format: 
        A dictionary showing the format of the axis label

        limits: 
        Axis limits. A tuple in the form of (min,max)

        set_minorticks: 
        Showing whether to set the minor ticks (True) or nor (False)

        minorticks_spacing: 
        Space between minor ticks

        majorticks_spacing: 
        Space between major ticks

         custom_ticks: 
         Custom ticks

         custom_ticklabels: 
         Custom tick labels

         ticklabels_format: 
         Format of the tick labels

         plot_gridlines: 
         Showing whether to plot the grid lines (True) or not (False)
            gridlines_fomrat: A dictionary containing the gridlines properties
                 show_legend: Whether to show the legend (True) or not (False)
               legend_format: A dictionary containing the legend properties

        spines_format: 
        Format of the axis lines. This a dictionary with allowed keys being 'top', 'bottom', 
        'left', 'right'. The values must be another dictionary with either of these three keys: 
       'linewidth', 'linestyle' and 'linecolor'

        invert: 
        Shows whether to invert an axis (True) or not (False)

        scale: 
        The scale of the axis. Allowed choices are 'linear', 'log', 'logit', 'symlog'
        """
        # Axis label
        self.label = label

        # label format
        self.label_format = label_format
        if 'fontname' not in self.label_format.keys():
            self.label_format['fontname'] = axes_label_default_fontname 
        if 'fontsize' not in self.label_format.keys():
            self.label_format['fontsize'] = axes_label_default_fontsize 
        if 'fontweight' not in self.label_format.keys():
            self.label_format['fontweight'] = figure_default_fontweight 
        if 'rotation' not in self.label_format.keys():
            self.label_format['rotation'] = 0 
        if 'distance_from_ticklabels' not in self.label_format.keys():
            self.label_format['distance_from_ticklabels'] = None 
        if 'position' not in self.label_format.keys():
            self.label_format['position'] = None 

        # ticklabels_format
        self.ticklabels_format = ticklabels_format 
        if 'fontname' not in self.ticklabels_format.keys():
            self.ticklabels_format['fontname'] = ticklabels_default_fontname
        if 'fontsize' not in self.ticklabels_format.keys():
            self.ticklabels_format['fontsize'] = ticklabels_default_fontsize
        if 'fontweight' not in self.ticklabels_format.keys():
            self.ticklabels_format['fontweight'] = figure_default_fontweight 
        if 'rotation' not in self.ticklabels_format.keys():
            self.ticklabels_format['rotation'] = 0 
        if 'string_format' not in self.ticklabels_format.keys():  # string_format is like %1.2f or %d
            self.ticklabels_format['string_format'] = None 
        if 'position' not in self.ticklabels_format.keys():  
            self.ticklabels_format['position'] = None 
        if 'horizontalalignment' not in self.ticklabels_format.keys(): 
            self.ticklabels_format['horizontalalignment'] = None 
        if 'verticalalignment' not in self.ticklabels_format.keys(): 
            self.ticklabels_format['verticalalignment'] = None 

        # Axis limits (a tuple with the first and second elements being the min and max)
        self.limits = limits

        # Set minor ticks (True or False)
        self.set_minorticks = set_minorticks

        # Specify minor tick spacing
        self.minorticks_spacing = minorticks_spacing

        # Specify major tack spacing
        self.majorticks_spacing = majorticks_spacing

        # Custom ticks (must be an array of integers or float) 
        self.custom_ticks = custom_ticks

        # Custom tick labels (must be an array of strings) 
        self.custom_ticklabels = custom_ticklabels

        # Grid lines
        self.plot_gridlines = plot_gridlines 
    
        # gridlines format
        if self.plot_gridlines:
            self.gridlines_format = gridlines_format
            gridlines_format_keys = [k.lower() for k in gridlines_format.keys()]
            if 'color' not in gridlines_format_keys:
                self.gridlines_format['color'] = 'k' 
            if 'linestyle' not in gridlines_format_keys:
                self.gridlines_format['linestyle'] = 'dashed' 
            if 'linewidth' not in gridlines_format_keys:
                self.gridlines_format['linewidth'] = 2 

        # Spines (axes and graph's border lines) format
        # Here main and opposite are the main axis spine and the one opoosite to it. For example, for 
        # an x axis main may refer to the one in the bottom and opposite to that on the top. Similarly for a y
        # axis, main and opoosite refer to the left and right spines. The definision, which one is main and which one
        # is the opposite is up to the user
        # Source: http://matplotlib.org/api/spines_api.html#matplotlib.spines.Spine.set_position
        self.spines_format = spines_format

        # Whether to use inverted axis
        self.invert = invert

        # Scale. Acceptable choices are 'linear', 'log', 'logit', 'symlog'
        self.scale = scale

    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
        -------
        attr_name: Attribute name
        attr_value: Attribute value
        """
        if attr_name.lower() == 'label' and not isinstance(attr_value,str):
            raise TypeError('label must be a string')
 
        if attr_name.lower() == 'limits' and attr_value != None and not isinstance(attr_value,tuple):
            raise TypeError('limits must be a tuple')
        else:
            if attr_name.lower() == 'limits' and attr_value != None and len(attr_value) != 2:
                raise ValueError('limits must be a tuple in the form of (min_value,max_value)')
            if attr_name.lower() == 'limits' and attr_value != None and attr_value[0] != None and attr_value[1] != None and attr_value[0] >= attr_value[1]:
                raise ValueError('The second entry in limits must be greater than the first one. The entered value for limits is: {}'.format(attr_value))
 
        if attr_name.lower() == 'set_minorticks' and not isinstance(attr_value,bool):
            raise TypeError('set_minorticks must be either True or False')
 
        if attr_name.lower() == 'set_majorticks' and not isinstance(attr_value,bool):
            raise TypeError('set_majorticks must be either True or False')
         
        if attr_name.lower() == 'label_format' and not isinstance(attr_value,dict):
            raise TypeError('label_format must be a dictionary')
        elif attr_name.lower() == 'label_format' and len([k.lower() for k in attr_value.keys() if k not in ['fontname','fontweight','fontsize','rotation','distance_from_ticklabels','position']]) > 0:
            raise ValueError('Unknown key(s) for label_format: {}'.format([k for k in attr_value.keys() if k.lower() not in ['fontname','fontweight','fontsize','distance_from_ticklabels','position']]))
 
        if attr_name.lower() == 'ticklabels_format' and not isinstance(attr_value,dict):
            raise TypeError('ticklabels_format must be a dictionary')
        elif attr_name.lower() == 'ticklabels_format' and len([k for k in attr_value.keys() if k.lower() not in ['fontname','fontweight','fontsize','rotation','string_format','position','horizontalalignment','verticalalignment','distance_from_ticklabels']]) > 0: 
            raise ValueError('Unknown key(s) for label_format: {}'.format([k for k in attr_value.keys() if k.lower() not in ['fontname','fontweight','fontsize','rotation','string_format','position','horizontalalignment','verticalalignment','distance_from_ticklabels']]))
 
        if attr_name.lower() == 'custom_ticks' and attr_value != None and not isinstance(attr_value,list):
            raise TypeError('custom_ticks must be an array of integers or floats')
        elif attr_name.lower() == 'custom_ticks' and attr_value != None and len([k for k in attr_value if not isinstance(k,int) and not isinstance(k,float)]): 
            raise ValueError('custom_ticks must be an array of integers or floats. Non-integer and non-float values were entered: {}'.format([k for k in attr_value if not isinstance(k,int) and not isinstance(k,float)]))
 
        if attr_name.lower() == 'custom_ticklabels' and attr_value != None and not isinstance(attr_value,list):
            raise TypeError('custom_ticklabels must be an array of strings')
        elif attr_name.lower() == 'custom_ticklabels' and attr_value != None and len([k for k in attr_value if not isinstance(k,str)]) > 0: 
            raise ValueError('custom_ticklabels must be an array of strings. Non-string values were entered: {}'.format([k for k in attr_value if not isinstance(k,str)]))
 
        if attr_name.lower() == 'invert' and not isinstance(attr_value,bool):
            raise TypeError('invert must be either True or False')
 
        if attr_name.lower() == 'scale' and attr_value != None and not isinstance(attr_value,str):
            raise TypeError('scale must be a string')
        elif attr_name.lower() == 'scale' and attr_value != None and attr_value not in ['linear', 'log', 'logit', 'symlog']: 
           raise ValueError('Invalid value for scale! Allowed choices are: linear, log, logit, symlog')

        # Gridlines
        if attr_name.lower() == 'plot_gridlines' and not isinstance(attr_value,bool):
            raise TypeError('plot_gridlines must be either True or False')
 
        if attr_name.lower() == 'gridlines_format' and not isinstance(attr_value,dict):
            raise TypeError('gridlines_format must be a dictionary')
        elif attr_name.lower() == 'gridlines_format' and len([k for k in attr_value.keys() if k.lower() not in ['color', 'linestyle', 'linewidth']]) > 0: 
            raise ValueError('Unknown key(s) for label_format: {}'.format([k for k in attr_value.keys() if k.lower() not in ['color', 'linestyle', 'linewidth']]))


        # Spines (axes and graph's border lines)
        if attr_name.lower() == 'spines_format':
            if  not isinstance(attr_value,dict):
                raise TypeError('spines_format must be a dictionary')
            if len([k for k in attr_value.keys() if k.lower() not in ['top','bottom','left','right']]) > 0:
                raise ValueError('Invalid key for spines_format: {}. Allowed keys are [top, bottom, left, right]'.format([k for k in attr_value.keys() if k.lower() not in ['top','bottom', 'left', 'right']]))

            for spines_key in attr_value.keys():
                if not isinstance(attr_value[spines_key],dict):
                    raise TypeError("spines_format['{}'] must be a dictionary".format(spines_key))
                if len([k for k in attr_value[spines_key].keys() if k.lower() not in ['linewidth','linestyle','linecolor']]) > 0:
                    raise ValueError("Invalid key for spines_format['{}']: {}. Allowed keys are [linestyle,linewidth,linecolor]".format(spines_key,[k for k in attr_value[spines_key].keys() if k.lower() not in ['linestyle','linewidth','linecolor']]))


        self.__dict__[attr_name] = attr_value

class color_bar(object):
    """
    A class holding the properties of the color bar (for related graphs)
    """
    def __init__(self, colormap = None, colorlimits = None, label = '', label_format = {'fontweight':figure_default_fontweight,'fontsize':axes_label_default_fontsize,'rotation':270,'distance_from_ticklabels':None, 'string_format': None}, set_minorticks = False, minorticks_spacing = None, majorticks_spacing = None, custom_ticks = None, custom_ticklabels = None, ticklabels_format = {'fontname':ticklabels_default_fontname, 'fontweight':figure_default_fontweight,'fontsize':ticklabels_default_fontsize,'rotation':0, 'string_format':None, 'axis_position':None, 'ticks_position':None, 'horizontalalignment': None, 'verticalalignment':None}):

        # colorbar
        if colormap != None:
            self.colormap = colormap
        else:
            self.colormap = colormap_default
        
        # colorlimits: A tuple in the form (min,max) indicating the Color range for the colormap. These values are 
        # used as vmin and vmax for pcolor, pcolormesh, matshow or imshow.  
        # As a general rule for pcolor and pcolormesh mesh vmin = data.min() and vmax = data.max(). 
        # If you like to have ticklabels placed in the middle of ticks in (1) a discerete colormap for pcolor 
        # and pcolormesh, or (2) in general for matshow and imshow use: 
        # vmin = data.min() - 0.5 and vmax = data.max() + 0.5. 
        self.colorlimits = colorlimits

        # Label
        self.label = label

        # label_format
        self.label_format = label_format
        if 'fontsize' not in self.label_format.keys():
            self.label_format['fontsize'] = axes_label_default_fontsize 
        if 'fontweight' not in self.label_format.keys():
            self.label_format['fontweight'] = figure_default_fontweight 
        if 'rotation' not in self.label_format.keys():
            self.label_format['rotation'] = 270 
        if 'distance_from_ticklabels' not in self.label_format.keys():
            self.label_format['distance_from_ticklabels'] = None 
        if 'string_format' not in self.label_format.keys():
            self.label_format['string_format'] = None 

        # Set minor ticks (True or False)
        self.set_minorticks = set_minorticks

        # Specify minor tick spacing
        self.minorticks_spacing = minorticks_spacing

        # Specify major tack spacing
        self.majorticks_spacing = majorticks_spacing

        # Custom ticks (must be an array of integers or float) 
        self.custom_ticks = custom_ticks

        # Custom tick labels (must be an array of strings) 
        self.custom_ticklabels = custom_ticklabels

        # ticklabels_format
        self.ticklabels_format = ticklabels_format 
        if 'fontname' not in self.ticklabels_format.keys():
            self.ticklabels_format['fontname'] = ticklabels_default_fontname
        if 'fontsize' not in self.ticklabels_format.keys():
            self.ticklabels_format['fontsize'] = ticklabels_default_fontsize
        if 'fontweight' not in self.ticklabels_format.keys():
            self.ticklabels_format['fontweight'] = figure_default_fontweight 
        if 'rotation' not in self.ticklabels_format.keys():
            self.ticklabels_format['rotation'] = 0 
        if 'string_format' not in self.ticklabels_format.keys():  # string_format is like %1.2f or %d
            self.ticklabels_format['string_format'] = None 
        if 'axis_position' not in self.ticklabels_format.keys(): # Tick labels on top/bottom/left/right acis 
            self.ticklabels_format['axis_position'] = None 
        if 'ticks_position' not in self.ticklabels_format.keys(): # Tick labels position w.r.t. to ticks. The only aalowed choice is 'middle' 
            self.ticklabels_format['ticks_position'] = None 
        if 'horizontalalignment' not in self.ticklabels_format.keys(): 
            self.ticklabels_format['horizontalalignment'] = None 
        if 'verticalalignment' not in self.ticklabels_format.keys(): 
            self.ticklabels_format['verticalalignment'] = None 


    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
        -------
        attr_name: Attribute name
        attr_value: Attribute value
        """
        if attr_name.lower() == 'label' and not isinstance(attr_value,str):
            raise TypeError('label must be a string')

        if attr_name.lower() == 'colorlimits' and attr_value != None and not isinstance(attr_value,tuple):
            raise TypeError('colorlimits must be a tuple of form (min,max)')
        if attr_name.lower() == 'colorlimits' and attr_value != None and len(attr_value) != 2:
            raise TypeError('colorlimits must be a tuple of size two (min,max)')
        elif attr_name.lower() == 'colorlimits' and attr_value != None and not attr_value[1] > attr_value[0]:
            raise TypeError('min value in colorlimits is greater than its max: {}'.format(attr_value))

        if attr_name.lower() == 'label_format' and not isinstance(attr_value,dict):
            raise TypeError('label_format must be a dictionary')
        elif attr_name.lower() == 'label_format' and len([k for k in attr_value.keys() if k.lower() not in ['fontweight','fontsize','rotation','distance_from_ticklabels','string_format']]) > 0: 
            raise ValueError('Unknown key(s) for label_format: {}'.format([k for k in attr_value.keys() if k.lower() not in ['fontweight','fontsize','rotation','distance_from_ticklabels','string_format']]))

        if attr_name.lower() == 'ticklabels_format' and not isinstance(attr_value,dict):
            raise TypeError('ticklabels_format must be a dictionary')
        elif attr_name.lower() == 'ticklabels_format' and len([k for k in attr_value.keys() if k.lower() not in ['fontname','fontweight','fontsize','rotation','string_format','axis_position', 'ticks_position','horizontalalignment','verticalalignment','distance_from_ticklabels']]) > 0: 
            raise ValueError('Unknown key(s) for label_format: {}'.format([k for k in attr_value.keys() if k.lower() not in ['fontname','fontweight','fontsize','rotation','string_format','axis_position', 'ticks_position','horizontalalignment','verticalalignment','distance_from_ticklabels']]))

        if attr_name.lower() == 'set_minorticks' and not isinstance(attr_value,bool):
            raise TypeError('set_minorticks must be either True or False')
 
        if attr_name.lower() == 'set_majorticks' and not isinstance(attr_value,bool):
            raise TypeError('set_majorticks must be either True or False')
 
        if attr_name.lower() == 'custom_ticks' and attr_value != None and not isinstance(attr_value,list) and not isinstance(attr_value, np.ndarray):
            raise TypeError('custom_ticks must be an array of integers or floats or a numpy ndarray')
        elif attr_name.lower() == 'custom_ticks' and attr_value != None and len([k for k in attr_value if not isinstance(k,int) and not isinstance(k,float)]): 
            raise ValueError('custom_ticks must be an array of integers or floats or a numpy ndarray. Non-integer and non-float values were entered: {}'.format([k for k in attr_value if not isinstance(k,int) and not isinstance(k,float)]))
 
        if attr_name.lower() == 'custom_ticklabels' and attr_value != None and not isinstance(attr_value,list):
            raise TypeError('custom_ticklabels must be an array of strings')
        elif attr_name.lower() == 'custom_ticklabels' and attr_value != None and len([k for k in attr_value if not isinstance(k,str)]) > 0: 
            raise ValueError('custom_ticklabels must be an array of strings. Non-string values were entered: {}'.format([k for k in attr_value if not isinstance(k,str)]))

        self.__dict__[attr_name] = attr_value

class plot(object):
    """
    Makes different types of plots including
    2D plot
    3D plot
    Heatmap plots
    """
    def __init__(self, title = '', title_format = {'fontname':title_default_fontname,'fontweight':figure_default_fontweight,'fontsize':title_default_fontsize, 'distance_from_graph':1.05}, xaxis = None, yaxis = None,zaxis = None, plot_gridlines = False, gridlines_format = {'color':'black', 'linestyle':'dashed', 'linewidth':1}, show_legend = False, legend_format = {'location':0, 'fontname':legend_default_fontname,'fontsize': figure_default_fontsize, 'fontweight':figure_default_fontweight, 'linewidth': 2, 'bbox_to_anchor': None}, fig_format = {'figsize':figsize_default, 'fontname': figure_default_fontname, 'fontsize':figure_default_fontsize,'fontweight':figure_default_fontweight,'use_tight_layout':True, 'dpi':None, 'use_LaTex': False, 'mathtext_fontname':None}, output_filename = ''):
        """
        Constructor takes only inputs that are general to all types of plots
        INPUTS:
        -------
                   title: Title of the figure
            title_format: A dictionary containing the title format properties
                   xaxis: An instance of class axis
                   yaxis: An instance of class axis
                   zaxis: An instance of class axis
          plot_gridlines: Showing whether to plot the grid lines (True) or not (False). If True,
                          gridlines are plotted for ALL axes by setting their plot_gridlines to True.
                          To plot gridlines only for specific axes, set this parameter to False here,
                          but set it to True for the desired axis.
        gridlines_fomrat: A dictionary containing the gridlines properties. This format is used for 
                          ALL axes. To choose different formats for different axes, set this parameter
                          for each axis separately (see class axis for details)
             show_legend: Whether to show the legend (True) or not (False)
           legend_format: A dictionary containing the legend properties
              fig_format: Figure format. Keys include:
                                   figsize: A tuple showing the dimensions of the figure
                                  fontname: Font name (font family) of the figure
                                  fontsize: Font size of the figure
                          use_tight_layout: Showing whether to to use a tight layout (True) or not (False)
                                       dpi: dpi of the output figure
                                 use_LaTex: use Latex to render the text elements in the figure
                         mathtext_fontname: Custom fontname for mathtext (build-in matplotlib LaTex)
         output_filename: Name of the output file to save the figure
    
        Ali R. Zomorrodi, Segre Lab @ BU
        Last updated: 02-29-2016
        """
        # Title
        self.title = title
    
        # Title format
        self.title_format = title_format
        title_format_keys = [k.lower() for k in title_format.keys()]
        if 'fontname' not in title_format_keys: 
            self.title_format['fontname'] = title_default_fontname 
        if 'fontsize' not in title_format_keys: 
            self.title_format['fontsize'] = title_default_fontsize 
        if 'fontweight' not in title_format_keys: 
            self.title_format['fontweight'] = figure_default_fontweight 
    
        # x axes
        if xaxis == None:
            self.xaxis = axis
        else:
            self.xaxis = xaxis
    
        # y axes
        if yaxis == None:
            self.yaxis = axis
        else:
            self.yaxis = yaxis

        # z axis  
        if zaxis != None:
            self.zaxis = zaxis
    
        # Grid lines and their format
        if plot_gridlines:
            gridlines_format_keys = [k.lower() for k in gridlines_format.keys()]
            if 'color' not in gridlines_format_keys:
                gridlines_format['color'] = 'k' 
            if 'linestyle' not in gridlines_format_keys:
                gridlines_format['linestyle'] = 'dashed' 
            if 'linewidth' not in gridlines_format_keys:
                gridlines_format['linewidth'] = 2 

            self.xaxis.plot_gridlines = True
            self.xaxis.gridlines_format = gridlines_format

            self.yaxis.plot_gridlines = True
            self.yaxis.gridlines_format = gridlines_format

            if zaxis != None:
                self.zaxis.plot_gridlines = True
                self.zaxis.gridlines_format = gridlines_format

        # Show legend
        self.show_legend = show_legend

        # Legend properties
        self.legend_format = legend_format
        # See http://matplotlib.org/api/legend_api.html for help on location
        legend_format_keys = [k.lower() for k in self.legend_format.keys()]
        if self.legend_format != None and 'location' not in legend_format_keys:
            self.legend_format['location'] = 'best'
        if self.legend_format != None and 'fontname' not in legend_format_keys:
            self.legend_format['fontname'] = legend_default_fontname
        if self.legend_format != None and 'fontsize' not in legend_format_keys:
            self.legend_format['fontsize'] = figure_default_fontsize
        if self.legend_format != None and 'fontweight' not in legend_format_keys:
            self.legend_format['fontweight'] = figure_default_fontweight 
        if self.legend_format != None and 'linewidth' not in legend_format_keys:
            self.legend_format['linewidth'] = 2 
        if self.legend_format != None and 'bbox_to_anchor' not in legend_format_keys:
            self.legend_format['bbox_to_anchor'] = None 
    
        # Figure format
        self.fig_format = fig_format
        fig_format_keys = [k.lower() for k in self.fig_format.keys()] 
        if 'figsize' not in fig_format_keys:
            self.fig_format['figsize'] = figsize_default
        if 'fontname' not in fig_format_keys:
            self.fig_format['fontname'] = figure_default_fontname
        if 'fontsize' not in fig_format_keys:
            self.fig_format['fontsize'] = figure_default_fontsize
        if 'fontweight' not in fig_format_keys:
            self.fig_format['fontweight'] = figure_default_fontweight 
        if 'use_tight_layout' not in fig_format_keys:
            self.fig_format['use_tight_layout'] = True
        if 'dpi' not in fig_format_keys: 
            self.fig_format['dpi'] = None
        if 'use_latex' not in fig_format_keys: 
            self.fig_format['use_LaTex'] = False
        if 'mathtext_fontname' not in fig_format_keys: 
            self.fig_format['mathtext_fontname'] = None
        
        # File name storing the figure
        self.output_filename = output_filename
    
        # A parameter showing the type of plot ('heatmap', '2d_line', '3d'), The value is assigned inside each function 
        self._plot_type = None
    
        # Update the matplotlib configuration parameters:
        # If use_LaTex - True and the fontweight is bold then one needs to load amsmath into the TeX preamble to enable
        # bold LaTex fonts, otherwise the fonts will not be bold.
        # Source: http://stackoverflow.com/questions/14324477/bold-font-weight-for-latex-axes-label-in-matplotlib
        if self.fig_format['fontweight'].lower() == 'bold' and self.fig_format['use_LaTex']:
            #matplotlib.rcParams.update({'font.family': self.fig_format['fontname'],'font.weight':self.fig_format['fontweight'],'font.size': self.fig_format['fontsize'], 'text.usetex': self.fig_format['use_LaTex'], 'text.latex.preamble':r"\usepackage{amsmath}"})
            matplotlib.rcParams.update({'font.family': self.fig_format['fontname'],'font.weight':self.fig_format['fontweight'],'font.size': self.fig_format['fontsize'], 'text.usetex': self.fig_format['use_LaTex'], 'text.latex.preamble':r'\boldmath'})
        else:
            matplotlib.rcParams.update({'font.family': self.fig_format['fontname'],'font.weight':self.fig_format['fontweight'],'font.size': self.fig_format['fontsize'], 'text.usetex': self.fig_format['use_LaTex']})

        # Custom font for mathtext
        if self.fig_format['mathtext_fontname'] != None:
            # The folloiwng lines changes the latext font to a custom font
            # Source: http://stackoverflow.com/questions/11367736/matplotlib-consistent-font-using-latex
            # http://matplotlib.org/users/mathtext.html#fonts
            matplotlib.rcParams.update({'mathtext.fontset':'custom', 'mathtext.rm':self.fig_format['mathtext_fontname'], 'mathtext.it':self.fig_format['mathtext_fontname'] + ':italic', 'mathtext.bf':self.fig_format['mathtext_fontname'] + ':bold'})


    def __setattr__(self,attr_name,attr_value):
        """
        Redefines funciton __setattr__
        INPUTS:
        -------
        attr_name: Attribute name
        attr_value: Attribute value
        """
        if attr_name.lower() == 'title' and not isinstance(attr_value,str):
            raise TypeError('title must be a string')
 
        if attr_name.lower() == 'title_format' and not isinstance(attr_value,dict):
            raise TypeError('title_format must be a dictionary')
        elif attr_name.lower() == 'title_format' and len([k for k in attr_value.keys() if k.lower() not in ['fontname','fontweight','fontsize','distance_from_graph']]) > 0: 
            raise ValueError('Unknown key(s) for label_format: {}'.format([k for k in attr_value.keys() if k.lower() not in ['fontname','fontweight','fontsize','distance_from_graph']]))
 
        if attr_name.lower() == 'xaxis' and attr_value != None and not isinstance(attr_value,axis):
            raise TypeError('xaxis must be either None or an instance of class axis')
        if attr_name.lower() == 'yaxis' and attr_value != None and not isinstance(attr_value,axis):
            raise TypeError('yaxis must be either None or an instance of class axis')
        if attr_name.lower() == 'zaxis' and attr_value != None and not isinstance(attr_value,axis):
            raise TypeError('zaxis must be either None or an instance of class axis')
 
        # Gridlines
        if attr_name.lower() == 'plot_gridlines' and not isinstance(attr_value,bool):
            raise TypeError('plot_gridlines must be either True or False')
 
        if attr_name.lower() == 'gridlines_format' and not isinstance(attr_value,dict):
            raise TypeError('gridlines_format must be a dictionary')
        elif attr_name.lower() == 'gridlines_format' and len([k for k in attr_value.keys() if k.lower() not in ['color', 'linestyle', 'linewidth']]) > 0: 
            raise ValueError('Unknown key(s) for label_format: {}'.format([k for k in attr_value.keys() if k.lower() not in ['color', 'linestyle', 'linewidth']]))

        # Legend
        if attr_name.lower() == 'show_legend' and not isinstance(attr_value,bool):
            raise TypeError('show_legend must be either True or False')

        if attr_name.lower() == 'legend_format' and attr_value != None and not isinstance(attr_value,dict):
           raise TypeError('legend_format must be a dicitionary')
        elif attr_name.lower() == 'legend_format' and attr_value != None and len([k for k in attr_value.keys() if k not in ['location','fontname','fontsize','fontweight','linewidth','bbox_to_anchor']]) > 0:
            raise ValueError('Invalid key(s) for legend_format: {}. Allowed keys are: [location, fontsize,fontweight,linewidth, bbox_to_anchor]'.format([k for k in attr_value.keys() if k not in ['location','fontname','fontsize','fontweight','linewidth','bbox_to_anchor']]))

        # Figure format
        if attr_name.lower() == 'fig_format' and not isinstance(attr_value,dict):
            raise TypeError('fig_format must be a dictionary')
        elif attr_name.lower() == 'fig_format' and len([k for k in attr_value.keys() if k.lower() not in ['figsize','fontname','fontsize','fontweight','use_tight_layout','dpi','use_latex','mathtext_fontname']]) > 0:
            raise ValueError('Unknown key(s) for fig_format: {}. Allowed keys include [figsize,use_tight_layout,dpi,use_latex]'.format([k for k in attr_value.keys() if k.lower() not in ['figsize','fontname','fontsize','fontweight','use_tight_layout','dpi','use_latex','mathtext_fontname']])) 
 
        # plot_type
        if attr_name.lower() == 'plot_type' and attr_value != None and not isinstance(attr_value,str):
            raise TypeError('plot_type must be a string')
        if attr_name.lower() == 'plot_type' and attr_value != None and attr_value.lower() not in ['heatmap', '2d_line', '3d']:
            raise ValueError('Invalid plot_type: {}. Allowed choices are: heatmap, 2d_line, 3d'.format(attr_value))
 
        if attr_name.lower() == 'output_filename' and not isinstance(attr_value,str):
            raise TypeError('output_filename must be a string')

        # Color bar
        if attr_name.lower() == 'clrbar' and not isinstance(attr_value, color_bar):
            raise TypeError('clorbar must be an instance of color_bar')

        self.__dict__[attr_name] = attr_value

    def plot2D(self, x, y, sort_data = False, label = '', line_format = {'color':None, 'width': 3, 'style': '-'}, show_markers = False, markers_format = {'style':'o','size':2,'facecolor':None,'edgewidth':None, 'edgecolor': None}, create_new_figure = True, save_current = True):
        """
        Creates a 2D line plot

        INPUTS:
        ------
                          x: x values (q list or a numpy array). NOTE: x values must be ordered
                          y: y values (a list or a numpy array)
                 sorte_data: If True, x is sorted and y is sorted according to the softed x
                      label: Label (legend text) for the current plotted line. The format of the label
                             is specified for the legend
                line_format: Line properties. This is a dicitonary with keys including 'style', 'color', 'width'
               show_markers: Indicates whether to show the markers
             markers_format: A dictionary containing the markers format information 
              spines_format: Format axes and graph's border lines
        save_current_figure: Indicates whether to make this plot in a new figure (True) or in a 
                             previously created one (False_
               save_current: If True runs customize_and_save at the end. This parameter can be set to 
                             False if this funciton is going to run for multiple datasets, where one 
                             wants to make several line plots in the same figure and save the final 
                             figure at the end. 
        """
        self._plot_type = '2d_line' 

        # x and y 
        if isinstance(x,list):
            x = np.array(x)
        elif not isinstance(x,np.ndarray):
            raise TypeError('x must be either a list or a numpy array')
        if isinstance(y,list):
            y = np.array(y)
        elif not isinstance(y,np.ndarray):
            raise TypeError('y must be either a list or a numpy array')
        # Make sure x andy are of the same size
        if x.shape != y.shape:
            raise userError('x and y are of the same size. size of x is: {}, size of of y is: {}'.format(x.shape,y.shape))
        # Sort x and accordingly y
        if sort_data:
            data_dict = dict([(x[i],y[i]) for i in range(x.shape[0])])
            X = np.array(sorted(data_dict.keys()))
            Y = np.array([data_dict[xx] for xx in X])
        else:
            X, Y = x, y 
        
        # label
        if not isinstance(label,str):
            raise TypeError('lable must be a string')

        # Whether to save the current figure
        if not isinstance(save_current,bool):
            raise TypeError('save_current must be either True or False')        
    
        # line format
        if not isinstance(line_format,dict):
           raise TypeError('line_format must be a dicitionary. A {} was provided'.format(type(line_format)))
        elif len([k for k in line_format.keys() if k not in ['color','width','style']]) > 0:
            raise ValueError('Invalid key(s) for line_format: {}. Allowed keys are: [color,width,style]'.format([k for k in line_format.keys() if k not in ['color','width','style']]))
    
        if 'width' not in line_format.keys():
            line_format['width'] = 2
        if 'color' not in line_format.keys():
            line_format['color'] = None
        if 'style' not in line_format.keys():
            line_format['style'] = '-'
    
        # markers
        if not isinstance(show_markers,bool):
            raise TypeError('show_markers must be either True or False')
        if not isinstance(markers_format,dict):
           raise TypeError('markers_format must be a dicitionary')
        elif len([k for k in markers_format.keys() if k not in ['style','size','facecolor','edgewidth','edgecolor']]) > 0:
            raise ValueError('Invalid key(s) for markers_format: {}. Allowed keys are: [style,size,facecolor,edgewidth,edgecolor]'.format([k for k in markers_format.keys() if k not in ['style','size','facecolor','edgewidth','edgecolor']]))
    
        if markers_format == None:
            markers_format = {'style': None,'size': None,'facecolor': None,'edgewidth': None,'edgecolor': None} 
        if 'style' not in markers_format.keys():
            markers_format['style'] = None
        if 'size' not in markers_format.keys():
            markers_format['size'] = None
        if 'facecolor' not in markers_format.keys():
            markers_format['facecolor'] = None
        if 'edgewidth' not in markers_format.keys():
            markers_format['edgewidth'] = None
        if 'edgecolor' not in markers_format.keys():
            markers_format['edgecolor'] = None

        # create_new_figure
        if not isinstance(create_new_figure, bool):
            raise TypeError('create_new_figure must be eithter True or False')

        # save_current_figure
        if not isinstance(save_current, bool):
            raise TypeError('save_current must be eithter True or False')

        if create_new_figure:
            self.fig, self.ax = plt.subplots(figsize = self.fig_format['figsize'])
        
        #--- Make the plot ----    
        if show_markers:
            self.ax.plot(X, Y, label = label, color = line_format['color'], linewidth = line_format['width'], linestyle = line_format['style'], marker = markers['style'], markersize = markers['size'], markerfacecolor = markers['facecolor'], markeredgewidth = markers['edgewidth'], markeredgecolor = markers['edgecolor'])  
        else:
            self.ax.plot(X, Y, label = label, color = line_format['color'], linewidth = line_format['width'], linestyle = line_format['style'])  

        # Save the current figure         
        if save_current:
            self.customize_and_save()


    def plot3D(self, x, y, z, sort_data = False, plot_func = 'plot_trisurf', line_format = {'color':None, 'width': 3, 'style': '-'}, clrbar = None, save_current = True):
        """
        Creates a 3D line plot

        INPUTS:
        ------
                     x: x values (q list or a numpy array). NOTE: x values must be ordered
                     y: y values (a list or a numpy array)
                     z: z values (a list or a numpy array)
            sorte_data: If True, x is sorted and y is sorted according to the softed x
             plot_func: 3D plot function to use. Allowed choices are plot, scatter, plot_wireframe,
                        plot_surface, plot_trisurf, contour, contourf,quiver  
           line_format: Line properties. This is a dicitonary with keys including 'style', 'color', 'width'
                clrbar: An instance of class colorbar 
        save_current: If True runs customize_and_save at the end. This parameter can be set to False if this funciton
                      is going to run for multiple datasets, where one wants to make several line plots in the same figure
                      and save the final figure at the end. 

        Sources:
        http://matplotlib.org/mpl_toolkits/mplot3d/tutorial.html
        http://stackoverflow.com/questions/21161884/plotting-a-3d-surface-from-a-list-of-tuples-in-matplotlib
        """
        from mpl_toolkits.mplot3d import Axes3D

        self._plot_type = '3d' 

        # x, y and z 
        if isinstance(x,list):
            x = np.array(x)
        elif not isinstance(x,np.ndarray):
            raise TypeError('x must be either a list or a numpy array')
        if isinstance(y,list):
            y = np.array(y)
        elif not isinstance(y,np.ndarray):
            raise TypeError('y must be either a list or a numpy array')
        if isinstance(z,list):
            z = np.array(z)
        elif not isinstance(z,np.ndarray):
            raise TypeError('z must be either a list or a numpy array')
        # Make sure x andy are of the same size
        if x.shape != y.shape or x.shape != z.shape or y.shape != z.shape:
            raise userError('x, y and z are not of the same size. size of x is: {}, size of of y is: {}, size of z is: {}'.format(x.shape,y.shape,z.shape))

        # Sort x and accordingly y
        if sort_data:
            data_dict = dict([((x[i],y[i]), z[i]) for i in range(x.shape[0])])
            # Here, we first sort x and then find the y values corresponding to the sorted x
            # Next, we sort y and find the z values corresponidng ot sorted y
            sorted_x_y = sorted(data_dict.keys() , key = lambda m:(m[0],m[1]))
            X = np.array([s[0] for s in sorted_x_y]) 
            Y = np.array([s[1] for s in sorted_x_y]) 
            Z = np.array([data_dict[xy] for xy in sorted_x_y])
        else:
            X, Y, Z = x, y, z 

        # plot_func
        if not isinstance(plot_func,str):
            raise TypeError('plot_func must be a string')
        elif plot_func not in ['plot', 'scatter', 'plot_wireframe', 'plot_surface', 'plot_trisurf', 'contour', 'contourf','quiver']:
            raise ValueError('Invalid plot_func value! Allowed choices are [plot, scatter, plot_wireframe, plot_surface, plot_trisurf, contour, contourf,quiver]')

        # line format
        if not isinstance(line_format,dict):
           raise TypeError('line_format must be a dicitionary. A {} was provided'.format(type(line_format)))
        elif len([k for k in line_format.keys() if k not in ['color','width','style']]) > 0:
            raise ValueError('Invalid key(s) for line_format: {}. Allowed keys are: [color,width,style]'.format([k for k in line_format.keys() if k not in ['color','width','style']]))
    
        if 'width' not in line_format.keys():
            line_format['width'] = 2
        if 'color' not in line_format.keys():
            line_format['color'] = None
        if 'style' not in line_format.keys():
            line_format['style'] = '-'
    
        # Color bar
        self.clrbar = clrbar
        if self.clrbar.colorlimits == None:
            clrbar.colorlimits = (z.min(),z.max())

        self.fig = plt.figure(figsize = self.fig_format['figsize'])
        self.ax = self.fig.gca(projection='3d')

        #------ plot_trisurf -------
        self.ax.plot_trisurf(X, Y, Z, linewidth = line_format['width'], linestyle = line_format['style'], cmap = clrbar.colormap, vmin = clrbar.colorlimits[0], vmax = clrbar.colorlimits[1])

        plt.draw()

    def heatmap(self, x, y, data, plot_func = 'pcolor', interpolate = False, clrbar = None):
        """
        Creates a heatmap of data 
     
        INPUTS:
        ------
                  x: horizontal axis data (list or numpy array)
                  y: Vertical axis data (list or numpy array)
               data: 2-D numpy array where ROWS correspond to data on VERTIACAL axis and COLUMNS correspond to data on
                     HORIZONTAL axis.
          plot_func: Eligible choices are 'pcolor', 'pcolormesh', 'matshow', 'imshow'
             clrbar: An instance of class colorbar 
        interpolate: Interpolate data (applicable only to matshow and imshow)
        """
        self.fig, self.ax = plt.subplots(figsize = self.fig_format['figsize'])
    
        self._plot_type = 'heatmap'
    
        if clrbar == None:
            clrbar = color_bar()
        self.clrbar = clrbar

        # Color bar
        if self.clrbar.colorlimits == None:
            self.clrbar.colorlimits = (data.min(),data.max())
 
        # If you like to have ticklabels placed in the middle of ticks (e.g., in a discerete colormap for pcolor 
        # and pcolormesh), use: vmin = data.min() - 0.5 and vmax = data.max() + 0.5. 
        if self.clrbar.ticklabels_format['ticks_position'] == 'middle':
            (vmin,vmax) = (self.clrbar.colorlimits[0] - 0.5, self.clrbar.colorlimits[1] + 0.5)
        else:
            (vmin,vmax) = (self.clrbar.colorlimits[0], self.clrbar.colorlimits[1])

        #------------ pcolor and pcolormesh ----------------
        if plot_func.lower() in ['pcolor','pcolormesh']:
            if plot_func.lower() == 'pcolor':
                pc = self.ax.pcolor(x,y,data, cmap = clrbar.colormap, vmin = vmin, vmax = vmax)
            elif plot_func.lower() == 'pcolormesh':
                pc = self.ax.pcolormesh(x,y,data, cmap = clrbar.colormap, vmin = vmin, vmax = vmax) 
    
            #-- colorbar --
            if clrbar.custom_ticklabels != None:
                # Source: http://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
                self.cbar = self.fig.colorbar(pc,ticks = np.arange(np.min(data),np.max(data)+1))
            else:
                self.cbar = self.fig.colorbar(pc, ticks = self.clrbar.custom_ticks)

        #------------ imshow and matshow ----------------
        elif plot_func.lower() in ['matshow','imshow']:
    
            if plot_func.lower() == 'matshow':
                if interpolate:
                    ms = self.ax.matshow(data, cmap = clrbar.colormap, vmin = vmin, vmax = vmax, interpolation = 'bilinear')
                else: 
                    ms = self.ax.matshow(data, cmap = clrbar.colormap, vmin = vmin, vmax = vmax)
    
            elif plot_func.lower() == 'imshow':
                if interpolate:
                    ms = self.ax.imshow(data, cmap = clrbar.colormap, vmin = clrbar.colorlimits[0], vmax = clrbar.colorlimits[1], interpolation = 'bilinear')
                else:
                    ms = self.ax.imshow(data, cmap = clrbar.colormap, vmin = clrbar.colorlimits[0], vmax = clrbar.colorlimits[1])

            #-- colorbar --
            if self.clrbar.custom_ticklabels != None:
                # Source: http://stackoverflow.com/questions/14777066/matplotlib-discrete-colorbar
                self.cbar = self.fig.colorbar(ms,ticks = np.arange(np.min(data),np.max(data)+1))
            else:
                self.cbar = self.fig.colorbar(ms, ticks = self.clrbar.custom_ticks)

            #--- Put the major ticks at the middle of each cell (don't do this for pcolor and pcolormesh) ---
            # Source: http://stackoverflow.com/questions/14391959/heatmap-in-matplotlib-with-pcolor
            self.ax.set_xticks(np.arange(data.shape[1]), minor = False)
            self.ax.set_yticks(np.arange(data.shape[0]), minor = False)
    
            # Custom tick labels
            if self.xaxis.custom_ticklabels == None:
                self.xaxis.custom_ticklabels = [str(i) for i in x]
            if self.yaxis.custom_ticklabels == None:
                self.yaxis.custom_ticklabels = [str(i) for i in y]
    
        # Customize other properties of the plot
        self.customize_and_save()
                
    def customize_and_save(self):
        """
        This function customizes axes
        """
        # Make the plot fill the picture
        if self.fig_format['use_tight_layout']:
            self.ax.set_aspect('auto')
            self.fig.tight_layout() 

        #--- Axes limits ---
        if self.xaxis.limits != None: 
            self.ax.set_xlim([self.xaxis.limits[0],self.xaxis.limits[1]])
        if self.yaxis.limits != None: 
            self.ax.set_ylim([self.yaxis.limits[0],self.yaxis.limits[1]])
        if self._plot_type.lower() == '3d' and self.zaxis.limits != None:
            self.ax.set_zlim([self.zaxis.limits[0],self.zaxis.limits[1]])

        #---- Ticks and tick labels ---- 
        #-- Set the position of tick labels (top, bottom, left, right) --
        # x axis
        if self.xaxis.ticklabels_format['position'] != None and self.xaxis.ticklabels_format['position'].lower() == 'bottom':
            self.ax.xaxis.tick_bottom()
        elif self.xaxis.ticklabels_format['position'] != None and self.xaxis.ticklabels_format['position'].lower() == 'top':
            self.ax.xaxis.tick_top()
        elif self.xaxis.ticklabels_format['position'] != None and self.xaxis.ticklabels_format['position'].lower() not in ['bottom','top']:
            raise userError('Invalid ticklabels position for the x axis: {}. Allowed choices are top and bottom'.format(self.xaxis.ticklabels_format['position']))

        # y axis
        if self.yaxis.ticklabels_format['position'] != None and self.yaxis.ticklabels_format['position'].lower() == 'left':
            self.ax.yaxis.tick_left()
        elif self.yaxis.ticklabels_format['position'] != None and self.yaxis.ticklabels_format['position'].lower() == 'right':
            self.ax.yaxis.tick_right()
        elif self.yaxis.ticklabels_format['position'] != None and self.yaxis.ticklabels_format['position'].lower() not in ['left','right']:
            raise userError('Invalid ticklabels position for the y axis: {}. Allowed choices are left and right'.format(self.yaxis.ticklabels_format['position']))

        # z axis
        if self._plot_type.lower() == '3d':
            if self.zaxis.ticklabels_format['position'] != None and self.zaxis.ticklabels_format['position'].lower() == 'left':
                self.ax.zaxis.tick_left()
            elif self.zaxis.ticklabels_format['position'] != None and self.zaxis.ticklabels_format['position'].lower() == 'right':
                self.ax.zaxis.tick_right()
            elif self.zaxis.ticklabels_format['position'] != None and self.zaxis.ticklabels_format['position'].lower() not in ['left','right']:
                raise userError('Invalid ticklabels position for the y axis: {}. Allowed choices are left and right'.format(self.zaxis.ticklabels_format['position']))

        #-- Set the alignment of tick labels ---
        # This can be used to make the ticklabels closer to the axes (useful for 3d plots)
        # Source: http://stackoverflow.com/questions/16061349/tick-label-positions-for-matplotlib-3d-plot  
        # x axis
        if self.xaxis.ticklabels_format['horizontalalignment'] != None:
            for ticklabel in self.ax.get_xmajorticklabels():
                ticklabel.set_horizontalalignment(self.xaxis.ticklabels_format['horizontalalignment'])     
        if self.xaxis.ticklabels_format['verticalalignment'] != None:
            for ticklabel in self.ax.get_ymajorticklabels():
                ticklabel.set_verticalalignment(self.xaxis.ticklabels_format['verticalalignment'])     

        # y axis
        if self.yaxis.ticklabels_format['horizontalalignment'] != None:
            for ticklabel in self.ax.get_ymajorticklabels():
                ticklabel.set_horizontalalignment(self.yaxis.ticklabels_format['horizontalalignment'])     
        if self.yaxis.ticklabels_format['verticalalignment'] != None:
            for ticklabel in self.ax.get_ymajorticklabels():
                ticklabel.set_verticalalignment(self.yaxis.ticklabels_format['verticalalignment'])     

        # z axis
        if self._plot_type.lower() == '3d':
            if self.zaxis.ticklabels_format['horizontalalignment'] != None:
                for ticklabel in self.ax.get_zmajorticklabels():
                    ticklabel.set_horizontalalignment(self.zaxis.ticklabels_format['horizontalalignment'])     
            if self.zaxis.ticklabels_format['verticalalignment'] != None:
                for ticklabel in self.ax.get_zmajorticklabels():
                    ticklabel.set_verticalalignment(self.zaxis.ticklabels_format['verticalalignment'])     

        #-- Set major and minor ticks and tick labels --
        # Source: http://matplotlib.org/examples/pylab_examples/major_minor_demo1.html
        #- x ticks --
        # Minor tick labels
        if self.xaxis.set_minorticks:
            if self.xaxis.minorticks_spacing != None:
                x_minorLocator = MultipleLocator(self.xaxis.minorticks_spacing)
                self.ax.xaxis.set_minor_locator(x_minorLocator)
                if self.xaxis.ticklabels_format['string_format'] != None: 
                    x_minorFormatter = FormatStrFormatter(self.xaxis.ticklabels_format['string_format'])
                    self.ax.xaxis.set_minor_formatter(x_minorFormatter)
            for ticklabel in self.ax.get_xminorticklabels():
                ticklabel.set_fontsize(self.xaxis.ticklabels_format['fontsize']) 
                ticklabel.set_fontweight(self.xaxis.ticklabels_format['fontweight']) 
                ticklabel.set_rotation(self.xaxis.ticklabels_format['rotation']) 
        # Major tick labels
        if self.xaxis.majorticks_spacing != None:
            x_majorLocator = MultipleLocator(self.xaxis.majorticks_spacing)
            self.ax.xaxis.set_major_locator(x_majorLocator)
            if self.xaxis.ticklabels_format['string_format'] != None: 
                x_majorFormatter = FormatStrFormatter(self.xaxis.ticklabels_format['string_format'])
                self.ax.xaxis.set_major_formatter(x_majorFormatter)
        for ticklabel in self.ax.get_xmajorticklabels():
            ticklabel.set_fontsize(self.xaxis.ticklabels_format['fontsize']) 
            ticklabel.set_fontweight(self.xaxis.ticklabels_format['fontweight']) 
            ticklabel.set_rotation(self.xaxis.ticklabels_format['rotation']) 
    
        #- y ticks --
        # Minor tick labels
        if self.yaxis.set_minorticks:
            if self.yaxis.minorticks_spacing != None:
                y_minorLocator = MultipleLocator(self.yaxis.minorticks_spacing)
                self.ax.yaxis.set_minor_locator(y_minorLocator)
                if self.yaxis.ticklabels_format['string_format'] != None: 
                    y_minorFormatter = FormatStrFormatter(self.yaxis.ticklabels_format['string_format'])
                    self.ax.yaxis.set_minor_formatter(y_minorFormatter)
            for ticklabel in self.ax.get_yminorticklabels():
                ticklabel.set_fontsize(self.yaxis.ticklabels_format['fontsize']) 
                ticklabel.set_fontweight(self.yaxis.ticklabels_format['fontweight']) 
                ticklabel.set_rotation(self.yaxis.ticklabels_format['rotation']) 
        # Major tick labels
        if self.yaxis.majorticks_spacing != None:
            y_majorLocator = MultipleLocator(self.yaxis.majorticks_spacing)
            self.ax.yaxis.set_major_locator(y_majorLocator)
            if self.yaxis.ticklabels_format['string_format'] != None: 
                y_majorFormatter = FormatStrFormatter(self.yaxis.ticklabels_format['string_format'])
                self.ax.yaxis.set_major_formatter(y_majorFormatter)
        for ticklabel in self.ax.get_ymajorticklabels():
            ticklabel.set_fontsize(self.yaxis.ticklabels_format['fontsize']) 
            ticklabel.set_fontweight(self.yaxis.ticklabels_format['fontweight']) 
            ticklabel.set_rotation(self.yaxis.ticklabels_format['rotation']) 

        #- z ticks --
        if self._plot_type.lower() == '3d':
            # Minor tick labels
            if self.zaxis.set_minorticks:
                if self.zaxis.minorticks_spacing != None:
                    z_minorLocator = MultipleLocator(self.zaxis.minorticks_spacing)
                    self.ax.zaxis.set_minor_locator(z_minorLocator)
                    if self.zaxis.ticklabels_format['string_format'] != None: 
                        z_minorFormatter = FormatStrFormatter(self.zaxis.ticklabels_format['string_format'])
                        self.ax.zaxis.set_minor_formatter(z_minorFormatter)
                for ticklabel in self.ax.get_zminorticklabels():
                    ticklabel.set_fontsize(self.zaxis.ticklabels_format['fontsize']) 
                    ticklabel.set_fontweight(self.zaxis.ticklabels_format['fontweight']) 
                    ticklabel.set_rotation(self.zaxis.ticklabels_format['rotation']) 
            # Major tick labels
            if self.zaxis.majorticks_spacing != None:
                z_majorLocator = MultipleLocator(self.zaxis.majorticks_spacing)
                self.ax.zaxis.set_major_locator(z_majorLocator)
                if self.zaxis.ticklabels_format['string_format'] != None: 
                    z_majorFormatter = FormatStrFormatter(self.zaxis.ticklabels_format['string_format'])
                    self.ax.zaxis.set_major_formatter(z_majorFormatter)
            for ticklabel in self.ax.get_zmajorticklabels():
                ticklabel.set_fontsize(self.zaxis.ticklabels_format['fontsize']) 
                ticklabel.set_fontweight(self.zaxis.ticklabels_format['fontweight']) 
                ticklabel.set_rotation(self.yaxis.ticklabels_format['rotation']) 

        #--- Cutom ticks ---
        if self.xaxis.custom_ticks != None:
            self.ax.set_xticks(self.xaxis.custom_ticks)
        if self.yaxis.custom_ticks != None:
            self.ax.set_yticks(self.yaxis.custom_ticks)
        if self._plot_type.lower() == '3d' and  self.zaxis.custom_ticks != None:
            self.ax.set_zticks(self.zaxis.custom_ticks)
    
        #-- Custom tick labels --
        # x axis
        if self.xaxis.custom_ticklabels != None:
            self.ax.set_xticklabels(self.xaxis.custom_ticklabels,fontsize = self.xaxis.ticklabels_format['fontsize'], weight = self.xaxis.ticklabels_format['fontweight'])

        # y axis
        if self.yaxis.custom_ticklabels != None:
            self.ax.set_yticklabels(self.yaxis.custom_ticklabels, fontsize = self.yaxis.ticklabels_format['fontsize'], weight = self.yaxis.ticklabels_format['fontweight'])

        # z axis
        if self._plot_type.lower() == '3d' and self.zaxis.custom_ticklabels != None:
            self.ax.set_zticklabels(self.zaxis.custom_ticklabels, fontsize = self.zaxis.ticklabels_format['fontsize'], weight = self.zaxis.ticklabels_format['fontweight'])

            self.ax.set_yticklabels([str(i) for i in y], fontsize = self.yaxis.ticklabels_format['fontsize'], weight = self.yaxis.ticklabels_format['fontweight'])

        #---- Title and Axes labels ----
        # Title
        # distance_from_graph represents the distance between the title and the graph
        # Source: http://stackoverflow.com/questions/12750355/python-matplotlib-figure-title-overlaps-axes-label-when-using-twiny
        if self.title != '':    
            self.ax.set_title(self.title,{'weight':self.title_format['fontweight'],'size':self.title_format['fontsize']}, y = self.title_format['distance_from_graph'])

        # x axis label
        if self.xaxis.label != '':    
            if self.xaxis.label_format['distance_from_ticklabels'] != None: 
                self.ax.set_xlabel(self.xaxis.label,{'fontname': self.xaxis.label_format['fontname'],'weight':self.xaxis.label_format['fontweight'],'size':self.xaxis.label_format['fontsize']}, labelpad = self.xaxis.label_format['distance_from_ticklabels'])
            else:
                self.ax.set_xlabel(self.xaxis.label,{'fontname': self.xaxis.label_format['fontname'],'weight':self.xaxis.label_format['fontweight'],'size':self.xaxis.label_format['fontsize']})

        # y axis label
        if self.yaxis.label != '':    
            if self.yaxis.label_format['distance_from_ticklabels'] != None: 
                self.ax.set_ylabel(self.yaxis.label,{'fontname': self.yaxis.label_format['fontname'], 'weight':self.yaxis.label_format['fontweight'],'size':self.yaxis.label_format['fontsize']}, labelpad = self.yaxis.label_format['distance_from_ticklabels'])
            else:
                self.ax.set_ylabel(self.yaxis.label,{'fontname': self.yaxis.label_format['fontname'],'weight':self.yaxis.label_format['fontweight'],'size':self.yaxis.label_format['fontsize']})

        # z axis label
        if self._plot_type.lower() == '3d' and self.zaxis.label != '':    
            if self.zaxis.label_format['distance_from_ticklabels'] != None: 
                self.ax.set_zlabel(self.zaxis.label,{'fontname': self.zaxis.label_format['fontname'], 'weight':self.zaxis.label_format['fontweight'],'size':self.zaxis.label_format['fontsize']}, labelpad = self.zaxis.label_format['distance_from_ticklabels'])
            else:
                self.ax.set_zlabel(self.zaxis.label,{'fontname': self.zaxis.label_format['fontname'],'weight':self.zaxis.label_format['fontweight'],'size':self.zaxis.label_format['fontsize']})

        #-- Position of axes labels --
        # x axis
        if self.xaxis.label_format['position'] != None and self.xaxis.label_format['position'].lower() == 'top':
            self.ax.xaxis.set_label_position('top') 
        elif self.xaxis.label_format['position'] != None and self.xaxis.label_format['position'].lower() == 'bottom':
            self.ax.xaxis.set_label_position('bottom') 
        elif self.xaxis.label_format['position'] != None and self.xaxis.label_format['position'].lower() not in ['top','bottom']:
            raise ValueError('Invalid label position for the x axis: {}. Allowed choices are top and bottom'.format(self.xaxis.label_format['position']))

        # y axis
        if self.yaxis.label_format['position'] != None and self.yaxis.label_format['position'].lower() == 'right':
            self.ax.yaxis.set_label_position('right') 
        elif self.yaxis.label_format['position'] != None and self.yaxis.label_format['position'].lower() == 'left':
            self.ax.yaxis.set_label_position('left') 
        elif self.yaxis.label_format['position'] != None and self.yaxis.label_format['position'].lower() not in ['left','right']:
            raise ValueError('Invalid label position for the y axis: {}. Allowed choices are right and left'.format(self.yaxis.label_format['position']))

        # z axis
        if self._plot_type.lower() == '3d':
            if self.zaxis.label_format['position'] != None and self.zaxis.label_format['position'].lower() == 'right':
                self.ax.zaxis.set_label_position('right') 
            elif self.zaxis.label_format['position'] != None and self.zaxis.label_format['position'].lower() == 'left':
                self.ax.zaxis.set_label_position('left') 
            elif self.zaxis.label_format['position'] != None and self.zaxis.label_format['position'].lower() not in ['left','right']:
                raise ValueError('Invalid label position for the z axis: {}. Allowed choices are right and left'.format(self.zaxis.label_format['position']))

        #---- Invert axis ----
        if self.xaxis.invert:
            self.ax.invert_xaxis()
        if self.yaxis.invert:
            self.ax.invert_yaxis()
        if self._plot_type.lower() == '3d' and self.zaxis.invert:
            self.ax.invert_zaxis()

        #--- Axes scales ----
        if self.xaxis.scale != None:
            self.ax.set_xscale(self.xaxis.scale)
        if self.yaxis.scale != None:
            self.ax.set_yscale(self.yaxis.scale)
        if self._plot_type.lower() == '3d' and self.zaxis.scale != None:
            self.ax.set_zscale(self.zaxis.scale)

        #---- Grid lines ----
        if self._plot_type.lower() != '3d':
            # Gridlines format
            # From: http://stackoverflow.com/questions/9127434/how-to-create-major-and-minor-gridlines-with-different-linestyles-in-python
            # From: http://stackoverflow.com/questions/17925545/adjusting-gridlines-on-a-3d-matplotlib-figure
            # gridlines_format = {'color':'k', 'linestyle':'dashed', 'linewidth':2}   
            if self.xaxis.plot_gridlines:
                self.ax.xaxis.grid(True, color = self.xaxis.gridlines_format['color'], linestyle = self.xaxis.gridlines_format['linestyle'], linewidth = self.xaxis.gridlines_format['linewidth'])
            if self.yaxis.plot_gridlines:
                self.ax.yaxis.grid(True, color = self.yaxis.gridlines_format['color'], linestyle = self.yaxis.gridlines_format['linestyle'], linewidth = self.yaxis.gridlines_format['linewidth'])

        else: 
            # for 1 3d plot
            # From: http://stackoverflow.com/questions/17925545/adjusting-gridlines-on-a-3d-matplotlib-figure
            if self.xaxis.plot_gridlines:
                self.ax.w_xaxis.gridlines.set_lw(self.xaxis.gridlines_format['linewidth'])
                self.ax.w_xaxis.gridlines.set_linestyle(self.xaxis.gridlines_format['linestyle'])
                self.ax.w_xaxis._axinfo.update({'grid': {'color': self.xaxis.gridlines_format['color']}})

            if self.yaxis.plot_gridlines:
                self.ax.w_yaxis.gridlines.set_linestyle(self.yaxis.gridlines_format['linestyle'])
                self.ax.w_yaxis.gridlines.set_lw(self.yaxis.gridlines_format['linewidth'])
                self.ax.w_yaxis._axinfo.update({'grid': {'color': self.yaxis.gridlines_format['color']}})

            if self.zaxis.plot_gridlines:
                self.ax.w_zaxis.gridlines.set_lw(self.zaxis.gridlines_format['linewidth'])
                self.ax.w_zaxis.gridlines.set_linestyle(self.zaxis.gridlines_format['linestyle'])
                self.ax.w_zaxis._axinfo.update({'grid': {'color': self.zaxis.gridlines_format['color']}})

        #--- Legend ---
        # Source: http://matplotlib.org/api/legend_api.html
        # and     http://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
        # When using bbox_to_anchor and loc toghether, beware that if you like your legend to appear on the right
        # you should specify "left" in loc and vice versa. For example, to have your legend appear in the center right
        # you should specify "cetner left" for loc.
        if self.show_legend:
            if self.legend_format['bbox_to_anchor'] != None:
                lgd = self.ax.legend(loc = self.legend_format['location'], prop = {'family': self.legend_format['fontname'], 'size': self.legend_format['fontsize'],'weight':self.legend_format['fontweight']},handlelength = self.legend_format['linewidth'], bbox_to_anchor = self.legend_format['bbox_to_anchor'])
            else:
                lgd = self.ax.legend(loc = self.legend_format['location'], prop = {'size': self.legend_format['fontsize'],'weight':self.legend_format['fontweight']},handlelength = self.legend_format['linewidth'])

        #--- Set axes and borders (spines) properties ---
        # Set zorder to a large value so it is drawn after gridnlines (otherwise if you change the color of spiens, gridlines
        # will show on top of the spines)
        # Source: http://matplotlib.org/api/artist_api.html#matplotlib.artist.Artist.set_zorder
        # Bottom
        self.ax.spines['bottom'].set_zorder(10)
        if 'bottom' in self.xaxis.spines_format.keys():
            if 'linewidth' in  self.xaxis.spines_format['bottom'].keys():
                self.ax.spines['bottom'].set_linewidth(self.xaxis.spines_format['bottom']['linewidth'])
            if 'linestyle' in  self.xaxis.spines_format['bottom'].keys():
                self.ax.spines['bottom'].set_linestyle(self.xaxis.spines_format['bottom']['linestyle'])
            if 'linecolor' in  self.xaxis.spines_format['bottom'].keys():
                self.ax.spines['bottom'].set_color(self.xaxis.spines_format['bottom']['linecolor'])

        # Top
        self.ax.spines['top'].set_zorder(10)
        if 'top' in self.xaxis.spines_format.keys():
            if 'linewidth' in  self.xaxis.spines_format['top'].keys():
                self.ax.spines['top'].set_linewidth(self.xaxis.spines_format['top']['linewidth'])
            if 'linestyle' in  self.xaxis.spines_format['top'].keys():
                self.ax.spines['top'].set_linestyle(self.xaxis.spines_format['top']['linestyle'])
            if 'linecolor' in  self.xaxis.spines_format['top'].keys():
                self.ax.spines['top'].set_color(self.xaxis.spines_format['top']['linecolor'])

        # Left
        self.ax.spines['left'].set_zorder(10)
        if 'left' in self.yaxis.spines_format.keys():
            if 'linewidth' in  self.yaxis.spines_format['left'].keys():
                self.ax.spines['left'].set_linewidth(self.yaxis.spines_format['left']['linewidth'])
            if 'linestyle' in  self.yaxis.spines_format['left'].keys():
                self.ax.spines['left'].set_linestyle(self.yaxis.spines_format['left']['linestyle'])
            if 'linecolor' in  self.yaxis.spines_format['left'].keys():
                self.ax.spines['left'].set_color(self.yaxis.spines_format['left']['linecolor'])

        # Right
        self.ax.spines['right'].set_zorder(10)
        if 'left' in self.yaxis.spines_format.keys():
            if 'linewidth' in  self.yaxis.spines_format['right'].keys():
                self.ax.spines['right'].set_linewidth(self.yaxis.spines_format['right']['linewidth'])
            if 'linestyle' in  self.yaxis.spines_format['right'].keys():
                self.ax.spines['right'].set_linestyle(self.yaxis.spines_format['right']['linestyle'])
            if 'linecolor' in  self.yaxis.spines_format['right'].keys():
                self.ax.spines['right'].set_color(self.yaxis.spines_format['right']['linecolor'])


        #--- Customize colorbar ---
        if hasattr(self,'cbar'): 
            # Coloarbar label format
            if self.clrbar.label_format['distance_from_ticklabels'] != None: 
                self.cbar.set_label(self.clrbar.label,size = self.clrbar.label_format['fontsize'], weight = self.clrbar.label_format['fontweight'], rotation = self.clrbar.label_format['rotation'], labelpad = self.clrbar.label_format['distance_from_ticklabels'])
            else:
                self.cbar.set_label(self.clrbar.label, size = self.clrbar.label_format['fontsize'], weight = self.clrbar.label_format['fontweight'], rotation = self.clrbar.label_format['rotation'])

            # Set the position of tick labels (top, bottom, left, right)
            if self.clrbar.ticklabels_format['axis_position'] != None and self.clrbar.ticklabels_format['axis_position'].lower() == 'left':
                self.cbar.ax.yaxis.tick_left()
            elif self.clrbar.ticklabels_format['axis_position'] != None and self.clrbar.ticklabels_format['axis_position'].lower() == 'right':
                self.cbar.ax.yaxis.tick_right()
            elif self.clrbar.ticklabels_format['axis_position'] != None and self.clrbar.ticklabels_format['axis_position'].lower() not in ['left','right']:
                raise userError('Invalid ticklabels position for the y axis: {}. Allowed choices are left and right'.format(self.clrbar.ticklabels_format['position']))

            #- Set major and minor ticks and tick labels -
            # ***The following disrupts the color limit for some reasons. Use custom ticks instead
            # Source: http://matplotlib.org/examples/pylab_examples/major_minor_demo1.html
            # Minor tick labels
            if self.clrbar.set_minorticks:
                if self.clrbar.minorticks_spacing != None:
                    y_minorLocator = MultipleLocator(self.clrbar.minorticks_spacing)
                    self.cbar.ax.yaxis.set_minor_locator(y_minorLocator)
                    if self.clrbar.ticklabels_format['string_format'] != None: 
                        y_minorFormatter = FormatStrFormatter(self.clrbar.ticklabels_format['string_format'])
                        self.cbar.ax.yaxis.set_minor_formatter(y_minorFormatter)
                for ticklabel in self.cbar.ax.get_yminorticklabels():
                    ticklabel.set_fontsize(self.clrbar.ticklabels_format['fontsize']) 
                    ticklabel.set_fontweight(self.clrbar.ticklabels_format['fontweight']) 
                    ticklabel.set_rotation(self.clrbar.ticklabels_format['rotation']) 
            # Major tick labels
            if self.clrbar.majorticks_spacing != None:
                y_majorLocator = MultipleLocator(self.clrbar.majorticks_spacing)
                self.cbar.ax.yaxis.set_major_locator(y_majorLocator)
                if self.clrbar.ticklabels_format['string_format'] != None: 
                    y_majorFormatter = FormatStrFormatter(self.clrbar.ticklabels_format['string_format'])
                    self.cbar.ax.yaxis.set_major_formatter(y_majorFormatter)
            for ticklabel in self.cbar.ax.get_ymajorticklabels():
                ticklabel.set_fontsize(self.clrbar.ticklabels_format['fontsize']) 
                ticklabel.set_fontweight(self.clrbar.ticklabels_format['fontweight']) 
                ticklabel.set_rotation(self.clrbar.ticklabels_format['rotation']) 

            # Custom ticks --> ** This doesn't work for some reasons. Use ticks when defining the colorbar instead 
            if self.clrbar.custom_ticks != None:
                self.cbar.ax.set_yticks(self.clrbar.custom_ticks)

            # Custom tick labels 
            if self.clrbar.custom_ticklabels != None:
                self.cbar.ax.set_yticklabels(self.clrbar.custom_ticklabels, fontsize = self.clrbar.ticklabels_format['fontsize'], weight = self.clrbar.ticklabels_format['fontweight'])

        #--- Save the figure ---
        if self.output_filename != '':
            # If bbox_inches = 'tight' is not included part of the axes labels may be cut off in the saved image
            # Similarly, using bbox_extra_artists is to save part of the legend that is cut off
            # Source: http://stackoverflow.com/questions/4700614/how-to-put-the-legend-out-of-the-plot
            if self.fig_format['dpi'] == None:
                if self.show_legend:
                    self.fig.savefig(self.output_filename, bbox_extra_artists=(lgd,), bbox_inches='tight')
                else:
                    self.fig.savefig(self.output_filename, bbox_inches='tight')
            else:
                if self.show_legend:
                    self.fig.savefig(self.output_filename, dpi = self.fig_format['dpi'], bbox_extra_artists=(lgd,), bbox_inches='tight')
                else:
                    self.fig.savefig(self.output_filename, dpi = self.fig_format['dpi'], bbox_inches='tight')

        # Use draw instead of show so the graph gets updated it this is run after customzie_and_save
        plt.draw()

 
    def add_text(self,text = '', text_format = {'x_pos':None, 'y_pos':None, 'fontsize':figure_default_fontsize, 'fontweight':figure_default_fontweight}):
        """
        Adds a text to the graph

        INPUTS:
        -------
               text: Text to be added to the graph
        text_format: A dictionary containing the text format information 
        """
        if not isinstance(text,str):
            raise TypeError('text must be a string')

        if not isinstance(text_format,dict):
           raise TypeError('text_format must be a dicitionary')
        elif len([k for k in text_format.keys() if k not in ['x_pos','y_pos','fontsize','fontweight']]) > 0:
            raise ValueError('Invalid key(s) for text_format: {}. Allowed keys are: [x_pos, y_pos, fontsize, fontweight]'.format([k for k in text_format.keys() if k not in ['x_pos','y_pos','fontsize','fontweight']]))

        if 'x_pos' not in text_format.keys():
            text_format['x_pos'] = None
        if 'y_pos' not in text_format.keys():
            text_format['y_pos'] = None
        if text_format != None and 'fontsize' not in text_format.keys():
            text_format['fontsize'] = figure_default_fontsize
        if text_format != None and 'fontweight' not in text_format.keys():
            text_format['fontweight'] = figure_default_fontweight

        self.ax.text(text_format['x_pos'], text_format['y_pos'], text, fontsize = text_format['fontsize'], fontweight = text_format['fontweight'])

        # Use draw instead of show so the graph gets updated it this is run after customzie_and_save
        plt.draw()

    def shade(self,x_range = (), y_range = (),  shade_format = {'color':'green', 'transparency':0.5}):
        """
        Shades the area between x1 and x2 and y1 and y2

        INPUTS:
        -------
             x_range: A tuple with two elements showing the x range between which should be shaded.
                      The second element must be strictly greater than the first one.
             y_range: A tuple with two elements showing the y range between which should be shaded.
                      The second element must be strictly greater than the first one.
        shade_format: The format of shaded area. Here, transparency must be an integer or float
                      between 0 (fully transparent) to 1 (fully opaque).i.
                      See: http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.axvspan

        See: http://stackoverflow.com/questions/23248435/fill-between-two-vertical-lines-in-matplotlib
        """
        # x_range 
        if not isinstance(x_range,tuple): 
            raise TypeError('x_range must be a tuple')
        elif len(x_range) == 2 and x_range[1] <= x_range[0]: 
            raise TypeError('The second entry in x_range must be strictly greater than the first one')
       
        # y_range 
        if not isinstance(y_range,tuple): 
            raise TypeError('y_range must be a tuple')
        elif len(y_range) == 2 and y_range[1] <= y_range[0]: 
            raise TypeError('The second entry in y_range must be strictly greater than the first one')

        # shade_format
        if not isinstance(shade_format,dict):
            raise TypeError('shade_format must be a dictionary')
        elif len([k for k in shade_format.keys() if k.lower() not in ['color','transparency']]) > 0:
            raise ValueError('Invalid key for shade_fomrat: {}. Allowed keys are [color, transparency]', [k for k in shade_format.keys() if k.lower() not in ['color','transparency']])
        elif shade_format['transparency'] < 0 or shade_format['transparency'] > 1:
           raise ValueError("shade_format['transparency'] must be between 0 and 1")
            
        if x_range != ():
            self.ax.axvspan(x_range[0], x_range[1], alpha = shade_format['transparency'], color = shade_format['color'])

        if y_range != ():
            self.ax.axhspan(y_range[0], y_range[1], alpha = shade_format['transparency'], color = shade_format['color'])

        # Use draw instead of show so the graph gets updated it this is run after customzie_and_save
        plt.draw()

#------------------------
if __name__ == '__main__':
    pass
