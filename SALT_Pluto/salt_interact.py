# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 15:14:33 2015

@author: throop
"""
# From http://stackoverflow.com/questions/9997869/interactive-plot-based-on-tkinter-and-matplotlib
# See http://stackoverflow.com/questions/16037494/python-code-is-it-comma-operator
#   for into on the extra comma (self.line,) - tuple unpacker: s, =x ==== s = x[0].

import Tkinter
import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
#import astropy

from HBTSaltPluto import HBTSaltPluto # See http://stackoverflow.com/questions/4534438/typeerror-module-object-is-not-callable
#from   astropy.convolution import convolve, Box1DKernel

import numpy as np

class App:

##########
# INIT CLASS
##########

    def __init__(self, master):
        # Create a container
        frame = Tkinter.Frame(master)

        # Set up how many lines and structure of the boxes
        
        self.num_lines = 22 # Total number of curves to plot. Should be even. 22 will include all of them.
        self.num_sets_listboxes = 0 # Number of curves using listboxes. These are faster to manipulate, but more real estate.
        self.num_sets_optionlists = self.num_lines - self.num_sets_listboxes

        self.path_settings = "/Users/throop/python/salt_interact_settings/"
        
        self.do_plot_inhibit = False

# Define the lists of dates
        
        NONE = '--'
        
        list_dates_pl_2015 = \
                        ['2015-05-03', '2015-05-04', '2015-05-06', 
                         '2015-06-30', '2015-07-03', '2015-08-02',
                         '2015-08-24', '2015-09-20']
                        
        list_dates_hd_2015 = \
                        ['2015-05-03', '2015-05-04', '2015-05-06',
                         '2015-06-29', '2015-06-30', '2015-07-03',
                         '2015-08-02', '2015-08-24']

        list_dates_pl_2014 = ['2014-06-21', '2014-06-23', '2014-06-24',
                              '2014-06-28', '2014-07-01', '2014-07-02',
                              '2014-08-16', '2014-08-18', '2014-08-22',
                              '2014-08-23', '2014-08-23', '2014-08-25',
                              '2014-08-26', '2014-09-05']
                              
        list_dates_hd_2014 = ['2014-06-21', '2014-08-12', '2014-08-16', 
                              '2014-08-26', '2014-09-05']
        
        list_dates_pl = list_dates_pl_2014  + list_dates_pl_2015
        list_dates_hd = list_dates_hd_2014  + list_dates_hd_2015
        
        self.list_dates_pl = list_dates_pl
        self.list_dates_hd = list_dates_hd

        xlim  = (3500, 9300) # Ideally we would set this from the data. 
        xdata = np.arange(3500,9300,1)
        ydata = np.sin(0.01 * np.arange(3500,9300,1))

        self.figsize = (15,4)

        self.binning = 1
        
        self.dy = 0.2 # Vertical offset between lines 

        self.ylim_pl_hd    = (0,    0.7 + self.dy * self.num_lines) # Default y axis range for Pl and HD data
        self.ylim_ratio = (0.25, 1.5 + self.dy * self.num_lines) # Default y axis range for ratio plot

        # Save the xrange of the dataset, so we can use it later for scaling
        # We could get this from the data itself, but for now, just hardcode it.

        self.xrange = (3500, 9300) 
        
        self.yrange_ratio  = self.ylim_ratio
        self.yrange_pl_hd  = self.ylim_pl_hd
        
        # Define the colors we will use
         
        self.colors = ['red', 'blue', 'green', 'purple', 'grey', 'orange', 'black', 'violet', 
                       'darkolivegreen', 'pink', 'cyan', 'lightgrey', 'mediumpurple', 'tan', 'sienna', 
                       'rosybrown', 'darkred', 'coral', 'thistle', 'navy', 'yellowgreen', 'darkkhaki',
                       'tomato', 'salmon', 'gold', 'lightblue', 'seagreen']

        # Create the sliders and labels
        
        self.scalex1  = Tkinter.Scale(master, from_=xlim[0], to=xlim[1], orient=Tkinter.HORIZONTAL, 
                                   command=self.do_scale) # Slider for central wavelength
        self.scalex2  = Tkinter.Scale(master, from_=1, to=40, resolution = 0.1, orient=Tkinter.HORIZONTAL, 
                                   command=self.do_scale) # Slider for zoom

        self.scaley1  = Tkinter.Scale(master, from_=0, to=self.num_lines * self.dy, resolution = 0.1, orient=Tkinter.HORIZONTAL, 
                                   command=self.do_scale) # Slider for central line (vertical)
        self.scaley2  = Tkinter.Scale(master, from_=1, to=5, resolution = 0.1, orient=Tkinter.HORIZONTAL, 
                                   command=self.do_scale) # Slider for zoom, Y direction
        
        self.scalebinning = Tkinter.Scale(master, from_=1, to=50, orient=Tkinter.HORIZONTAL, command=self.do_plot,
                                          label = 'Binning')
                           
        self.label1 = Tkinter.Label(master, text = 'Date Pluto', fg = self.colors[0])
        self.label2 = Tkinter.Label(master, text = 'Date HD',    fg = self.colors[0])

        self.label3 = Tkinter.Label(master, text = 'Date Pluto', fg = self.colors[1])
        self.label4 = Tkinter.Label(master, text = 'Date HD',    fg = self.colors[1])
   
        # Create the buttons
   
        self.button_save = Tkinter.Button(master, text = 'Save Settings [Name->]', command=self.save_settings)
        self.button_load = Tkinter.Button(master, text = 'Load Settings [Comment->]', command=self.load_settings)
        
        # Create the text entry boxes
        
        self.text_entry1 = Tkinter.Entry(master)
        self.text_entry1.insert(0, 'A default configuration')
        
        self.text_entry2 = Tkinter.Entry(master)
        self.text_entry2.insert(0, 'salt_set1')
        
        # Read in all data from disk
        
        self.read_all_data()
        
        # Set the xrange based on the data
        
        
# Set up Plot 1
        
# add_subplot: " kwargs are legal Axes kwargs. The Axes instance will be returned."
# http://matplotlib.org/api/figure_api.html
# http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.plot
     
        fig1 = Figure(figsize = self.figsize)

        self.lines1 = []         

        self.ax1 = fig1.add_subplot(1,1,1, 
                                    xlabel = 'Wavelength [A]', ylabel = 'Flux Pluto', label = 'Pluto', 
                                    ylim = self.ylim_pl_hd, xlim = xlim) # Return the axes

# Loop: Add six lines, each one specified by self.lines1 (that is, line in plot_panel #1)
        
        for i in range(self.num_lines):
          self.lines1.append( (self.ax1.plot(xdata, ydata,color = self.colors[i]))[0] )
                
        self.canvas1 = FigureCanvasTkAgg(fig1,master=master)
        self.canvas1.show()
        
# Set up Plot 2
       
        fig2 = Figure(figsize = self.figsize)
        self.ax2 = fig2.add_subplot(1,1,1, 
                                     xlabel = 'Wavelength [A]', ylabel = 'Flux HD', label = 'HD', 
                                     ylim = self.ylim_pl_hd, xlim = xlim) # Return the axes
        self.lines2 = []        
        for i in range(self.num_lines):
            self.lines2.append( (self.ax2.plot(xdata, ydata, color=self.colors[i]))[0] )
  
        self.canvas2 = FigureCanvasTkAgg(fig2,master=master)
        self.canvas2.show()
        
# Set up Plot 3
       
        fig3 = Figure(figsize = self.figsize)
        self.ax3 = fig3.add_subplot(1,1,1, xlabel = 'Wavelength [A]', ylabel = 'Flux Pluto/HD', label = 'Pluto/HD', 
                                    ylim = self.ylim_ratio, xlim = xlim) # Return the axes
        self.lines3 = []
        for i in range(self.num_lines):
            self.lines3.append( (self.ax3.plot(xdata, ydata, color=self.colors[i]))[0] )

        self.canvas3 = FigureCanvasTkAgg(fig3,master=master)
        self.canvas3.show()

# Set up the Pluto 'Listboxes' (ie, menus, but with all options explicitly shown)
# *** Disadvantage of these is that I have not debugged the code to set their values after reading from a file.
# Plus, they take up a lot more real estate. So it is common to disable these entirely (num_sets_listboxes = 0).
        
        self.lists_pl = []
        for i in range(self.num_sets_listboxes):
            self.lists_pl.append(Tkinter.Listbox(master, height=len(list_dates_pl), exportselection=0))   
            for item in list_dates_pl:
                self.lists_pl[-1].insert(Tkinter.END, item)
            self.lists_pl[-1].selection_set(0) # Set the default position 
 
# Set up the Pluto 'Listboxes' (ie, menus, but with all options explicitly shown)
   
        self.lists_hd = []
        for i in range(self.num_sets_listboxes):
            self.lists_hd.append(Tkinter.Listbox(master, height=len(list_dates_hd), exportselection=0))
            for item in list_dates_hd:
                self.lists_hd[-1].insert(Tkinter.END, item)
            self.lists_hd[-1].selection_set(0)

# Set up a bunch more 'OptionMenu' widgets, which are more compact ways of choosing more dates

        self.optlists_hd_var = []
        self.optlists_pl_var = []

        self.optlists_pl = []
        self.optlists_hd = []
 
# For each of the Pluto 'Option Lists' (aka menus): Create the variable for it, create the menu, and set color
       
        for i in range(self.num_sets_optionlists):
            index_default = np.mod(i + self.num_sets_listboxes,len(list_dates_pl))
            print repr(i) + ': default pl[' + repr(index_default) + ']'
            self.optlists_pl_var.append(Tkinter.StringVar(master))
            self.optlists_pl_var[-1].set(list_dates_pl[index_default]) # Set the value for the one we just added
            self.optlists_pl.append(Tkinter.OptionMenu(master, self.optlists_pl_var[-1], *list_dates_pl))
            self.optlists_pl[-1].config(width=12, bg = self.colors[i+self.num_sets_listboxes])

# For each of the HD 'Option Lists' (aka menus): Create the variable for it, create the menu, and set color

        for i in range(self.num_sets_optionlists):
            self.optlists_hd_var.append(Tkinter.StringVar(master))
            index_default = np.mod(i + self.num_sets_listboxes,len(list_dates_hd))
            print repr(i) + ': default hd[' + repr(index_default) + ']'
            self.optlists_hd_var[-1].set(list_dates_hd[index_default]) # Set the value for the one we just added
            self.optlists_hd.append(Tkinter.OptionMenu(master, self.optlists_hd_var[-1], *list_dates_hd))
            self.optlists_hd[-1].config(width=12, bg = self.colors[i+self.num_sets_listboxes])

# Put objects onto appropriate gridlines
        
        self.canvas1.get_tk_widget().grid(row=1, column=0, columnspan=6, sticky = 'nsew')
        self.canvas2.get_tk_widget().grid(row=2, column=0, columnspan=6, sticky = 'nsew')
        self.canvas3.get_tk_widget().grid(row=3, column=0, columnspan=6, sticky = 'nsew')
   
        self.text_entry1.grid(row=4, column=1, columnspan=1, sticky='ew')
        self.text_entry2.grid(row=5, column=1, columnspan=1, sticky='ew')
        
        self.scalex1.grid(row=4, column=3, columnspan=1, sticky = 'ew')
        self.scalex2.grid(row=5, column=3, columnspan=1, sticky = 'ew')

        self.scaley1.grid(row=4, column=4, columnspan=1, sticky = 'ew')
        self.scaley2.grid(row=5, column=4, columnspan=1, sticky = 'ew')
        
        self.scalebinning.grid(row=4, column=2)
        
        self.button_load.grid(row=4, column=0)
        self.button_save.grid(row=5, column=0)
        
        self.label1.grid(row=8, column=0)
        self.label2.grid(row=8, column=1)
        self.label3.grid(row=8, column=2)
        self.label4.grid(row=8, column=3)

# Put the Listboxes on a grid
        if (self.num_sets_listboxes == 2):
            self.lists_pl[0].grid(row=10, column=0, sticky = 'n')
            self.lists_hd[0].grid(row=10, column=1, sticky = 'n')
            self.lists_pl[1].grid(row=10, column=2, sticky = 'n')
            self.lists_hd[1].grid(row=10, column=3, sticky = 'n')

        for i in range(self.num_sets_optionlists): # 0, 1, 2, 3  if 4 lines left              
            self.optlists_pl[i].grid(row=12 + np.int(i/3), column=np.mod(i,3)*2)
            self.optlists_hd[i].grid(row=12 + np.int(i/3), column=np.mod(i,3)*2+1)
#            self.optlists_hd[i].grid(row=12 + np.int(i/3), column=np.mod(i,3)*3+2)
              
# Set weights for various rows and columns

        frame.grid_rowconfigure(0, weight=1)
        frame.grid_rowconfigure(1, weight=1)
        frame.grid_rowconfigure(2, weight=1)
        frame.grid_columnconfigure(1, weight=1)        
        frame.grid_columnconfigure(5, weight=1)

        frame.grid_rowconfigure(3, weight=1)
        frame.grid_rowconfigure(4, weight=1)
        frame.grid_rowconfigure(5, weight=1)

# Test: Try to make another sub-grid to the right
# XXX not working -- not sure why not

        t = Tkinter.Frame(frame)
        t1 = Tkinter.Label(t, text='Text1')
        t2 = Tkinter.Label(t, text='Text2')
        t1.grid(column=1, row=1)
        t2.grid(column=1,  row=3)
        t.grid(row=4, column=2, sticky='nswe')
               
# Set up callbacks

# Bind to events. The normal way in Python is to have a callback, but for events (like listbox), we use bind.

        for l in self.lists_pl:
            l.bind('<<ListboxSelect>>', self.do_plot)
        for l in self.lists_hd:
            l.bind('<<ListboxSelect>>', self.do_plot)
             
        for v in self.optlists_pl_var:
            v.trace("w", self.do_plot)
        for v in self.optlists_hd_var:
            v.trace("w", self.do_plot)
            
# Call do_plot, to make a plot of current settings.

        self.do_plot()

##########        
### END OF __INIT__
##########

##########
# Rescale / zoom all of the plots
##########

# NB: Not well documented, but since this is called by dragging a slider, it is passed two args, not one.
  
    def do_scale(self, event):

        # Get the slider positions, for the four controls

        xcenter = self.scalex1.get() # Retrieve central value (in angstroms)
        xzoom   = self.scalex2.get() # Retreive zoom level (1 = unzoomed; 10 = 10x zoomed in)
        
        ycenter= self.scaley1.get()
        yzoom  = self.scaley2.get()

        # Get the xrange of the         
        x0      = self.xrange[0]    # Xmin of the full dataset 
        x1      = self.xrange[1]    # Xmax of the full dataset
        x0p     = xcenter - (xcenter - x0) / xzoom  # x0' : xmin of scaled dataset
        x1p     = xcenter + (x1 - xcenter) / xzoom  # x1' : xmax of scaled dataset

        # Set yrange for Pluto and HD curves

        y0            = self.yrange_pl_hd[0]
        y1            = self.yrange_pl_hd[1]
        y0p_pl_hd     = ycenter - (ycenter-y0)/yzoom
        y1p_pl_hd     = ycenter + (y1-ycenter)/yzoom
 
        y0            = self.yrange_ratio[0]
        y1            = self.yrange_ratio[1]
        y0p_ratio     = ycenter - (ycenter-y0)/yzoom
        y1p_ratio     = ycenter + (y1-ycenter)/yzoom       
        
        self.ax1.set_xlim(x0p, x1p)
        self.canvas1.draw()
        self.ax2.set_xlim(x0p, x1p)
        self.canvas2.draw()
        self.ax3.set_xlim(x0p, x1p)
        self.canvas3.draw()
        
        self.ax1.set_ylim(y0p_pl_hd, y1p_pl_hd)
        self.canvas1.draw()
        self.ax2.set_ylim(y0p_pl_hd, y1p_pl_hd)
        self.canvas2.draw()
        self.ax3.set_ylim(y0p_ratio, y1p_ratio)
        self.canvas3.draw()
        
##########
# Read in all of the data from disk. We do this ahead of time, rather than read when needed.
##########

    def read_all_data(self):
        
        print "Starting to read files..."
        self.flux_pl       = []
        self.wavelength_pl = []
        self.flux_hd       = []
        self.wavelength_hd = []
        
        for date in self.list_dates_pl:
            handle = HBTSaltPluto('spect_pluto_merged_' + date + '.txt')
            (w,f) = handle.readspect()
            self.wavelength_pl.append(w)
            self.flux_pl.append(f)

        for date in self.list_dates_hd:
            handle = HBTSaltPluto('spect_hd_merged_' + date + '.txt')
            (w,f) = handle.readspect()
            self.wavelength_hd.append(w)
            self.flux_hd.append(f)
            
        print "Finishing reading files..."

##########        
# Render all of the plots
##########

# This do_plot() function gets called a) explicitly, b) as a Listbox event callback, and c) as a Menuoption trace.
# These pass 1, 2, and 5 variables, respectively (mostly the calling event info, etc). But I don't care: in each
# case I just want to plot it. However, to avoid having problems, I must set up the func so it can accept up to 5
# arguments -- even if I ignore them.
# # http://stackoverflow.com/questions/6554805/getting-a-callback-when-a-tkinter-listbox-selection-is-changed
        
    def do_plot(self, junk=0, junk1=0, junk2=0, junk3=0):    
        
        # do_plot is called automatically when variables are changed. But sometimes (like when loading) we want to
        # inhibit the auto-redraw. So in this case, we check a flag, and if set, return immediately.
        
        if (self.do_plot_inhibit): 
            return
            
        dates_pl_selected = [] # This is the list of selected PLuto dates: from listboxes, and from optionboxes
        dates_hd_selected = []

# Look up the current date of each Listbox

        for l in self.lists_pl:
            dates_pl_selected.append(l.get(l.curselection()))
            
        for l in self.lists_hd:
            dates_hd_selected.append(l.get(l.curselection()))

        for l in self.optlists_pl_var:
            dates_pl_selected.append(l.get())

        for l in self.optlists_hd_var:
            dates_hd_selected.append(l.get())
        
        binning = self.scalebinning.get()
        
# Draw the lines for plot #1 (Pluto)
# We skip over the missing data using NaN or None. NaN is probably better. 
#   http://stackoverflow.com/questions/17534106/what-is-the-difference-between-nan-and-none
        
        i=0
        for line in self.lines1:
            index_pl = self.list_dates_pl.index(dates_pl_selected[i]) # Find index of dates_pl_selected within list_dates_pl
            ydata = np.copy(self.flux_pl[index_pl])
            ydata[ydata == 0] = None                  # Convert zero -> None so they will not plot. Cool! Not same as NaN I think.
            ydata += self.dy*(self.num_lines - i - 1) # Offset the lines from each other
            ydata = np.squeeze(ydata)                 # Remove any 1D axes from array
            ydata = np.convolve(ydata, 1./binning + np.zeros(binning), mode = 'same')
            line.set_ydata(ydata)
            line.set_xdata(self.wavelength_pl[index_pl])
            i+=1
        self.canvas1.draw()

# Draw the lines for plot #2 (HD)
        
        i=0
        for line in self.lines2:
            index_hd = self.list_dates_hd.index(dates_hd_selected[i])
            ydata = np.copy(self.flux_hd[index_hd])
            ydata[ydata == 0] = None
            ydata += self.dy*(self.num_lines - i - 1)
            ydata = np.squeeze(ydata)
            ydata = np.convolve(ydata, 1./binning + np.zeros(binning), mode = 'same')
            line.set_ydata(ydata)
            line.set_xdata(self.wavelength_hd[index_hd])
            i+=1
        self.canvas2.draw()

# Draw the lines for plot #3 (Ratio)
        
        i=0
        for line in self.lines3:
            index_hd = self.list_dates_hd.index(dates_hd_selected[i])
            index_pl = self.list_dates_pl.index(dates_pl_selected[i])
            ydata = np.copy(self.flux_pl[index_pl] / self.flux_hd[index_hd])
            ydata[ydata == 0] = None
            ydata += self.dy*(self.num_lines - i - 1)
            ydata = np.squeeze(ydata)
            ydata = np.convolve(ydata, 1./binning + np.zeros(binning), mode = 'same')
            line.set_xdata(self.wavelength_pl[index_pl])
            line.set_ydata(ydata)
            i+=1
        self.canvas3.draw()

##########
# SAVE SETTINGS.
# Write the current settings (date for each Pluto and HD) to a text file
# We should also allow a one-line comment, and a specified filename
##########

    def save_settings(self):
        dates_hd = []
        dates_pl = []
        
        # First read all of the current settings: Pluto and HD dates, from listbox
        
        for l in self.lists_pl:
            dates_pl.append(l.get(l.curselection()))
            
        for l in self.lists_hd:
            dates_hd.append(l.get(l.curselection()))

        # Then read all of the current settings: Pluto and HD dates, from optlist

        for l in self.optlists_pl_var:
            dates_pl.append(l.get())

        for l in self.optlists_hd_var:
            dates_hd.append(l.get())
        
        # Set the output filename
        
        filename_out = self.path_settings + self.text_entry2.get()
        
        # Prepare the output array
        
        arr_out = []
        for i in range(len(dates_pl)):
            arr_out.append(repr(i) + ', ' + dates_pl[i] + ', ' + dates_hd[i])
            print repr(i) + ": " + arr_out[-1]
                    
        header = self.text_entry1.get() + "\n"
        header = header + "Num, Pluto, HD, " + repr(self.num_lines) + ', ' + repr(self.num_sets_listboxes) + \
            ', ' + repr(self.num_sets_optionlists)

        # And write the file to disk
        
        np.savetxt(filename_out, np.array(arr_out), header=header, fmt="%s")
        print 'Wrote: ' + filename_out

##########
# Load Settings
##########
        
    def load_settings(self):  
            
        filename_in = self.path_settings + self.text_entry2.get()
        
        f = open(filename_in)
        lines = f.readlines()
        f.close()

#        d = np.loadtxt(filename_in, delimiter=',', dtype = 'S50')
        header_comment = (lines[0])[2:-1] # remove the hashtag and newline
        header_columns = (lines[1])[:-1]
        
        lines = lines[2:]
        
        print "Reading: " + filename_in
        
        print "d[1] = " + repr(d[1])
            
        (index, sa, sb, num_lines, num_sets_listboxes, num_sets_optionlists) = header_columns.split(',')
        num_lines = int(num_lines)
        num_sets_listboxes = int(num_sets_listboxes)
        num_sets_optionlists = int(num_sets_optionlists)
        
        if (self.num_lines == num_lines) and (self.num_sets_listboxes == num_sets_listboxes):
            for i in range(num_sets_listboxes):
                print "Unpacking line: " + lines[i]
                lines[i].strip()    # Remove whitespace, \n, etc.
                (num, date_pl, date_hd) = lines[i].strip().split(',')
                date_hd = date_hd.strip()
                date_pl = date_pl.strip()
                print "listbox: attempting to set date_pl = '" + date_pl + "'; date_hd = '" + date_hd + "'"   
                # For Listbox, we set it by passing the *index* to the list iteem

                self.lists_pl[i].selection_set(self.list_dates_pl.index(date_pl))
                self.lists_hd[i].selection_set(self.list_dates_hd.index(date_hd))
                
            for i in range(num_sets_optionlists):
                (num, date_pl, date_hd) = lines[i+num_sets_listboxes].strip().split(', ')

                # For optionbox, we set it by passing the *value* of the list item
                # There is a small bug here, that setting the value (by design) calls the update routine.
                # I actually don't want this to happen!
                
                self.do_plot_inhibit = True
                
                self.optlists_hd_var[i].set(date_hd)
                self.optlists_pl_var[i].set(date_pl)
                
                self.do_plot_inhibit = False
            self.text_entry1.delete(0, Tkinter.END)     # Clear current setting
            self.text_entry1.insert(0, header_comment) # Set the new value
            
            print "Loaded settings."
#                print "setting optionlist " + repr(i) + ' -> ' + date_hd
            self.do_plot() # Refresh screen
        else:
            print "Incompatible format -- can't load"
                
                
                
##########
# Start of main program
##########
        
date = '2015-08-24'
root = Tkinter.Tk()
app = App(root)

#w = 500 # width for the Tk root
#h = 1050 # height for the Tk root

# get screen width and height
ws = root.winfo_width() # width of the screen
hs = root.winfo_height() # height of the screen

# calculate x and y coordinates for the Tk root window
x = 20 # (ws/2) - (w/2)
y = 20 # (hs/2) - (h/2)

# set the dimensions of the screen 
# and where it is placed
#root.geometry('%dx%d+%d+%d' % (1160, 1400, 30, 30))

root.mainloop()
