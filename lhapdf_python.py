#! /usr/bin/python2

import lhapdf
import numpy
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
#import matplotlib.colors as colors


### SET DIRECTORIES ###

current_dir = os.getcwd() + "/"


### SET PDF ###

PDF_set_name = "high_x_var"
pdfs = lhapdf.mkPDFs(PDF_set_name)


### SET RANGE ###

xmin = 1e-2
xmax = 9.908496e-01
n_points = 100

xs = [x for x in numpy.logspace(numpy.log10(xmin), numpy.log10(xmax), n_points)]
#qs = [q for q in np.logspace(1, 4, 4)]
q = 5e+3


### GET VALUES FROM PDF GRIDS ###

#uV_list = numpy.empty([len(pdfs), len(xs)])
#dV_list = numpy.empty([len(pdfs), len(xs)])
#ubar_list = numpy.empty([len(pdfs), len(xs)])
#dbar_list = numpy.empty([len(pdfs), len(xs)])


#for PDF_mem in range(len(pdfs)):
    #for ix, x in enumerate(xs):
        #uV_list[PDF_mem, ix] = pdfs[PDF_mem].xfxQ(2, x, q)
        #ubar_list[PDF_mem, ix] = pdfs[PDF_mem].xfxQ(-2, x, q)
        #dV_list[PDF_mem, ix] = pdfs[PDF_mem].xfxQ(1, x, q)
        #dbar_list[PDF_mem, ix] = pdfs[PDF_mem].xfxQ(-1, x, q)


#dV_over_uV_list = numpy.empty([len(pdfs), len(xs)])
#dbar_minus_ubar_list = numpy.empty([len(pdfs), len(xs)])

dV_over_uV_list = numpy.empty([34, len(xs)])
dbar_minus_ubar_list = numpy.empty([34, len(xs)])

for PDF_mem in range(len(pdfs)):
    for ix, x in enumerate(xs):
        dV_over_uV_list[PDF_mem, ix] = pdfs[PDF_mem].xfxQ(1, x, q) / pdfs[PDF_mem].xfxQ(2, x, q)
        dbar_minus_ubar_list[PDF_mem, ix] = pdfs[PDF_mem].xfxQ(-1, x, q) - pdfs[PDF_mem].xfxQ(-2, x, q)


### SPECIAL COMBINATIONS ###

### VARIATION -0.3
PDF_mem += 1
for ix, x in enumerate(xs):
    dV_over_uV_list[PDF_mem, ix] = pdfs[7].xfxQ(1, x, q) / pdfs[19].xfxQ(2, x, q)
    dbar_minus_ubar_list[PDF_mem, ix] = pdfs[1].xfxQ(-1, x, q) - pdfs[13].xfxQ(-2, x, q)

### VARIATION +0.3
PDF_mem += 1
for ix, x in enumerate(xs):
    dV_over_uV_list[PDF_mem, ix] = pdfs[8].xfxQ(1, x, q) / pdfs[20].xfxQ(2, x, q)
    dbar_minus_ubar_list[PDF_mem, ix] = pdfs[2].xfxQ(-1, x, q) - pdfs[14].xfxQ(-2, x, q)

### VARIATION -0.5
PDF_mem += 1
for ix, x in enumerate(xs):
    dV_over_uV_list[PDF_mem, ix] = pdfs[9].xfxQ(1, x, q) / pdfs[21].xfxQ(2, x, q)
    dbar_minus_ubar_list[PDF_mem, ix] = pdfs[3].xfxQ(-1, x, q) - pdfs[15].xfxQ(-2, x, q)

### VARIATION +0.5
PDF_mem += 1
for ix, x in enumerate(xs):
    dV_over_uV_list[PDF_mem, ix] = pdfs[10].xfxQ(1, x, q) / pdfs[22].xfxQ(2, x, q)
    dbar_minus_ubar_list[PDF_mem, ix] = pdfs[4].xfxQ(-1, x, q) - pdfs[16].xfxQ(-2, x, q)

### VARIATION -1.0
PDF_mem += 1
for ix, x in enumerate(xs):
    dV_over_uV_list[PDF_mem, ix] = pdfs[11].xfxQ(1, x, q) / pdfs[23].xfxQ(2, x, q)
    dbar_minus_ubar_list[PDF_mem, ix] = pdfs[5].xfxQ(-1, x, q) - pdfs[17].xfxQ(-2, x, q)

### VARIATION +1.0
PDF_mem += 1
for ix, x in enumerate(xs):
    dV_over_uV_list[PDF_mem, ix] = pdfs[12].xfxQ(1, x, q) / pdfs[24].xfxQ(2, x, q)
    dbar_minus_ubar_list[PDF_mem, ix] = pdfs[6].xfxQ(-1, x, q) - pdfs[18].xfxQ(-2, x, q)


### MAXIMAL VARIATION ###

### MINUS (quark -1.0, antiquark +1.0)
PDF_mem += 1
for ix, x in enumerate(xs):
    dV_over_uV_list[PDF_mem, ix] = pdfs[11].xfxQ(1, x, q) / pdfs[23].xfxQ(2, x, q)
    dbar_minus_ubar_list[PDF_mem, ix] = pdfs[6].xfxQ(-1, x, q) - pdfs[18].xfxQ(-2, x, q)

### PLUS plus (quark +1.0, antiquark -1.0)
PDF_mem += 1
#print(PDF_mem)
for ix, x in enumerate(xs):
    dV_over_uV_list[PDF_mem, ix] = pdfs[12].xfxQ(1, x, q) / pdfs[24].xfxQ(2, x, q)
    dbar_minus_ubar_list[PDF_mem, ix] = pdfs[5].xfxQ(-1, x, q) - pdfs[17].xfxQ(-2, x, q)

tot_pdfs = PDF_mem + 1


### PLOT RESULTS ###

### matplotlib settings
matplotlib.rc("figure", dpi=600)
matplotlib.rc("savefig", pad_inches=0)
matplotlib.rc("xtick.minor", visible=True)
matplotlib.rc("ytick.minor", visible=True)
matplotlib.rc("axes.formatter", useoffset=False)

x_values = xs

#y_dV_over_uV_min = numpy.min(dV_over_uV_list)
#y_dV_over_uV_max = numpy.max(dV_over_uV_list)

#y_dbar_minus_ubar_min = numpy.min(dbar_minus_ubar_list)
#y_dbar_minus_ubar_max = numpy.max(dbar_minus_ubar_list)


y_dV_over_uV_min = 0.0
y_dV_over_uV_max = 0.9

y_dbar_minus_ubar_min = 0.00
y_dbar_minus_ubar_max = 0.04


#for PDF_mem in range(len(pdfs)):
for PDF_mem in range(tot_pdfs):
    
    plot_name = current_dir + "high_x_var_PDFs_" + str(PDF_mem) + ".pdf"

    fig, (ax1, ax2) = plt.subplots(2)
    
    fig.suptitle("PDF combination #" + str(PDF_mem))

    y_values = dV_over_uV_list[PDF_mem]
    ax1.set(xlabel=r"x", xscale="log", xlim=[xmin, xmax], ylabel=r"$d_V / u_V$", ylim=[y_dV_over_uV_min, y_dV_over_uV_max])
    #ax1.set_title("axes title")
    ax1.text(xmin * (1.1), 0.02, "Q = " + str(format(q, ".0f")) + " GeV", fontsize=10, verticalalignment="bottom", horizontalalignment="left")
    ax1.plot(x_values, y_values)

    y_values = dbar_minus_ubar_list[PDF_mem]
    ax2.set(xlabel=r"x", xscale="log", xlim=[xmin, xmax], ylabel=r"$\bar{d} - \bar{u}$", ylim=[y_dbar_minus_ubar_min, y_dbar_minus_ubar_max])
    ax2.text(xmax / (1.1), y_dbar_minus_ubar_max / (1.1), "Q = " + str(format(q, ".0f")) + " GeV", fontsize=10, verticalalignment="top", horizontalalignment="right")
    ax2.plot(x_values, y_values)

    plt.savefig(plot_name)
    
    plt.cla()
    plt.close(fig)


### MERGE RESULTS ###

### merge all plots in unique pdf
output_file_name = "high_x_var_PDFs_q_" + str(format(q/1e+3, ".0f")) + ".pdf"

file_list = ""
for PDF_mem in range(tot_pdfs):
    file_list += "high_x_var_PDFs_" + str(PDF_mem) + ".pdf"
    file_list += " "

command = "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dAutoRotatePages=/None -sOutputFile=" + output_file_name + " " + file_list

os.chdir(current_dir)
os.system(command)

### delete all sub plots files
for PDF_mem in range(tot_pdfs):
    file_name = "high_x_var_PDFs_" + str(PDF_mem) + ".pdf"
    if os.path.exists(current_dir + file_name):
            os.remove(current_dir + file_name)
