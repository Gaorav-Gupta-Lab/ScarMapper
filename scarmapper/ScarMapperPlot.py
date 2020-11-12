"""
Draws main figures for ScarMapper

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2020
"""

import argparse
from contextlib import suppress
import pandas
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
import collections
import numpy
import Valkyries.Tool_Box as ToolBox

__author__ = 'Dennis A. Simpson'
__version__ = '0.2.0'

from Valkyries import Tool_Box


def scarmapperplot(args, datafile=None, sample_name=None, plot_data_dict=None, label_dict=None):

    try:
        datafile = args.DataFile
    except AttributeError:
        pass

    if sample_name:
        output_file_name = "{}_{}.{}".format(args.Job_Name, sample_name, args.FigureType)
    else:
        output_file_name = args.OutFile

    output_file = "{}{}".format(args.WorkingFolder, output_file_name)
    # df = pandas.read_csv("{}{}".format(args.WorkingFolder, datafile), sep='\t', skiprows=8)

    # Define colors for scar types and labels
    nhej_color = "royalblue"
    tmej_color = "red"
    insertion_color = "olive"
    non_mh_del_color = "mediumorchid"
    ins_del_color = "cyan"
    '''
    color_dict = \
        {'TMEJ': tmej_color, 'NHEJ': nhej_color, 'Non-MH Deletion': non_mh_del_color, 'Insertion': insertion_color,
         "Marker": "black", "TMEJ_Not-PolQ": "green", "Ins": ins_del_color}
    '''
    # fig, ax = plt.subplots()
    fig = plt.figure()
    fig.set_size_inches(8.5, 11.0)
    gs = fig.add_gridspec(4, hspace=0)
    ax = gs.subplots(sharex=True, sharey=True)

    # set background color of subplots
    ax[0].set_facecolor('whitesmoke')
    ax[1].set_facecolor('whitesmoke')
    ax[2].set_facecolor('whitesmoke')
    ax[3].set_facecolor('whitesmoke')

    # Turn off y-ticks
    ax[0].set_yticks([])
    ax[1].set_yticks([])
    ax[2].set_yticks([])
    ax[3].set_yticks([])

    # fig, (ax1, ax2) = plt.subplots(2, sharey='all', sharex='all')
    # plot_data_dict = build_plot_data_dict(df, color_dict)
    opacity = 1

    # [Bar Width, lft_del_plot_value, rt_del_plot_value, lft_ins_plot_value, rt_ins_plot_value, left ins width, right ins width, y-value]

    # With HR some of the scartypes are empty causing a plot error.
    scar_list = ['NHEJ', 'TsEJ', 'Non-MH Deletion', 'Insertion']
    for scartype in scar_list:
        if scartype not in plot_data_dict:
            plot_data_dict[scartype] = [0, 0, 0, 0, 0, 0, 0, 0]

    # Common NHEJ
    width_nhej = plot_data_dict['NHEJ'][0]
    x_lft_del_nhej = plot_data_dict['NHEJ'][1]
    x_rt_del_nhej = plot_data_dict['NHEJ'][2]
    x_lft_ins_nhej = plot_data_dict['NHEJ'][3]
    x_rt_ins_nhej = plot_data_dict['NHEJ'][4]
    l_ins_width_nhej = plot_data_dict['NHEJ'][5]
    r_ins_width_nhej = plot_data_dict['NHEJ'][6]
    y_nhej = plot_data_dict['NHEJ'][7]

    # Common TMEJ
    width_tmej = ""
    x_lft_del_tmej = ""
    x_rt_del_tmej = ""
    x_lft_ins_tmej = ""
    x_rt_ins_tmej = ""
    l_ins_width_tmej = ""
    r_ins_width_tmej = ""
    y_tmej = ""

    # Common TsEJ
    width_tmej = plot_data_dict['TsEJ'][0]
    x_lft_del_tmej = plot_data_dict['TsEJ'][1]
    x_rt_del_tmej = plot_data_dict['TsEJ'][2]
    x_lft_ins_tmej = plot_data_dict['TsEJ'][3]
    x_rt_ins_tmej = plot_data_dict['TsEJ'][4]
    l_ins_width_tmej = plot_data_dict['TsEJ'][5]
    r_ins_width_tmej = plot_data_dict['TsEJ'][6]
    y_tmej = plot_data_dict['TsEJ'][7]


    # Common non-MH deletion
    width_non_mh_del = plot_data_dict['Non-MH Deletion'][0]
    x_lft_del_non_mh_del = plot_data_dict['Non-MH Deletion'][1]
    x_rt_del_non_mh_del = plot_data_dict['Non-MH Deletion'][2]
    x_lft_ins_non_mh_del = plot_data_dict['Non-MH Deletion'][3]
    x_rt_ins_non_mh_del = plot_data_dict['Non-MH Deletion'][4]
    l_ins_width_non_mh_del = plot_data_dict['Non-MH Deletion'][5]
    r_ins_width_non_mh_del = plot_data_dict['Non-MH Deletion'][6]
    y_non_mh_del = plot_data_dict['Non-MH Deletion'][7]

    '''
    # Common TMEJ not PolQ
    width_not_polq = plot_data_dict['TMEJ_Not-PolQ'][0]
    x_lft_not_polq = plot_data_dict['TMEJ_Not-PolQ'][1]
    x_rt_not_polq = plot_data_dict['TMEJ_Not-PolQ'][2]
    y_tmej_not_polq = plot_data_dict['TMEJ_Not-PolQ'][5]
    color_not_polq = color_dict['TMEJ_Not-PolQ']
    '''
    # Common Insertions
    width_ins = plot_data_dict['Insertion'][0]
    x_lft_del_insertion = plot_data_dict['Insertion'][1]
    x_rt_del_insertion = plot_data_dict['Insertion'][2]
    x_lft_ins_insertion = plot_data_dict['Insertion'][3]
    x_rt_ins_insertion = plot_data_dict['Insertion'][4]
    l_ins_width_insertion = plot_data_dict['Insertion'][5]
    r_ins_width_insertion = plot_data_dict['Insertion'][6]
    y_ins = plot_data_dict['Insertion'][7]

    # Set the limits of the x-axis for all plots
    ax[3].set_xlim(plot_data_dict['Marker'][0], plot_data_dict['Marker'][1])

    '''Plot order is ascending order starting at top of page'''
    # Insertion plots
    ax[0].barh(y_ins, x_lft_del_insertion, width_ins, align='center', linewidth=0, color=insertion_color)
    ax[0].barh(y_ins, x_rt_del_insertion, width_ins, align='center', linewidth=0, color=insertion_color)
    ax[0].barh(y_ins, x_lft_ins_insertion, l_ins_width_insertion, alpha=opacity, align='center', linewidth=0,
               color=ins_del_color)
    ax[0].barh(y_ins, x_rt_ins_insertion, r_ins_width_insertion, alpha=opacity, align='center', linewidth=0,
               color=ins_del_color)

    ax[1].barh(y_non_mh_del, x_lft_del_non_mh_del, width_non_mh_del, align='center', linewidth=0,
               color=non_mh_del_color)
    ax[1].barh(y_non_mh_del, x_rt_del_non_mh_del, width_non_mh_del, align='center', linewidth=0,
               color=non_mh_del_color)
    ax[1].barh(y_non_mh_del, x_lft_ins_non_mh_del, l_ins_width_non_mh_del, alpha=opacity, align='center', linewidth=0,
               color=ins_del_color)
    ax[1].barh(y_non_mh_del, x_rt_ins_non_mh_del, r_ins_width_non_mh_del, alpha=opacity, align='center', linewidth=0,
               color=ins_del_color)

    ax[2].barh(y_tmej, x_lft_del_tmej, width_tmej, align='center', linewidth=0, color=tmej_color)
    ax[2].barh(y_tmej, x_rt_del_tmej, width_tmej, align='center', linewidth=0, color=tmej_color)
    ax[2].barh(y_tmej, x_lft_ins_tmej, l_ins_width_tmej, alpha=opacity, align='center', linewidth=0,
               color=ins_del_color)
    ax[2].barh(y_tmej, x_rt_ins_tmej, r_ins_width_tmej, alpha=opacity, align='center', linewidth=0,
               color=ins_del_color)

    ax[3].barh(y_nhej, x_lft_del_nhej, width_nhej, align='center', linewidth=0, color=nhej_color)
    ax[3].barh(y_nhej, x_rt_del_nhej, width_nhej, align='center', linewidth=0, color=nhej_color)
    ax[3].barh(y_nhej, x_lft_ins_nhej, l_ins_width_nhej, alpha=opacity, align='center', linewidth=0,
               color=ins_del_color)
    ax[3].barh(y_nhej, x_rt_ins_nhej, r_ins_width_nhej, alpha=opacity, align='center', linewidth=0,
               color=ins_del_color)

    '''
    ax[2].barh(y_tmej_not_polq, x_lft_not_polq, width_not_polq, align='center', linewidth=0, color=color_not_polq)
    ax[2].barh(y_tmej_not_polq, x_rt_not_polq, width_not_polq, align='center', linewidth=0, color=color_not_polq)
    ax[2].barh(y_tmej_not_polq, x_lft_not_polq, plot_data_dict['TMEJ_Not-PolQ'][3], alpha=opacity, align='center', 
    linewidth=0, color=ins_del_color)
    ax[2].barh(y_tmej_not_polq, x_rt_not_polq, plot_data_dict['TMEJ_Not-PolQ'][4], alpha=opacity, align='center', 
    linewidth=0, color=ins_del_color)
    '''

    # Add the center line to each plot
    ax[0].axvline(x=0, ls='-', lw=0.2, color='black')
    ax[1].axvline(x=0, ls='-', lw=0.2, color='black')
    ax[2].axvline(x=0, ls='-', lw=0.2, color='black')
    ax[3].axvline(x=0, ls='-', lw=0.2, color='black')
    # ax[4].axvline(x=0, ls='-', lw=0.25, color='black')

    # Add labels to X-axis and each plot.
    ax[3].set_xlabel('INDEL Size')

    ax[0].annotate('Insertion {}'.format(round(label_dict['Insertion'], 3)),
                   xy=(ax[0].get_xlim()[1] * -0.98, ax[0].get_ylim()[1] * 0.9), color=insertion_color, fontsize=14)
    ax[1].annotate('Non-MH Deletion {}'.format(round(label_dict['Non-MH Deletion'], 3)),
                   xy=(ax[1].get_xlim()[1] * -0.98, ax[1].get_ylim()[1] * 0.9), color=non_mh_del_color, fontsize=14)
    # ax[2].annotate('Non-PolQ MH', xy=(ax[2].get_xlim()[1] * -0.98, ax[2].get_ylim()[1] * 0.9))
    ax[2].annotate('TMEJ {}'.format(round(label_dict['TMEJ'], 3)),
                   xy=(ax[2].get_xlim()[1] * -0.98, ax[2].get_ylim()[1] * 0.9), color=tmej_color, fontsize=14)
    ax[3].annotate('NHEJ {}'.format(round(label_dict['NHEJ'], 3)),
                   xy=(ax[3].get_xlim()[1] * -0.98, ax[2].get_ylim()[1] * 0.9), color=nhej_color, fontsize=14)
    fig.suptitle(sample_name)

    if args.FigureType == "pdf":
        with PdfPages(output_file) as pdf:
            pdf.savefig(fig)
            plt.close()
    else:
        plt.savefig(output_file, dpi=800)
        plt.close()
        # plt.show()


def build_plot_data_dict(df, color_dict):
    """
    Sort data into a dictionary suitable for visualization.
    :param df:
    :param color_dict:
    :return:
    """
    scale_value = 1
    plot_data_dict = collections.defaultdict(list)

    for data_pair in df.values:
        ins_left_value = data_pair[0]
        ins_right_value = data_pair[0]
        if data_pair[2] != 0:
            ins_left_value = data_pair[0]*0.5
        if data_pair[3] != 0:
            ins_right_value = data_pair[0]*0.5

        if not data_pair[6] in plot_data_dict:
            plot_data_dict[data_pair[6]] = \
                [[data_pair[0]], [ins_left_value], [ins_right_value], [data_pair[2]], [data_pair[3]], [data_pair[4]],
                 [data_pair[5]], [color_dict[data_pair[6]]], [data_pair[0]*0.5]]
        else:
            plot_data_dict[data_pair[6]][0].append(data_pair[0])
            plot_data_dict[data_pair[6]][1].append(ins_left_value)
            plot_data_dict[data_pair[6]][2].append(ins_right_value)

            plot_data_dict[data_pair[6]][3].append(data_pair[2])
            plot_data_dict[data_pair[6]][4].append(data_pair[3])
            plot_data_dict[data_pair[6]][5].append(data_pair[4])

            plot_data_dict[data_pair[6]][6].append(data_pair[5])
            plot_data_dict[data_pair[6]][7].append(color_dict[data_pair[6]])
            count = len(plot_data_dict[data_pair[6]][0])

            if count > 1:
                previous = plot_data_dict[data_pair[6]][0][count - 2]
                plot_data_dict[data_pair[6]][8]\
                    .append(plot_data_dict[data_pair[6]][8][count - 2] + 0.0007 + (0.5*previous) + data_pair[0] * 0.5)

    return plot_data_dict


# This is here to run the module as a stand-alone.
if __name__ == '__main__':
    ToolBox.debug_messenger("Standing Alone")

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--options_file', action='store', dest='options_file', required=True,
                        help='File containing program parameters.')

    options_parser = ToolBox.options_file(parser)
    args = options_parser.parse_args()

    scarmapperplot(args)
