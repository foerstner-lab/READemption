import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from bokeh.charts import Bar
from bokeh.layouts import column
from bokeh.models import (HoverTool, BoxZoomTool, ResetTool)
from bokeh.resources import CDN
from bokeh.embed import file_html


class AlignViz(object):

    def processing_viz(self, input_process, output_path):
        with open(input_process) as process_data:
            process_data = json.load(process_data)
        # Manage plotting of processing stats
        self._get_plotting_data_bokeh_his1_his2(process_data, output_path)
        self._get_plotting_data_pp_his1(process_data, output_path)
        self._get_plotting_data_pp_his2(process_data, output_path)
        
    def alignment_viz(self, input_align, output_path):
        with open(input_align) as align_data:
            alignment_data = json.load(align_data)
        # Manage plotting of alignment data
        self._get_plotting_data_pp_his3(alignment_data, output_path)
        self._get_plotting_data_bokeh_his3(alignment_data, output_path)

    def alignment_processing_overview(
            self, input_process, input_align, output_path):
        with open(input_align) as align_data:
            alignment_data = json.load(align_data)
        with open(input_process) as process_data:
            process_data = json.load(process_data)
        # Manage plotting alignment overview
        self._get_plotting_data_both_his4(
            process_data, alignment_data, output_path)

    def _get_plotting_data_bokeh_his1_his2(self, process_data, output_path):
        data_list_bokeh_his1 = []
        data_list_bokeh_his2 = []
        data_list_bokeh_his1.append(
            self._create_data_set_bokeh_his1(process_data))
        for condition, sub_dict in process_data.items():
            data_list_bokeh_his2.append(self._create_data_set_bokeh_his2(
                condition, sub_dict))
        self._generate_his1_his2_bokeh(
            data_list_bokeh_his1, data_list_bokeh_his2, output_path)

    def _create_data_set_bokeh_his1(self, process_data):
        read_type = []
        samples = []
        read_nr = []
        for condition, sub_dict in process_data.items():
            read_type.extend([
                'poly(A) removed', 'single(A) removed', 'unmodified'])
            samples.extend([condition] * 3)
            read_nr.extend([
                sub_dict['polya_removed'], sub_dict['single_a_removed'],
                sub_dict['unmodified']])
        data_bokeh_his1 = {}
        data_bokeh_his1['read_type'] = read_type
        data_bokeh_his1['condition'] = samples
        data_bokeh_his1['No Reads [nt]'] = read_nr
        return(data_bokeh_his1)

    def _create_data_set_bokeh_his2(self, condition, sub_dict):
        read_length_bokeh_his2 = []
        read_freq_bokeh_his2 = []
        for length, freq in sub_dict[
                'read_length_after_processing_and_freq'].items():
            read_length_bokeh_his2.append(float(length))
            read_freq_bokeh_his2.append(freq)
        if not read_length_bokeh_his2:
            read_length_bokeh_his2.append(0)
        if not read_freq_bokeh_his2:
            read_freq_bokeh_his2.append(0)
        data_bokeh_his2 = {}
        data_bokeh_his2['samples'] = ([condition] * int(
            self._count_his2(sub_dict)))
        if not data_bokeh_his2['samples']:
            data_bokeh_his2['samples'].append(0)
        data_bokeh_his2[
            'Read length [nt] after processing'] = read_length_bokeh_his2
        data_bokeh_his2['No Reads [nt]'] = read_freq_bokeh_his2
        return(data_bokeh_his2)

    def _count_his2(self, sub_dict):
        nr_of_items = []
        for length, freq in sub_dict[
                'read_length_after_processing_and_freq'].items():
                    nr_of_items.extend([length])
        return(len(nr_of_items))
    
    def _generate_his1_his2_bokeh(
            self, data_list_bokeh_his1, data_list_bokeh_his2, output_path):
        bars = []
        for data_set in data_list_bokeh_his1:
            bar = Bar(
                data_set, values='No Reads [nt]', label='condition',
                stack='read_type', agg='sum', legend='top_right',
                title='Total No of Input Reads after Processing', palette=[
                    'DarkBlue', 'DarkCyan', 'Black'], tools=[
                        HoverTool(tooltips=[("No of reads", "@height")]),
                        BoxZoomTool(), ResetTool()])
            bar.legend.background_fill_alpha = 0.5
            bars.append(bar)
        for data_set in data_list_bokeh_his2:
            bar = Bar(
                data_set, values='No Reads [nt]', stack='samples',
                label='Read length [nt] after processing', agg='sum',
                legend='top_left', width=1500, bar_width=1.0,
                color='DarkSlateGray', tools=[
                    HoverTool(tooltips=[(
                        "Read length", "@x"), ("No of reads", "@height")])],
                title='Length & Frequency of Processed Reads')
            bar.legend.background_fill_alpha = 0.5
            bars.append(bar)
        plots = column(*bars)
        html = file_html(plots, CDN, 'Processing Stats')
        with open(
                output_path + '/Process Overview + Lengths & Freqs.html',
                'w') as output_bokeh_his1_his2:
            output_bokeh_his1_his2.write(html)
            
    def _get_plotting_data_pp_his1(self, process_data, output_path):
        poly_A_rem = []
        single_A_rem = []
        unmodified = []
        condition_list = []
        for condition, sub_dict in process_data.items():
            poly_A_rem.append(sub_dict['polya_removed'])
            unmodified.append(sub_dict['unmodified'])
            single_A_rem.append(sub_dict['single_a_removed'])
            condition_list.append(condition)
        self._generate_his1_pp(poly_A_rem, single_A_rem, unmodified,
                               condition_list, output_path)

    def _generate_his1_pp(self, poly_A_rem, single_A_rem,
                          unmodified, condition_list, output_path):
        with PdfPages(
                '%s/Read Processing Overview.pdf' % (output_path)) as pdf:
            matplotlib.style.use('ggplot')
            bar_width = 0.5
            pos = np.arange(len(condition_list))
            fig, ax = plt.subplots()
            plt.bar(pos, poly_A_rem, bar_width, label='Poly(A) removed',
                    color='b')
            plt.bar(pos, single_A_rem, bar_width, label='Single(A) removed',
                    color='w', edgecolor='b', hatch='//', bottom=poly_A_rem)
            plt.bar(pos, unmodified, bar_width, label='Unmodified', color='k',
                    bottom=np.array(single_A_rem) + np.array(poly_A_rem))
            plt.ylabel('No reads [nt]')
            plt.title('Total No of Input Reads after Processing')
            plt.xticks(pos, tuple(condition_list), rotation=10)
            leg = plt.legend(fancybox=True, loc='best')
            leg.get_frame().set_alpha(0.5)
            pdf.savefig()
            plt.close()
    
    def _get_plotting_data_pp_his2(self, process_data, output_path):
        pp = PdfPages(
            '%s/Length & Frequency of Processed Reads.pdf' % (output_path))
        for condition, sub_dict in process_data.items():
            title = condition
            fig = self._generate_his2_pp(title, sub_dict)
            pp.savefig()
            plt.close(fig)
        pp.close()

    def _generate_his2_pp(self, title, sub_dict):
        fig = plt.figure()
        length_list = [int(length) for length in sub_dict[
            'read_length_after_processing_and_freq'].keys()]
        if not length_list:
            lengths = np.array([0])
        else:
            lengths = np.array(length_list)
        freq_list = [int(freq) for freq in sub_dict[
            'read_length_after_processing_and_freq'].values()]
        if not freq_list:
            freqs = np.array([0])
        else:
            freqs = np.array(freq_list)
        ax = fig.add_subplot(111)
        plt.title(title)
        plt.xlabel('Read length [nt]')
        plt.ylabel('No of reads')
        font = {'size': 8}
        matplotlib.rc('font', **font)
        ax.xaxis.set_ticks_position('bottom')
        plt.xticks(np.arange(0, max(lengths)+1, 10.0))
        plt.bar(lengths, freqs, align='center', color='black',
                edgecolor='none')
        plt.xlim([0, max(lengths)+1])
        return(fig)

    def _get_plotting_data_bokeh_his3(self, alignment_data, output_path):
        data_list_bokeh_his3 = []
        for condition, sub_dict in alignment_data.items():
            for replicon, values in (sub_dict['stats_per_reference']).items():
                data_list_bokeh_his3.append(
                    self._create_data_set_bokeh_his3(
                        condition, replicon, sub_dict, values))
        self._generate_his3_bokeh(data_list_bokeh_his3, output_path)

    def _create_data_set_bokeh_his3(
            self, condition, replicon, sub_dict, values):
        alignment_length = []
        alignment_freq = []
        for length, freq in values['alignment_length_and_freqs'].items():
            alignment_length.append(float(length))
            alignment_freq.append(freq)
        if not alignment_length:
            alignment_length.append(0)
        if not alignment_freq:
            alignment_freq.append(0)
        data_bokeh_his3 = {}
        data_bokeh_his3['sample'] = list([
            condition + ' (' + replicon + ')'] * int(self._count_his3(values)))
        if not data_bokeh_his3['sample']:
            data_bokeh_his3['sample'].append(0)
        data_bokeh_his3['Length [nt] of aligned reads'] = alignment_length
        data_bokeh_his3['No Reads [nt]'] = alignment_freq
        return(data_bokeh_his3)

    def _generate_his3_bokeh(self, data_list_bokeh_his3, output_path):
        bars = []
        for data_set in data_list_bokeh_his3:
            bar = Bar(
                data_set, values='No Reads [nt]', stack='sample',
                label='Length [nt] of aligned reads', agg='sum',
                legend='top_left', width=1500, bar_width=1.0,
                color='Black', tools=[
                    HoverTool(tooltips=[(
                        "Read length", "@x"), ("No of reads", "@height")])],
                title='Length & Frequency of Aligned Reads')
            bar.legend.background_fill_alpha = 0.5
            bars.append(bar)
        plots = column(*bars)
        html = file_html(plots, CDN, 'Alignment Stats')
        with open(output_path + '/Length & Frequency of Aligned Reads.html',
                  'w') as output_bokeh_his3:
            output_bokeh_his3.write(html)

    def _count_his3(self, values):
        nr_of_items = []
        for read_type, value in values.items():
            if read_type == 'alignment_length_and_freqs':
                for length, freq in value.items():
                    nr_of_items.extend([length])
        return(len(nr_of_items))

    def _get_plotting_data_pp_his3(self, alignment_data, output_path):
        pp = PdfPages(
            '%s/Length & Frequency of Aligned Reads.pdf' % (output_path))
        for condition, sub_dict in alignment_data.items():
            for replicon, values in sub_dict['stats_per_reference'].items():
                title = (condition + '(' + replicon + ')')
                fig = self._generate_his3_pp(title, values)
                pp.savefig()
                plt.close(fig)
        pp.close()
        
    def _generate_his3_pp(self, title, values):
        fig = plt.figure()
        length_list = [int(length) for length in values[
            'alignment_length_and_freqs'].keys()]
        if not length_list:
            lengths = np.array([0])
        else:
            lengths = np.array(length_list)
        freq_list = [int(freq) for freq in values[
            'alignment_length_and_freqs'].values()]
        if not freq_list:
            freqs = np.array([0])
        else:
            freqs = np.array(freq_list)
        ax = fig.add_subplot(111)
        plt.title(title)
        plt.xlabel('Read length [nt]')
        plt.ylabel('No of reads')
        font = {'size': 8}
        matplotlib.rc('font', **font)
        ax.xaxis.set_ticks_position('bottom')
        plt.xticks(np.arange(0, max(lengths)+1, 10.0))
        plt.bar(lengths, freqs, align='center', color='k',
                edgecolor='none')
        plt.xlim([0, max(lengths)+1])
        return(fig)

    def _get_plotting_data_both_his4(
            self, process_data, alignment_data, output_path):
        data_list_bokeh_his4 = []
        alignment_data_mod = {}
        for replicon, sub_dict in alignment_data.items():
            del(sub_dict["stats_per_reference"])
            alignment_data_mod[replicon] = sub_dict["stats_total"]
        both_dicts = [process_data, alignment_data_mod]
        files_merged = {k: [d[k] for d in both_dicts] for k in both_dicts[0]}
        data_list_bokeh_his4.append(
            self._create_data_set_bokeh_his4(files_merged))
        self._generate_histogram_bokeh_his4(data_list_bokeh_his4, output_path)
        self._get_plotting_data_pp_his4(files_merged, output_path)

    def _create_data_set_bokeh_his4(self, files_merged):
        read_type = []
        read_types = "No input reads (total)", "No input reads (long enough)", "No aligned reads", "No uniquely aligned reads"
        read_type.extend(read_types * int(self._count_his4(files_merged)))
        conditions = []
        read_nrs = []
        for library, dict_list in files_merged.items():
            conditions.extend([library] * 4)
            for dicts in dict_list:
                for read, value in dicts.items():
                    if read == 'total_no_of_reads':
                        read_nrs.append(value)
                    if read == 'long_enough':
                        read_nrs.append(value)
                    if read == 'no_of_aligned_reads':
                        read_nrs.append(value)
                    if read == 'no_of_uniquely_aligned_reads':
                        read_nrs.append(value)
        read_nrs_sorted = []
        for index in range(0, len(read_nrs), 4):
            sub_list = sorted(read_nrs[index:index+4], reverse=True)
            read_nrs_sorted = read_nrs_sorted + sub_list
        data_bokeh_his4 = {}
        data_bokeh_his4['Library'] = conditions
        data_bokeh_his4['read_type'] = read_type
        data_bokeh_his4['No reads [nt]'] = read_nrs_sorted
        return(data_bokeh_his4)
    
    def _generate_histogram_bokeh_his4(
            self, data_list_bokeh_his4, output_path):
        for data_set in data_list_bokeh_his4:
            bar = Bar(data_set, values='No reads [nt]', label='Library',
                      group='read_type', title='Read Alignment Overview',
                      palette=[
                          'DarkBlue', 'DarkSlateGray', 'Black', 'DarkCyan'],
                      legend='top_right', tools=[
                          HoverTool(tooltips=[("No of reads", "@height")])])
            bar.legend.background_fill_alpha = 0.5
        html = file_html(bar, CDN, 'Mapping Overview')
        with open(output_path +
                  '/Alignment Overview.html', 'w') as output_bokeh_his4:
            output_bokeh_his4.write(html)
    
    def _get_plotting_data_pp_his4(self, files_merged, output_path):
        pp_conditions = []
        pp_total_read_nr = []
        pp_long_enough_reads = []
        pp_aligned_reads = []
        pp_uniq_aligned_reads = []
        for condition, dict_list in files_merged.items():
            pp_conditions.append(condition)
            for dicts in dict_list:
                for read, value in dicts.items():
                    if read == 'total_no_of_reads':
                        pp_total_read_nr.append(value)
                    if read == 'long_enough':
                        pp_long_enough_reads.append(value)
                    if read == 'no_of_aligned_reads':
                        pp_aligned_reads.append(value)
                    if read == 'no_of_uniquely_aligned_reads':
                        pp_uniq_aligned_reads.append(value)
        self._generate_his4_pp(pp_conditions, pp_total_read_nr,
                               pp_long_enough_reads, pp_aligned_reads,
                               pp_uniq_aligned_reads, output_path)

    def _generate_his4_pp(self, pp_conditions, pp_total_read_nr,
                          pp_long_enough_reads, pp_aligned_reads,
                          pp_uniq_aligned_reads, output_path):
        with PdfPages('%s/Alignment Overview.pdf' % (output_path)) as pdf:
            matplotlib.style.use('ggplot')
            bar_width = 0.2
            pos1 = np.arange(len(pp_conditions))
            fig, ax = plt.subplots()
            plt.bar(pos1, pp_total_read_nr, bar_width, linewidth=0.2,
                    label='No input reads (total)', color='k')
            plt.bar((pos1 + bar_width), pp_long_enough_reads, bar_width,
                    linewidth=0.1, label='No input reads (long enough)',
                    color='b')
            plt.bar((pos1 + bar_width * 2), pp_aligned_reads, bar_width,
                    linewidth=0.1, label='No aligned reads', color='c')
            plt.bar((pos1 + bar_width * 3), pp_uniq_aligned_reads, bar_width,
                    color='w', edgecolor='c', hatch='//', linewidth=0.1,
                    label='No uniquely aligned reads')
            plt.xticks(pos1, pp_conditions, rotation=10)
            leg = plt.legend(fancybox=True, loc='best')
            leg.get_frame().set_alpha(0.5)
            plt.ylabel('No reads [nt]')
            plt.title('Read Alignment Overview')
            pdf.savefig()
            plt.close()

    def _count_his4(self, new_dic):
        Nr_of_conditions = []
        for condition, sub_dict in new_dic.items():
            Nr_of_conditions.append(condition)
        return(len(Nr_of_conditions))

