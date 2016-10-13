import pandas as pd
from bokeh.charts import Bar, output_file, show
from bokeh.layouts import column
from bokeh.models import (HoverTool, BoxSelectTool, BoxZoomTool, ResetTool,
                          PanTool, WheelZoomTool, ResizeTool)


class AlignViz2(object):

    def alignment_stats(self, input_path_file, output_path_file_csv,
                        output_path_file_html):
        df = pd.read_json(input_path_file, orient='index')
        df.to_csv(output_path_file_csv + '.csv',
                  index=True, sep='\t')

        samples = []
        alignment_length = []
        alignment_freq = []
        for index, row in df.iterrows():
            for key, value in row['stats_per_reference'].items():
                samples.extend([
                    index + '(' + key + ')'] * int(self.nr_items_align(value)))
                for k, v in value.items():
                    if k == 'alignment_length_and_freqs':
                        for keys, values in v.items():
                            alignment_length.extend([float(keys)])
                            alignment_freq.extend([values])

        data1 = {}
        data1['samples'] = samples
        data1['aligned read length'] = alignment_length
        data1['frequency'] = alignment_freq

        bar1 = Bar(
            data1, values='frequency', label='aligned read length',
            stack='samples', agg='sum',
            title="Alignment read length and frequency",
            legend='top_left', width=1200, bar_width=1.0,
            palette=['Blue', 'Aqua', 'SeaGreen', 'SpringGreen', 'Brown',
                     'Peru', 'Purple', 'Violet'],
            tools=[
                HoverTool(tooltips=[(
                    "Read length", "@x"), ("Frequency", "@y")]),
                PanTool(), BoxSelectTool(), BoxZoomTool(),
                WheelZoomTool(), ResizeTool(), ResetTool()])

        output_file(output_path_file_html + '.html')
        show(bar1)
        
    def nr_items_align(self, value):
        item_nr_align = []
        for k, v in value.items():
            if k == 'alignment_length_and_freqs':
                for keys, values in v.items():
                    item_nr_align.extend([keys])
        return(len(item_nr_align))

    def process_stats(self, input_path_file, output_path_file_csv,
                      output_path_file_html):
        df = pd.read_json(input_path_file, orient='index')
        df.to_csv(output_path_file_csv + '.csv', index=True, sep='\t')

        # plotting read processing(poly/single A removed or unmodified)
        conditions = []
        for row in df.iterrows():
            conditions.extend([
                'poly(A) removed', 'single(A) removed', 'unmodified'])

        samples1 = []
        for index, row in df.iterrows():
            samples1.extend([index] * 3)

        read_nr = []
        for index, row in df.iterrows():
            read_nr.extend([
                row['polya_removed'], row['single_a_removed'],
                row['unmodified']])

        data1 = {}
        data1['condition'] = conditions
        data1['sample'] = samples1
        data1['Nr. of reads'] = read_nr

        bar1 = Bar(
            data1, values='Nr. of reads', label='sample',
            stack='condition',
            agg='sum', title="Input read types", legend='top_right',
            palette=['darkseagreen', 'salmon', 'darkslateblue'],
            tools=[
                HoverTool(tooltips=[(
                    "Sample", "@x"), ("Nr of reads", "@y")]),
                PanTool(), BoxSelectTool(), BoxZoomTool(),
                WheelZoomTool(), ResizeTool(), ResetTool()])

        samples2 = []
        read_length = []
        frequency = []
        for index, row in df.iterrows():
            for key, value in row[
                    'read_length_after_processing_and_freq'].items():
                samples2.extend([index] * int(self.nr_items_pro(key)))
                read_length.extend([float(key)])
                frequency.extend([value])

        data2 = {}
        data2['samples'] = samples2
        data2['read length'] = read_length
        data2['frequency'] = frequency

        bar2 = Bar(
            data2, values='frequency', label='read length',
            stack='samples', agg='sum',
            title="Input read length and frequency",
            legend='top_left', palette=[
                'darkseagreen', 'salmon', 'darkslateblue', 'olive'],
            width=1200, bar_width=1.0,
            tools=[
                HoverTool(tooltips=[(
                    "Read length", "@x"), ("Frequency", "@y")]),
                PanTool(), BoxSelectTool(), BoxZoomTool(),
                WheelZoomTool(), ResizeTool(), ResetTool()])

        bar = column(bar1, bar2)
        output_file(output_path_file_html + '.html')
        show(bar)
        
    def nr_items_pro(self, key):
        item_nr_pro = []
        item_nr_pro.extend([key])
        return(len(item_nr_pro))
    

