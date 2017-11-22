from reademptionlib.visdeseq2 import DESeq2Vis
from reademptionlib.helpers import Helpers
from reademptionlib.paths import Paths


class DESeq2Controller(object):
    
    def __init__(self, args):
        self._args = args
        self._paths = Paths(args)
        self._helpers = Helpers(args)

    def vis_deseq2(self):
        """Generate plots based on the DESeq2 analysis"""
        conditions = list(set(self._args.conditions.split(',')))
        comparison_path_template_1 = "{}/DESeq2_comp_{}_vs_{}_with_annotation_and_countings.csv".format(
            self._paths.deseq2_extended_folder, conditions[0], conditions[1])
        comparison_1 = "{}_vs_{}".format(conditions[0], conditions[1])
        deseq2_vis = DESeq2Vis(
            comparison_path_template_1, self._paths.vis_deseq2_base_folder,
            self._args.padj_cutoff, comparison_1, self._args.alpha,
            self._args.color_sig, self._args.color_non_sig, self._args.shape,
            self._args.glyph_size)
        deseq2_vis.read_and_modificate_input()
        comparison_path_template_2 = "{}/DESeq2_comp_{}_vs_{}_with_annotation_and_countings.csv".format(
            self._paths.deseq2_extended_folder, conditions[1], conditions[0])
        comparison_2 = "{}_vs_{}".format(conditions[1], conditions[0])
        deseq2_vis = DESeq2Vis(
            comparison_path_template_2, self._paths.vis_deseq2_base_folder,
            self._args.padj_cutoff, comparison_2, self._args.alpha,
            self._args.color_sig, self._args.color_non_sig, self._args.shape,
            self._args.glyph_size)
        deseq2_vis.read_and_modificate_input()
        
