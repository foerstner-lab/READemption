#!/usr/bin/env python

"""RAPL - A RNA-seq Analysing PipeLine"""

import argparse
from libs.controller import Controller

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="commands")

    # Arguments for project creation
    start_project_parser = subparsers.add_parser(
        "start", help="Start a project")
    start_project_parser.add_argument(
        "project_path", default=".", help="Name/path of the project.")
    start_project_parser.set_defaults(func=start_project)

    # Parameters for read mapping
    read_mapping_parser = subparsers.add_parser(
        "map", help="Run read mappings")
    read_mapping_parser.add_argument(
        "project_path", default=".", nargs="?", 
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    read_mapping_parser.add_argument(
        "--min_read_length", "-l", default=12, type=int,
        help="Minimal read length after clipping")
    read_mapping_parser.add_argument(
        "--threads", "-t", default=1, type=int,
        help="Number of threads that should be used.")
    read_mapping_parser.add_argument(
        "--force", "-f", default=False, action="store_true",
        help="Overwrite existing files.")
    read_mapping_parser.add_argument(
        "--segemehl_accuracy", "-a", default=95.0, type=float,
        help="Segemehl's minimal accuracy (in %%) (default 95).")
    read_mapping_parser.add_argument(
        "--segemehl_evalue", "-e", default=5.0, type=float,
        help="Segemehl's maximal e-value. (default 5.0)")
    read_mapping_parser.add_argument(
        "--segemehl_bin", "-s", default="segemehl",
        help="segemehl's binary path.")
    read_mapping_parser.set_defaults(func=map_reads)
    read_mapping_parser.add_argument(
        "--samtools_bin", "-b", default="samtools",
        help="samtools' binary path.")

    # Parameters for coverage file building
    coverage_creation_parser = subparsers.add_parser(
        "coverage", help="Create coverage (WIGGLE) files")
    coverage_creation_parser.add_argument(
        "project_path", default=".", nargs="?", 
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    coverage_creation_parser.add_argument(
        "--unique_only", "-u", default=False, action="store_true",
        help="Use uniquely mapped reads only.")
    coverage_creation_parser.add_argument(
        "--force", "-f", default=False, action="store_true",
        help="Overwrite existing files.")
    coverage_creation_parser.set_defaults(func=create_coverage_files)
    coverage_creation_parser.add_argument(
        "--threads", "-t", default=1, type=int,
        help="Number of threads that should be used.")
    coverage_creation_parser.add_argument(
        "--samtools_bin", "-b", default="samtools",
        help="samtools' binary path.")

    # Parameters for annotation overlap searches
    annotation_overlap_parser = subparsers.add_parser(
        "annotate", help="Search annoation overlaps")
    annotation_overlap_parser.add_argument(
        "project_path", default=".", nargs="?", 
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    annotation_overlap_parser.add_argument(
        "--force", "-f", default=False, action="store_true", 
        help="Overwrite existing files.")
    annotation_overlap_parser.set_defaults(func=search_annoation_overlaps)

    # Parameters for annotation overlap overviews
    annotation_overlap_parser = subparsers.add_parser(
        "annotation_overview", help="Create annotation overlap overviews")
    annotation_overlap_parser.add_argument(
        "project_path", default=".", nargs="?", 
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    annotation_overlap_parser.add_argument(
        "--unique_only", "-u", default=False, action="store_true",
        help="Use uniquely mapped reads only.")
    annotation_overlap_parser.add_argument(
        "--min_overlap", "-o", default=10, type=int,
        help="Minimal read-annotation-overlap (in nt) (default 10)")
    annotation_overlap_parser.add_argument(
        "--force", "-f", default=False, action="store_true", 
        help="Overwrite existing files.")
    annotation_overlap_parser.set_defaults(func=search_annoation_overlaps)

    # Parameters for report generation
    generate_report_parser = subparsers.add_parser(
        "report", help="Generate a report")
    generate_report_parser.add_argument(
        "project_path", default=".", nargs="?", 
        help="Path of the project folder. If none is given the current "
        "directory is used.")
    generate_report_parser.add_argument(
        "--force", "-f", default=False, action="store_true", 
        help="Overwrite existing files.")
    generate_report_parser.set_defaults(func=generate_report)
    
    args = parser.parse_args()
    controller = Controller(args)
    args.func(controller)

def start_project(controller):
    controller.start_project()

def map_reads(controller):
    controller.map_reads()

def create_coverage_files(controller):
    controller.create_coverage_files()

# OBSOLTE
# def create_gr_files(controller):
#     controller.create_gr_files()

def search_annoation_overlaps(controller):
    controller.search_annotation_overlaps()

def generate_report(controller):
    controller.generate_report()

main()
