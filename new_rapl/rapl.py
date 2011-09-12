#!/usr/bin/env python

"""RAPL - A RNA-seq Analysing PipeLine"""

import argparse
from libs.controller import Controller

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='commands')

    start_project_parser = subparsers.add_parser(
        'start', help='Start a project')
    start_project_parser.add_argument(
        'project_name', action='store', help='Name of the new project')
    start_project_parser.set_defaults(func=start_project)

    read_mapping_parser = subparsers.add_parser(
        'map', help='Run read mappings')
    read_mapping_parser.add_argument(
        '--force', '-f', default=False, action='store_true',
        help='Overwrite existing files.')
    read_mapping_parser.add_argument(
        '--threads', '-t', default=None, action='store',
        help='Number of threads that should be used.', type=int)
    read_mapping_parser.set_defaults(func=map_reads)

    gr_creation_parser = subparsers.add_parser(
        'gr', help='Create GR files')
    gr_creation_parser.add_argument(
        '--force', '-f', default=False, action='store_true',
        help='Overwrite existing files.')
    gr_creation_parser.set_defaults(func=create_gr_files)

    annotation_overlap_parser = subparsers.add_parser(
        'annotate', help='Search annoation overlaps')
    annotation_overlap_parser.add_argument(
        '--force', '-f', default=False, action='store_true', 
        help='Overwrite existing files.')
    annotation_overlap_parser.set_defaults(func=search_annoation_overlaps)

    generate_report_parser = subparsers.add_parser(
        'report', help='Generate a report')
    generate_report_parser.add_argument(
        '--force', '-f', default=False, action='store_true', 
        help='Overwrite existing files.')
    generate_report_parser.set_defaults(func=generate_report)
    
    args = parser.parse_args()
    controller = Controller()
    args.func(args, controller)

def start_project(args, controller):
    controller.start_project(args)

def map_reads(args, controller):
    controller.map_reads()

def create_gr_files(args, controller):
    controller.create_gr_files()

def search_annoation_overlaps(args, controller):
    controller.search_annotation_overlaps()

def generate_report(args, controller):
    controller.generate_report()

main()
