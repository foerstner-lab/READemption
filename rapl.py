#!/usr/bin/env python

"""RAPL - A RNA-seq Analysing PipeLine"""

import argparse
from rapl.rapl import Rapl

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='commands')

    start_project_parser = subparsers.add_parser(
        'startproject', help='Start a project')
    start_project_parser.add_argument(
        'project_name', action='store', help='Name of the new project')
    start_project_parser.set_defaults(func=start_project)

    read_mapping_parser = subparsers.add_parser(
        'mapreads', help='Run a read mapping')
    read_mapping_parser.add_argument(
        '--force', '-f', default=False, action='store_true',
        help='Overwrite existing files.')
    read_mapping_parser.add_argument(
        '--threads', '-t', default=None, action='store',
        help='Number of threads that should be used.', type=int)
    read_mapping_parser.set_defaults(func=map_reads)

    annotation_overlap_parser = subparsers.add_parser(
        'annotationsearch', help='Search annoation overlaps')
    annotation_overlap_parser.add_argument(
        '--force', '-f', default=False, action='store_true', 
        help='Overwrite existing files.')
    annotation_overlap_parser.set_defaults(func=search_annoation_overlaps)
    
    args = parser.parse_args()
    rapl = Rapl()
    args.func(args, rapl)

def start_project(args, rapl):
    print("Start it!")
    rapl.start_project(args)

def map_reads(args, rapl):
    print("Map it!")
    rapl.map_reads(args)

def search_annoation_overlaps(args):
    print("Search it!")

main()
