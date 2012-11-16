import argparse

__description__=""

def main():
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('input_fasta_file', help='A fasta file.')
    parser.add_argument('--output_file', '-o', dest='output_file',
                        help='An output file. If not give the output ' +
                        'is written to stdout.')
    args = parser.parse_args()
    clipper = Clipper(args.input_fasta_file, output_file=args.output_file)

class Clipper(object):

    def __init__(self, input_fasta_file, output_file=None):
        self.input_fasta_file = input_fasta_file
        self.output_file = output_file

    def clip(self):
        pass

    

if __name__ == '__main__': 
    main()



