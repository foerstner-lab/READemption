import json

class RawStatDataWriter(object):

    def __init__(self, pretty=False):
        self._pretty = pretty

    def write(self, input_object, output_path):
        output_fh = open(output_path, "w")
        self._write_json(input_object, output_fh)
        output_fh.close()

    def _write_json(self, input_object, output_fh):
        if self._pretty is True:
            indent = 4
        else:
            indent = None
        output_fh.write(json.dumps(input_object, indent=indent))

class RawStatDataReader(object):

    def read(self, input_file):
        input_fh = open(input_file)
        data = self._read(input_fh)
        input_fh.close()
        return(data)

    def _read(self, input_fh):
        return(json.loads(input_fh.read()))
