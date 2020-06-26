class ParameterLog(object):
    def __init__(self, log_file):
        self.log_file = log_file
        self.delimiter = ": "


class ParameterLogger(ParameterLog):
    def add_paramter(self, parameter, value):
        log_fh = open(self.log_file, "a")
        log_fh.write(self.delimiter.join([parameter, value]) + "\n")
        log_fh.close()


class ParameterLogReader(ParameterLog):
    def read_parameters(self):
        parameters = {}
        for line in open(self.log_file):
            paramter, value = line[:-1].split(self.delimiter)
            parameters[parameters] = value
