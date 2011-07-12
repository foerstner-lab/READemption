class Parameters(object):
    segemehl_accuracy = 85
    segemehl_hit_strategy = "1"
    segemehl_max_e_value = 5.0
    segemehl_number_of_threads = 1
    python_number_of_threads = 3

    # This define how expection should be treated exspecilly the once
    # that run in parellel. The option are "report" or "crash"
    # - report: The exceptions are written to stderr but the program
    #           should continue
    # - crash: the exception is raised
    exception_handling = "report"
    # Filtering
    min_seq_length = 12
    max_a_content = 70.0
    min_overlap = 1
    min_read_overlap_perc = 0.0
    uniquely_mapped_reads_only = True

