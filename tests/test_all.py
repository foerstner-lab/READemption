import unittest

all_tests = unittest.TestLoader().discover("./tests")
unittest.TextTestRunner(verbosity=1).run(all_tests)
