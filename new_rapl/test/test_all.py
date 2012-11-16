import unittest

all_tests = unittest.TestLoader().discover("./test")
unittest.TextTestRunner(verbosity=1).run(all_tests)
