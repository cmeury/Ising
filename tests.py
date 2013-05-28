import unittest
import TdW_Ising_test
import sys

if __name__ == "__main__":
	suite = unittest.TestLoader().loadTestsFromModule(TdW_Ising_test)
	unittest.TextTestRunner().run(suite)
