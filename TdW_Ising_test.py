import unittest
from numpy import array, testing
from TdW_Ising import Spin

class SpinTest(unittest.TestCase):

	def setUp(self):
		self.spin = Spin(array([1,1,-1,-1]),0)

	def test_flip_first(self):
		flipped_spins = self.spin.flip(0)
		testing.assert_array_equal([-1,1,-1,-1], flipped_spins)

	def test_flip_last(self):
		flipped_spins = self.spin.flip(3)
		testing.assert_array_equal([1,1,-1,1], flipped_spins)

	def test_flip_middle(self):
		flipped_spins = self.spin.flip(1)
		testing.assert_array_equal([1,-1,-1,-1], flipped_spins)

	def test_flip_single(self):
		spin = Spin(array([1]),0)
		flipped_spins = spin.flip(0)
		testing.assert_array_equal([-1], flipped_spins)

	def test_flip_empty(self):
		spin = Spin(array([]),0)
		flipped_spins = self.spin.flip(0)
		testing.assert_array_equal([], flipped_spins)

	def test_flip_outofbounds(self):
		self.assertRaises(IndexError, self.spin.flip, 5)

