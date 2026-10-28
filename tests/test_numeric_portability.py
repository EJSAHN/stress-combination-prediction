"""Numerical serialization checks; real changes must still fail validation."""
import sys
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from stress_prediction.replay import compare_json_vector


class NumericVectorTests(unittest.TestCase):
    def test_last_digit_serialization_difference_passes(self):
        count,error=compare_json_vector('[116.42883699734087,69.12990093033048]',
                                        '[116.42883699734104,69.12990093033045]')
        self.assertEqual(count,2)
        self.assertLess(error,1e-10)

    def test_changed_values_are_rejected(self):
        with self.assertRaises(AssertionError):
            compare_json_vector('[1,2,3]','[1,2,3.001]')

    def test_invalid_numeric_arrays_are_rejected(self):
        for other in ['[1,2]', '[1,2,"3"]', '[1,2,NaN]', '[true,2,3]', '[]', 'not JSON']:
            with self.subTest(other=other),self.assertRaises((ValueError,AssertionError)):
                compare_json_vector('[1,2,3]',other)


if __name__=='__main__': unittest.main()
