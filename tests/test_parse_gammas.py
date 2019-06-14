import unittest
import parse_gammas

class TestParseGammas(unittest.TestCase):
    
    def setUp(self):
        pass
        
    def tearDown(self):
        pass
        
    def test_strip_gamma(self):
        test_data = [
                        ("gamma(Y;Y,Z,Z)", "Y;Y,Z,Z")
                    ]
        for test_input, expected in test_data:
            test_output = parse_gammas.strip_gamma(test_input)
            self.assertTrue(test_output == expected, msg=str(test_output))

    
    def test_parse_gamma(self):
        test_data = [
                        ("@ gamma(Y;Y,Z,Z)       -190.6259", ((1, 1, 2, 2), -190.6259)),
                        ("@ gamma(Z;Z,Z,Z)       6205.8187", ((2, 2, 2, 2), 6205.8187)),
                        ("@ gamma(X;X,X,X)        -18.3970", ((0, 0, 0, 0), -18.3970))
                    ]
        for test_input, expected in test_data:
            self.assertTrue(parse_gammas.parse_gamma(test_input) == expected, msg="got {0}, expected {1}".format(test_input, expected))

    def test_convert_to_ints(self):
        test_data = [
                        ("Y;Y,Z,Z", (1, 1, 2, 2))
                    ]
        for test_input, expected in test_data:
            test_output = parse_gammas.convert(test_input)
            self.assertEqual(test_output, expected, msg="output: {0}, type: {1}, expected: {2}".format(test_output, type(test_output), expected))
