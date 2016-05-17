import unittest
import lsv

dpsi_file="/home/cradens/home_base/splicing/spurlock_ranzani_sullivan/dpsi_txt_results/Ra_T2_Ra_T17.deltapsi_quantify_deltapsi.txt"

# Here's our "unit tests".
class LSVTests(unittest.TestCase):

    # Tests whether user-friendly mode works: this is annoying to paste, so I only tested it once.
    # def testUserFriendly(self):
    #     imported_lsv = lsv.LSV(user_friendly=True)

    def testUserFriendlyANDpathProvided(self):
        expected_error = "User friendly mode is set to True, but a "
        expected_error+= "non-empty string argument was provided for file_path."
        with self.assertRaises(IOError) as context:
            lsv.LSV(user_friendly=True,file_path=dpsi_file)
            self.assertTrue(expected_error in context.exception)

    def testFileExists(self):
        expected_error = "Please only provide an extant file path."
        with self.assertRaises(IOError) as context:
            lsv.LSV(file_path="/home/cradens/Downloads/nothing_here.txt")
            self.assertTrue(expected_error in context.exception)

    def testExtantFileIsntVoila(self):
        expected_error = "That's not a voila text file..."
        with self.assertRaises(ValueError) as context:
            lsv.LSV(file_path="/home/cradens/Downloads/RadensCaleb_CV (1).pdf")
            self.assertTrue(expected_error in context.exception)

    def testImport(self):
        imported_lsv = lsv.LSV(file_path=dpsi_file)

    def testLenOverride(self):
        imported_lsv = lsv.LSV(file_path=dpsi_file)
        assert(len(imported_lsv)==1397)

    def test_lsv_id_check(self):
        imported_lsv = lsv.LSV(file_path=dpsi_file, match_ids = "all")
        self.assertTrue(imported_lsv.lsv_id_check("ENSG00000178035:49065047-49065575:source"))
        self.assertFalse(imported_lsv.lsv_id_check("this isn't an LSV ID"))
        assert(len(imported_lsv)==1397)

    def testMatchID_real_ID(self):
        print "cow"
        cow = lsv.LSV(file_path=dpsi_file, match_ids = ["ENSG00000178035:49065047-49065575:source"])
        self.assertTrue(len(cow) == 1)

#     def testMatchID_fake_ID(self):
#         expected_error = "Not all IDs found in file..."
#         with self.assertRaises(ValueError) as context:
#             test = lsv.LSV(file_path=dpsi_file, match_ids = ["this isn't an LSV ID"])
#             print context.exception
#             self.assertTrue(expected_error in context.exception)
# 
#     def testMatchIDsArg_isList(self):
#         expected_error = "If not a list, match_ids needs to be 'all'"
#         with self.assertRaises(ValueError) as context:
#             imported_lsv = lsv.LSV(file_path=dpsi_file, match_ids = "ENSG00000178035:49065047-49065575:source")
#             self.assertTrue(expected_error in context.exception)



    # def testSubset(self):
    #     imported_lsv = lsv.LSV(file_path=dpsi_file)
    #     subset = imported_lsv.subset(["ENSG00000178035:49065047-49065575:source"])
    #     assert len(subset) == 1

def main():
    unittest.main()

if __name__ == '__main__':
    main()