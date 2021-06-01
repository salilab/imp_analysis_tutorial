import IMP
import unittest
import os
import shutil
import sys
import subprocess

test_dir = os.path.dirname(os.path.abspath(__file__))

class Tests(unittest.TestCase):
    
    def test_python_modeling_script_serial(self):
        os.chdir(os.path.join(test_dir, "..", "..", "rnapolii", "modeling/"))
        p = subprocess.check_call([sys.executable, "modeling.py", "test", "1", "10"])
        
        # Ensure that a stat file was created in the correct directory
        self.assertTrue(os.path.exists("test1/stat.0.out"))

        # clean up the directory
        shutil.rmtree("./test1")

        # Go back to the test home
        os.chdir(test_dir)

if __name__ == '__main__':
    unittest.main()
