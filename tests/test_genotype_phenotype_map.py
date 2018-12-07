#!/usr/bin/env python

__author__ = "Hoi Kin Cheng"
__copyright__ = "Copyright 2019-, The SMEAGL Project"
__credits__ = ["Jesse Zaneveld","Hoi Kin Cheng"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Hoi Kin Cheng"
__email__ = "hoikin@uw.edu"
__status__ = "Development"

from numpy import array
import unittest
from scripts.genotype_phenotype_map import GenotypePhenotypeMap

"""
Tests for simulate_microbiome.py
"""

class TestGenotypePhenotypeMap(unittest.TestCase):
    """Test of the microbiome simulation"""

    def setUp(self):
        """Set up data that will be used in multiple tests"""
        self.SmallGenome = array([0.0,0.0,0.5,1.0])
        self.Map = GenotypePhenotypeMap(genome_length=4)


    def test_GenotypePhenotypeMap_calculates_phenotypes_with_correct_data(self):
        """Additive phenotype map calculates phenotype from genotype"""
       
        min_val = 0.0
        max_val = 1000.0
        phenotype = "bigness"
        
        self.Map.add_phenotype(phenotype,min_val,max_val,\
          percent_of_sites=1.0)
       
        genome = self.Map.Phenotypes[phenotype]["max_trait_array"]
        #The genome should therefore have the max value for this trait
        
        result = self.Map.get_phenotypes(genome) 

        
        obs = result["bigness"]
        exp = 1000.0
        self.assertEqual(obs,exp)
    
    def test_get_phenotype_raises_array_size_exception(self):
        min_val = 0.0
        max_val = 1000.0
        phenotype = "bigness"
        
        self.Map.add_phenotype(phenotype,min_val,max_val,\
          percent_of_sites=1.0)
        
        # expected array size is 4 but passing in size 5 array to raise exception
        wrong_size_array = array([0.0,0.0,0.5,1.0,0.4])
        
        with self.assertRaises(ValueError):
            self.Map.get_phenotype("bigness", wrong_size_array)
            
    def test_get_phenotype_raises_name_not_found_exception(self):
        min_val = 0.0
        max_val = 1000.0
        phenotype = "bigness"
        
        self.Map.add_phenotype(phenotype,min_val,max_val,\
          percent_of_sites=1.0)
        
        
        with self.assertRaises(ValueError):
            # expected name is "bigness" but passing in "no name" to raise exception
            self.Map.get_phenotype("no name", self.SmallGenome)
            
    def test_min_value_equals_max_value(self):
        min_val = 1000.0
        max_val = 1000.0
        phenotype = "bigness"
        
        self.Map.add_phenotype(phenotype,min_val,max_val,\
          percent_of_sites=1.0)
        
        genome = self.Map.Phenotypes[phenotype]["max_trait_array"]
        #The genome should therefore have the max value for this trait
        
        result = self.Map.get_phenotypes(genome) 

        
        obs = result["bigness"]
        exp = 1000.0
        self.assertEqual(obs, exp)
        
    
    # test percent_of_sites working correctly
    def test_percent_of_sites(self):
        min_val = 0.0
        max_val = 1000.0
        phenotype = "bigness"
        
        self.Map.add_phenotype(phenotype,min_val,max_val,\
          random_seed = 100, percent_of_sites=0.5)
        
        result = self.Map.get_phenotypes(self.SmallGenome)
        
        obs = result["bigness"]
        exp = 783.20337361305371
        
        self.assertEqual(obs, exp)

if __name__ == "__main__":
    unittest.main()
