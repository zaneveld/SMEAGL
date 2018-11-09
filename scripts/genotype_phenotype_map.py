import numpy as np
import numpy.ma as ma

class GenotypePhenotypeMap():

    def __init__(self, genome_length):
        self.GenomeLength = genome_length
        self.Phenotypes = {}
    
    # temporarily set percent_of_sites to 1.0 for now
    def add_phenotype(self, name, min_value, max_value, percent_of_sites = 1.0):
        """
        Adds a phenotype to the dictionary of phenotypes.

        :param name: the name of the phenotype (e.g., temperature tolerance, edge stickiness, etc.)
        :param min_value: minimum value that phenotype can have (for example the lowest temperature tolerace could be 5 degrees Celsius)
        :param max_value: maximum value that phenotype can have (for example the lowest temperature tolerace could be 300 degrees Celsius)
        :param percent_of_sites: optional parameter for the % of sites in the genome that contributes to the
        phenotype. (0.0 = mask all; 1.0 = no masks at all)
        :return: none
        """
        
        # random values between 0 and 1.
        max_trait_array = np.random.rand(1, self.GenomeLength)
        
        # masks percent of sites of max_trait_array
        if (percent_of_sites != 0.0):
            n_zeros = round(percent_of_sites * self.GenomeLength)
            n_ones = self.GenomeLength - n_zeros
            ones_array = np.ones(n_ones)
            zeros_array = np.zeros(n_zeros)
            ones_and_zeros_array = np.append(ones_array, zeros_array)
            np.random.shuffle(ones_and_zeros_array)
            max_trait_array = ma.array(max_trait_array, mask = ones_and_zeros_array)
            
        phenotype = {"min_value": min_value,"max_value": max_value,\
        "max_trait_array": max_trait_array}
        print (max_trait_array)
        
        # put {name : phenotype} into self.Phenotypes dictionary
        self.Phenotypes[name] = phenotype
        
    def get_phenotype(self, name, genome):
        """
        Returns the scaled trait value (optimum) of one specified phenotype 
        (e.g., temperature tolerance, edge stickiness, etc.)

    
        :param name: the name of the phenotype (e.g., temperature tolerance, edge stickiness, etc.)
        :param genome: the genome that is being examined
        :return: scaled trait value (optimal value)
        """
        # This checks whether genome has the same length as max_trait_array
        if (genome.size != self.Phenotypes[name]["max_trait_array"].size):
            # Raise error here when size of genome array does not match size of max_trait_array
            raise ValueError("The size of genome must equal the size of max_trait_array, which is " +\
            str(self.Phenotypes[name]["max_trait_array"].shape))

        # These steps calculate the scaled_trait_value
        diff = abs(genome - self.Phenotypes[name]["max_trait_array"])
        mean_diff = np.mean(diff)
        slope = self.Phenotypes[name]["min_value"] - self.Phenotypes[name]["max_value"]
        scaled_trait_value = slope * mean_diff + self.Phenotypes[name]["max_value"]
        return scaled_trait_value
        
    # Is this one supposed to be static method? no.
    def get_phenotypes(self, genome):
        """
        Returns a dictionary of phenotypes (e.g., temperature tolerance, edge stickiness, etc.)
        and their corresponding scaled trait value.

    
        :param genome: the genome that is being examined
        :return: dictionary of phenotypes and their corresponding scaled trait value.
        """
        result = {}
        for phenotype in self.Phenotypes:
            result[phenotype] = self.get_phenotype(phenotype, genome)
        return result
        
def main():
    genome = np.array([1.0,0.4,0.25,0.13,0.0, 0.5, 0.7, 0.8, 0.3, 0.2])
    gpm = GenotypePhenotypeMap(10)
    gpm.add_phenotype("temp", 5, 300)
    gpm.add_phenotype("edge_stickiness", 0.0, 1.0)
    phenotypes = gpm.get_phenotypes(genome)
    print (phenotypes)
    
    genome2 = np.array([0.5,0.1,0.8,0.13,0.0, 0.5, 0.4, 0.8, 0.3, 0.2])
    phenotypes2 = gpm.get_phenotypes(genome2)
    print (phenotypes2)

if __name__ == "__main__":
    main()
