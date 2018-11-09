from numpy import array,full,zeroes
"""
Simulate a niche based microbial community dynamics,
and microbial genome evolution
"""

class Environment(object):
    def __init___(self,x_size,y_size,resources={"Glucose":0.01},
        temperature={"mean":32.0,"type" = "uniform"},pressure={"mean":1.0,"type":"uniform"}):
        """Initialize an environment
        x_size -- environment dimensions
        y_size -- environment dimensions
        resource -- dictionary of resources
          keys are resource names,
           values are resource
          concentrations (mols/L)
        temperature -- integer
        temperature_type -- string options are 'uniform' (more will be added later)
        """
        self.Dimensions = (x_size,y_size)
        self.Resources = {}
        if temperature_type = "uniform":
            self.Temperatures = full(dimensions,temperature)

        #Initialize maps of resource concentrations
        #setting each to the starting concentration
        for resource,concentration in resources.items():

            resource_array =\
              full(dimensions,concentration,dtype=float)
            self.Resources[resource] = resource_array

class Reaction(object)
    def __init__(self,name,EC= None,reactants,products,delta_G):
        """
        name -- a string (example: "Alcoholic Fermentation")
        EC -- An enzyme commission id
        Reactants -- a dict of reactants and their stoichiometry
        example: {"C6H1206":1}
        products -- a dict of products and their stoichiometry
        """

        self.Name = name
        self.EC = EC
        self.Reactants = reactants
        self.Products = products
        self.ChemicalPotential = delta_G

    def calculte_DeltaG(self,environment,location_array):
        """Return an array of real (non-standard conditon) delta Gs across the Environment
            environment -- a reference to an Environment object


            deltaG = deltaGstandardState + gas_constant_R * tempearature * ln([productA]**productAstoichiometry*[productB]**productBstoichiometry)/[reactantA]**reactantAstoichiometry*[reactantB]**reactantBstoichiometry)
        Strategy:
        We will calculate an array of deltaG across the entire landscape array
        represented by the environment
        """

        gas_constant_R = 8.3144598
        temperatures_K = environment.Temperatures + 273.15
          #NOTE we add 273.15 because this must be in Kelvin, not Celsius!!!
        delta_G_SS = self.ChemicalPotential

        stoichiometry_numerator =\
         [environment[p]**self.Products[p] for p in self.Products.keys()]
        stoichiometry_denominator =\
         [environment[r]**self.Reactants[r] for r in self.Reactants.keys()]
        deltaG = delta_G_SS + gas_constant_R * temperatures * ln(stoichiometry_numerator/stoichiometry_demoninator)

class Strain(object):
    def __init__(self,environment,metabolic_capabilities={"Glucose":0.0001}):
        self.CellCount = zeroes(environment.Dimensions,dtype=float)
        self.Environment = environment
        self.

def main():
    """ """
    #test out the environment class
    #let's try modeling Rainbow vent
    #https://en.wikipedia.org/wiki/Rainbow_Vent_Field
    #Table 1 in German 2010 has flux values for compounds
    #https://www.sciencedirect.com/science/article/pii/S0967063710000105

    name = "Rainbow Hydrothermal Field"
    temperature = 365 #°C
    env = Environment(100,100,temperature=365.0)
    pH = 2.8
    fluxes_mols_second =\
     {"3He":7.6×10−9,
      "CH4":1.0,
      "Mn":0.9,
      "Fe":9.6,
      "Cu":0.11,
      "Zn":0.05,
      "Cd":5.2*10**−4,
      "Ag":1.1*10**−4,
      "Pb":4.5×10−4,
      "P":0.94,
      "V":0.036,
      "Co":2.5×10−3,
      "U":4.9×10−6,
      "Y":2.6×10−4,
      "La":1.2×10−4,
      "Ce":4.9×10−5,
      "Pr":2.4×10−5,
      "Nd":9.6×10−5,
      "Sm":1.8×10−5,
      "Eu":1.1×10−5,
      "Gd":1.9×10−5,
      "Tb":3.3×10−6,
      "Dy":2.0×10−5,
      "Ho":	4.3×10−6,
      "Er":	1.2×10−5,
      "Tm":	1.5×10−6,
      "Yb":	9.5×10−6,
      "Lu"	1.3×10−6}


if __name__ == "__main__":
    main()
