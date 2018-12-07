"""
A simulation of microbial ecology in a hydrothermal vent
"""
import numpy as np
from matplotlib.pyplot import imshow,show,figure,colorbar,get_cmap
from matplotlib.image import imread
from matplotlib.colors import ListedColormap

import matplotlib.animation as animation
from optparse import OptionParser
from optparse import OptionGroup
from random import randint,random
from copy import copy


def additive_phenotype_map(genome,phenotype='temperature_optimum',site_mask=None,\
  high_trait_value_map=None,min_phenotype=-20,max_phenotype=100.0):
    """Determine a property linearly using the sum of alleles at certain loci
    
    genome - an array with float values between 0 and 1
    site_mask - an array of 0 or 1 values that determines which sites in the genome
    will affect the property (e.g. temperature optimum)

    high_trait_values_map - an array of values (typically random) that would correspond
    to a maximum trait value. The idea of this property is that if each phenotype has
    a different random high_trait_values_map, then intrinsically there will be tradeoffs
    between traits. site_masks (which assert that only certain 'genes' affect certain
    traits) can be used to simulate the case where there are no tradeoffs between 
    phenotypes. 


    min_phenotype - minimal value of the phenotype
    max_phenotype - maximum value of the phenotype
    """

    if site_mask is None:
        site_mask = np.ones(genome.shape)
        #Ignore the allele at irrelevant loci (they are 0s in site mask)
        masked_genome = genome * site_mask 
    
    obs_value = np.sum(masked_genome)
    min_value = 0.0
    max_value = np.sum(site_mask) #If all sites are 1 the total will = sum of site mask
    trait_range = max_value - min_value

    relative_trait_value = (obs_value - min_value)/(trait_range)

    obs_phenotype_value = float(min_value + relative_trait_value*(trait_range))
    
    return {phenotype:obs_phenotype_value}

class Microbe():
    """A strain of microbe
    
    The class handles a microbial strain in the simulation
    It keeps track of its own distribution based on the map
    """

    def __init__(self,distribution=None,traits={},name=None,genome=None,max_genome_size=1000,
        spread_rates = {"up":0.025,"down":0.15,"right":0.05,"left":0.05},\
        migration_zones=[((0.0,1.0),(0.0,1.0))],migration_rate=10**2,\
        migration_probability=0.1,migration_sites=5):
        """Initialize a Microbe

        distribution - a numpy array showing the count of the microbe on the map
        traits - a dict of valued traits (e.g. optimal temperature)
        name - the microbes name. Will be randomly generated if None
        genome - a 1d numpy array, scaled between 0 and 1, representing the genome
        max_genome_size - if no genome is provided, use max_genome_size to initialize
         a random genome
        migration_zones - a list of sites describing proportionally where in the map to place migrants
        migration_rate - number of migrants added per event
        migration_probability - chances of migration at each possible migration site
        migration_sites - number of random sites per timestep where migration may occur
        """

        self.Name = name
        #Rates of migration are currently hardcoded but could be traits that evolve
        self.SpreadRates = spread_rates
        if self.Name == None:
            self.Name = "microbe_{}".format(randint(10000000))
        self.Distribution = distribution
        self.Genome = np.zeros((1,max_genome_size),dtype=float)
        self.Traits = self.get_traits_from_genome()

        #Each microbe can be added to the map in different areas at different rates
        self.MigrationZones = migration_zones #where to add migrants
        self.MigrationRate = migration_rate #how many to add
        self.MigrationProbability = migration_probability #chance to add migrants per timestep 
        self.MigrationSites = migration_sites 
        #Traits to make:
        #intrinsic growth rate
        #temperature optimum

        #Properties that probably don't evolve
        #max_temp_mortality
        #intinsic mortality

        #If the user specified some other traits, update them
        self.Traits.update(traits)
    
    def get_traits_from_genome(self,phenotype_map_fn=additive_phenotype_map):
        """Get a new dictionary of traits by examining the genome
        
        This function is intented to serve as a bridge between the genome
        and organism traits. For different applications this can be very simple (e.g.
        the 'genome' is a 1-element array and the value translates into a temperature
        optimum), or arbitrarily complex (the genome is a huge array and metabolic
        modeling is used to determine seed compounds).

        phenotype_map_fn -- a function that generates a dict of traits given a genome
        The function must return a dict of properties (e.g. temperature tolerance)
        that are used by a Simulation object.
        """
        
        updated_traits = phenotype_map_fn(self.Genome)
        return updated_traits
    
    def simulate_migration(self):
        """Add migrant microbes to self.Distribution
        
        This function depens on the following properties of the Microbe:

        self.MigrationZones -- portions of the map in which new migrants can
          appear (they will randomly appear in this zone)
        
        self.MigrationSites -- how many squares should 
          migrants be placed in?
        
        self.MigrationProbability -- per-site chances that migration occurs
        
        self.MigrationNumber -- number of migrant microbes per square
        """
        for zone in self.MigrationZones: 
            (min_y,max_y),(min_x,max_x) = zone
            height,width = self.Distribution.shape
            for i in range(int(self.MigrationSites)):
                if random() < self.MigrationProbability:
                    new_microbe_x = randint(int(min_x*width-1),int(max_x*width-1))
                    new_microbe_y = randint(int(min_y*height-1),int(max_y*height-1))
                    self.Distribution[new_microbe_y,new_microbe_x] += self.MigrationRate
    
    def simulate_spread(self,geometry=None,edge_map=None,edge_bonus=0.10):
        """Spread microbes to neighboring cells

        geometry -- a map of same dimensions as self.Distribution. 0 values indicate
          places a microbe can't spread

        edge_map -- a map of edges. 1 values mean the edge bonus is applied to
         spread rate (reflecting faster micorbial movement along edges)

        edge_bonus -- if >0.0, allow microbes on edges to spread faster
          (in any direction)
        """
        #Check that the spread does not increase cell numbers
        for direction in ["up","down","left","right"]:
            spread_rate = self.SpreadRates[direction]
            offset_array = generate_offset_array(copy(self.Distribution),direction=direction)
            
            #Microbes can't spread into solid rock
            #Optionally apply geometry to the offset array
            #any value times 0 = 0, so black (0 value) areas should not be included
            if geometry is not None:
                offset_array = offset_array * geometry
            
            #Optionally adjust the spread rate to be higher/lower on edges
            if edge_map is not None:
                adj_spread_rate = spread_rate + edge_bonus * edge_map
            else:
                adj_spread_rate = spread_rate
            
            #Add microbes to new cells and subtract them from the old cell
            new_data = self.Distribution + offset_array * adj_spread_rate  
            new_data = new_data - self.Distribution * adj_spread_rate#conservation of microbes

            self.Distribution = new_data



     
    def simulate_growth(self,growth_rate):
        """All micorbes divide"""
        #For now assume 'good' environments have both higher
        #growth rate and carrying capacity.

        random_array = self.get_uniform_random_array()
        growth_map = growth_rate
        growth_map = growth_map * random_array + 0.5 #add 50% stochasticity to growth rates
        growth_map = np.maximum(1.0,growth_map) #no mortalit in this step
        growth_map = np.minimum(2.0,growth_map) #no growth faster than 2.0
        self.Distribution = self.Distribution*growth_map
   
    def simulate_mortality(self,temperature_map,intrinsic_mortality=0.001,max_temp_mortality=0.005):
        """Simulate mortality due to temperature
        Microbes are assumed to have an optimum temperature and an extent of temperature
        sensitivity which reperesents the slope of increased mortality as they move
        away from their temperature optimum
        """

        temperature_optimum  = self.Traits["temperature_optimum"]
        temperature_sensitivity  = max_temp_mortality

        #The difference in temperature from the microbe's thermal optimum
        #is scaled by it's temperature sensistivity
        survival_map = 1.0 - ( np.absolute(temperature_map - temperature_optimum) * temperature_sensitivity)
        
        survival_map -= intrinsic_mortality
        survival_map[survival_map<0]=0 #Can't have negative survival
        self.Distribution = survival_map*self.Distribution
 
 
    def get_uniform_random_array(self):
        """Return a uniform random array between 0 and 1 """
        random_array = np.random.random(self.Distribution.shape)
        return random_array

      
class Simulation():
    def __init__(self,geometry=None,width=500,height=500,base_carrying_capacity=10**12,\
      edge_carrying_capacity=0,upwards_flow_map = None,edge_map=None,temperature_map = None,\
      microbes={},starting_microbes=2):
       
        #Set the simulation to timestep 0 
        self.TimeStep = 0
        self.Width = width
        self.Height = height
        
        self.Geometry = self.load_map_data(geometry)
        #If a geometry map was specified, override
        # user-specied height and width
        height,width = geometry.shape
        self.Width = width
        self.Height = height
        
        #Load maps        
        self.Temperatures = self.load_map_data(temperature_map)
        self.UpwardsFlowMap = self.load_map_data(upwards_flow_map)
        self.EdgeMap = self.load_map_data(edge_map)


        #Set up the starting population of microbes
        self.Microbes = microbes
        #If we don't have enough prespecified microbes
        #randomly generate more
        if len(microbes) < starting_microbes:
            n_microbes_to_generate = starting_microbes - len(microbes)
            self.make_random_microbial_strains(n_microbes_to_generate)
 
        #Calcuate maps that are a combination of input data and input parameters
        self.CarryingCapacityMap = base_carrying_capacity + (self.EdgeMap * edge_carrying_capacity) 
        
        #Define parameters that are currently hardcoded
        self.MigrationRate = 10**4
    
    def load_map_data(self,image,dtype=float):
        """Load map data from an image
        if the image is None, generate an empty map
        """

        if image is None:
            return np.zeroes(self.Height,self.Width,dtype=dtype)
        else:
            return image 
    
    def make_random_microbial_strains(self,n_microbes):
        """Make new random unrelated microbes, but keep their distribution empty"""
        for i in range(n_microbes):
            
            starting_distribution =  np.zeros((self.Height,self.Width),dtype=float)
            #Each 'microbe' is a specific strain
            #It has a name, genome and distribution (where it is on the map)

            new_microbe =\
              Microbe(name="random_microbe_%i"%i,\
              distribution=starting_distribution)

            self.Microbes[new_microbe.Name]=new_microbe
 
    def update(self,verbose=True,super_verbose=False):
        """Run one round of the simulation"""
        
        if self.TimeStep % 10 == 0 and verbose:
            print("Simulating timestep {}".format(self.TimeStep))
        
        for microbe_name,microbe in self.Microbes.items():
            #Stuff the environment 'does to' the microbe
            

            self.simulate_vent_flow(target=microbe)

            if super_verbose:
                print("Total population ({}) - {} phase:{}".format(microbe_name,"Vent Flow",\
                  np.sum(microbe.Distribution)))    
            
            microbe.simulate_migration()
            
            if super_verbose:
                print("Total population ({}) - {} phase:{}".format(microbe_name,"Migration",\
                  np.sum(microbe.Distribution)))    
            
            #Microbe objects know how to grow,die,etc    
            microbe.simulate_growth(growth_rate=2.0)
            if super_verbose:
                print("Total population ({}) - {} phase:{}".format(microbe_name,"Growth",\
                  np.sum(microbe.Distribution)))    
            microbe.simulate_mortality(temperature_map=self.Temperatures) #from temperature for now
           
            if super_verbose:
                print("Total population ({}) - {} phase:{}".format(microbe_name,"Mortality",\
                  np.sum(microbe.Distribution)))   
 
            microbe.simulate_spread(geometry=self.Geometry,edge_map=self.EdgeMap) #along surfaces
        
            if super_verbose:
                print("Total population ({}) - {} phase:{}".format(microbe_name,"Spread",\
                  np.sum(microbe.Distribution)))  
 
            if super_verbose:
                print("Total population ({}):{}".format(microbe_name,np.sum(microbe.Distribution)))    
        #Carrying capacity adds up counts of *all* microbes,
        #so it is outside the for loop    
        self.apply_carrying_capacity() 
        
        self.TimeStep += 1


    def total_microbes(self):
        """Return an array of the total number of microbes per square"""
        result = None
        for microbe_name,microbe in self.Microbes.items():
            if result is None:
                result = microbe.Distribution
            else:
                result += microbe.Distribution

        return result

    def apply_carrying_capacity(self):
        """Remove microbes above the per-square carrying capacity
        
        Steps:
        1. Figure out total microbes per square and 
           carrying capacity (max micorbes per square)

        2. Subtract them to figure out how many microbes need to die.
        3. Remove a proportional amount from each strain
            (NOTE: could also use weighted random number generation
            here to make this stochastic)
        """
        total_microbes_per_square = self.total_microbes()

        number_to_kill_per_square = total_microbes_per_square - self.CarryingCapacityMap
        #The max is there to make sure we don't kill negative numbers of microbes
        number_to_kill_per_square[number_to_kill_per_square < 0] = 0

        epsilon = 1/10**10
        for microbe_name,microbe in self.Microbes.items():
            proportion_per_square = microbe.Distribution/(total_microbes_per_square+epsilon) 
            
            #If a microbe is 70% of the population in a square, and we are 100 microbes
            #above the carrying capacity, 70 microbes from that strain should die.
            new_distribution = microbe.Distribution -\
             (number_to_kill_per_square * proportion_per_square) 

            #Ensure no cell counts are negative
            new_distribution[new_distribution < 0] = 0

            new_distribution *= self.Geometry 
            #self.Geometry should be 0 or 1, so any microbes in solid rock are removed
            
            #Round the result so we don't generate fractional cell counts
            microbe.Distribution = np.round(new_distribution,0)

    def simulate_vent_flow(self,target,vent_velocity=1.00,edge_stickiness=0.70):
        """Push target upwards due to the vent flow

        target - an object with a Distribution property (usually a microbe)

        vent_velocity - the percentage of counts in each cell that are pushed 1 square
         upwards per timestep (if not on an edge)

        edge_stickiness - if target is on an edge, subtract edge_stickiness % of 
          counts from being pushed upwards (e.g. higher numbers = microbes stick to edges)

        """
        offset_array = generate_offset_array(copy(target.Distribution),direction="up")
        #Flow occurs at different rates, and microbes on surfaces 'stick'
        #we subtract the edgemap to simulate this
        random_array = self.get_uniform_random_array()+0.50
        spread_rate = (self.UpwardsFlowMap - self.EdgeMap*edge_stickiness) * vent_velocity*random_array 
        offset_array = offset_array * self.Geometry * self.UpwardsFlowMap
        #Microbes can't spread into solid rock
        #apply geometry to the offset array (handled by multiplying by the geometry map)
        #any value times 0 = 0, so black areas should not be included

        new_data = target.Distribution + offset_array * spread_rate
        new_data = new_data - target.Distribution * spread_rate #conservation of microbes
        target.Distribution = new_data
    
    def get_uniform_random_array(self):
        """Return a uniform random array between 0 and 1 """
        random_array = np.random.random((self.Height,self.Width))
        return random_array

    def get_data(self,focal_microbe,scaling ="log10"):
        """Return simulation data as a numpy array
        focal_microbe -- microbe for which we are getting data
        scaling -- if set, rescale the data for visualization
        """
     
        if scaling == 'linear':
            return self.Microbes[focal_microbe].Distribution
        elif scaling == 'log2':
            eps = 0.1/self.CarryingCapacityMap
            #add a tiny epsilon value to avoid divide by zero errors
            return np.log2(self.Microbes[focal_microbe].Distribution + eps)
        elif scaling == 'log10':
            eps = 0.1/self.CarryingCapacityMap
            #add a tiny epsilon value to avoid divide by zero errors
            return np.log10(self.Microbes[focal_microbe].Distribution + eps)


def generate_offset_array(input_array,direction="up"):
    """generate an offset array
    input_array - an input 2d array
    direction - string representing one of four directions: up, down,left,and right
    Strategy: we want to get all the cells adjactent to every 
    cell on the map - all at once. One way to do this is to shift the
    whole board down, then superimpose that times some constant
    on the old board. That has the effect of adding something to one
    cell above it's old position.
    """
 
    if direction == "down":
        offset_array = np.pad(input_array,[(1,0),(0,0)],mode='constant',constant_values=0)
        offset_array = offset_array[:-1,:]
    
    elif direction == "up":
        offset_array = np.pad(input_array,[(0,1),(0,0)],mode='constant',constant_values=0)
        offset_array = offset_array[1:,:]

    elif direction == "right":
        offset_array = np.pad(input_array,[(0,0),(1,0)],mode='constant',constant_values=0)
        offset_array = offset_array[:,:-1]
        
    elif direction == "left":
        offset_array = np.pad(input_array,[(0,0),(0,1)],mode='constant',constant_values=0)
        offset_array = offset_array[:,1:]

    return offset_array

def make_plot_updater(plots,simulation,width=500,height=500):
    """Make an updatefig function for the specified plot
    plot -- a dict of matplotlib plots, keyed by microbe
    simulation -- a class object with an get_data method
      that returns a width * height numpy array of numbers
    width -- width in pixels
    height -- height in pixels

    """
    def update_fig(frame,*args):
        simulation.update()
        for microbe,plot in plots.items():
            new_data = simulation.get_data(focal_microbe=microbe)
            #Use plot.set_data to update the data in the plot
            #without making a new plot
            plot.set_data(new_data)
            #The comma is needed to make this a tuple and
            #therefore iterable
        return plots.values()

    return update_fig

def load_map_image(path,expected_height=None,expected_width=None,convert_to_grayscale=True):
    """Load an image and return a numpy array

    path -- path to a .png image file

    expected_height -- the expected height of the image in pixels. If provided,
      raise an error if the image is not this height.

    expected_width -- the expected width of the image in pixels. If provided,
      raise an error if the image is not this width.

    convert_to_grayscale -- if True, convert the image to a single number per cell
    using the blue channel. For true black / true white cells (which is all we
    expect), it doesn't matter which RGB channel we use. 
    """
    result = imread(path) #load the image an array
    if np.ndim(result)==3: #if three dimensional, this is an RGB image
        if convert_to_grayscale == True:
            result = result[:,:,0] 
            #Now we have a 2d (grayscale) array instead of a 3d array
            #with (R,G,B) values for each cell
    
    height,width = result.shape
    if not expected_height:
        return result
    
    #If we have expected dimensions for the map, check them
    #(and maybe raise an error) before returning the array
    if  height != expected_height or width !=expected_width:
        err_msg = "Image {} had dimesions {} x {} not {} x {}".format(\
         path,expected_width,expected_height,width,height)
        raise ValueError(err_msg)
    return result


def initiate_simulation_plot(geometry,visual_overlay,\
  microbe_color_maps,geometry_colormap='viridis',min_heatmap_value=0.0,max_heatmap_value=10**12):
    """Return a figure and heatmap subplot for animation
    
    geometry - a numpy array reflecting the geometry of the map
    
    visual_overlay - a visual overlay that is plotted on top of the data
    microbe_color_maps - a dict of color maps to use for each microbe
      each color map will say what color set to use to convert raw numbers (e.g. how many microbes)
      into a color to plot. Examples include 'jet','magma' and 'viridis'. 
      For valid values, see https://matplotlib.org/examples/color/colormaps_reference.html

    min_heatmap_value - what number should map to the first color on the color scheme?
    max_heatmap_value - what number should map to the last color on the color scheme?
    """
    #Initiate the base plot
    fig = figure()
    vent_plots = {}
    fig.patch.set_facecolor('royalblue')
    for microbe_name,color_scheme in microbe_color_maps.items():
        vent_plot = imshow(geometry,cmap=color_scheme,\
          vmin=min_heatmap_value,vmax=max_heatmap_value,animated = True)
        vent_plots[microbe_name]=vent_plot
        curr_colorbar =colorbar(shrink=0.25)
        curr_colorbar.set_label(microbe_name)
        
    imshow(visual_overlay)
    
    #fig is the overall figure
    #vent_plot is the specific heatmap part (we need both)
    return fig,vent_plots
 
def make_transparent_cmap(cmap):
    """make a transparent version of a matplotlib colormap
    
    From the stackoverflow answer here: \
      https://stackoverflow.com/questions/37327308/add-alpha-to-an-existing-matplotlib-colormap
    """
    # Get the colormap colors
    my_cmap = cmap(np.arange(cmap.N))

    # Set alpha
    my_cmap[:,-1] = np.linspace(0.3, 0.5, cmap.N)

    # Create new colormap
    my_cmap = ListedColormap(my_cmap)
    return my_cmap

    
	
def make_option_parser():
    """Return an optparse OptionParser object"""
    parser = OptionParser(usage = "SMEAGL (Simulated Microbial Evolution and " + 
    "Genomic LGT) is designed to simulate microbial ecology and evolution. " + 
    "The initial application is to explore the hypothesis that high " + 
    "abundances of DNA mobility elements in some deep sea hydrothermal " + 
    "vent microbiomes may be more strongly selected in low-richness communities")

    optional_options = OptionGroup(parser, "Optional options")

    optional_options.add_option('-b', '--base_carrying_capacity',\
    default=10**4, type="int",\
    help="How many microbes fit in a square (not an edge)? [default:%default]")
    
    optional_options.add_option('-e', '--edge_carrying_capacity',\
    default=10**9, type="int",\
    help="How many microbes fit in a square if the square is an edge? [default:%default]")
    
    optional_options.add_option('-l', '--simulation_length',\
    default=1200, type="int",\
    help="How many timesteps shall we run the simulation? [default:%default]")
    
    optional_options.add_option('-i', '--input_image',\
    default="../data/geometry.png", type="string",\
    help="Input directory for geometry. Defines the geometry of the map. [default:%default]")
    
    optional_options.add_option('-t', '--temperature_image',\
    default="../data/temperature_soft.png", type="string",\
    help="Input directory for temperature image. Defines which parts are hot. [default:%default]")
    
    optional_options.add_option('-f', '--flow_image',\
    default="../data/vent_flow_soft.png", type="string",\
    help="Input directory for flow image. Defines regions of upwards vent flow. [default:%default]")
    
    optional_options.add_option('-d', '--edge_image',\
    default="../data/vent_edges.png", type="string",\
    help="Input directory for edge image. Defines the edges of the " + 
    "map. (should match geometry) [default:%default]")
    
    optional_options.add_option('-v', '--visual_overlay_image',\
    default='../data/visual_overlay.png', type="string",\
    help="Input directory for visual overlay image. Purely graphical " + 
    "overlay to make things pretty (no effect on simulation). [default:%default]")
    
    optional_options.add_option('-c', '--color_scheme',\
    default='viridis' , type="string",\
    help="Color scheme. [default:%default]")
    
    optional_options.add_option('-o', '--output_movie_file',\
    default='./simulation_video.mp4', type="string",\
    help="Output directory for simulation video. [default:%default]")
    
    optional_options.add_option('-s', '--random_seed',\
    default=100, type="int",\
    help="The random seed number. [default:%default]")
    
    parser.add_option_group(optional_options)
    return parser


def main():
    """Run the simulation and save an output movie"""
    parser = make_option_parser()
    opts, args = parser.parse_args()
	
    ###Simulation parameters

    base_carrying_capacity = opts.base_carrying_capacity #how many microbes fit in a square?
    edge_carrying_capacity = opts.edge_carrying_capacity #how many more if the square is an edge?
    simulation_length = opts.simulation_length  #how many timesteps shall we run the simulation?

    
    input_image = opts.input_image  #defines the geometry of the map
    temperature_image = opts.temperature_image #defines which parts are hot
    flow_image = opts.flow_image #defines regions of upwards vent flow
    edge_image = opts.edge_image #defines the edges of the map (should match geometry)
    
    ### Visual and output parameters

    visual_overlay_image = opts.visual_overlay_image #purely graphical overlay to make things pretty (no effect on simulation)
    color_scheme = opts.color_scheme
    output_movie_file = opts.output_movie_file

    ### Load all the maps used in the simulation and check that they match

    geometry=load_map_image(input_image)  #read in the .png image file as a numpy array   
    height,width = geometry.shape #All other maps must have this shape!

    #Load the other maps, checking that they *EXACTLY* match in height,width
    temperature_map = load_map_image(temperature_image,height,width)
    temperature_map = temperature_map*300-80 
    flow_map = load_map_image(flow_image,height,width)
    edge_map = load_map_image(edge_image,height,width)
    
    #I don't use the load_map_image function here 
    #because I want to keep this RGB, not grayscale
    visual_overlay = imread(visual_overlay_image)

    #Set up the range of values to show in the heatmap 
    #If we change how data is output (e.g. log2 vs. log10) we have to change this
    min_heatmap_value = 0.0
    max_heatmap_value =np.log10(edge_carrying_capacity+base_carrying_capacity)
    
    #Manually set up three microbes
    microbes = {}

    starting_distribution =  np.zeros((height,width),dtype=float)

    name = "Planktonic microbe"
    planktonic_migration_zone = ((0.0,0.20),(0.0,1.0))
    traits = {"temperature_optimum":35.0}
    planktonic_microbe = Microbe(name=name,traits=traits,
      migration_zones=[planktonic_migration_zone],\
      distribution=starting_distribution,migration_rate=10**4)
    microbes[name]=planktonic_microbe

    name="Deep microbe"
    vent_migration_zone = ((0.8,1.0),(0.0,1.0))
    traits = {"temperature_optimum":200.0}
    deep_microbe = Microbe(name=name,traits=traits,
      migration_zones=[vent_migration_zone],\
      distribution=starting_distribution,migration_rate=10**4)
    microbes[name] = deep_microbe
    
    name="Deep microbe (lower temp)"
    vent_migration_zone = ((0.8,1.0),(0.0,1.0))
    traits = {"temperature_optimum":100.0}
    deep_microbe = Microbe(name=name,traits=traits,
      migration_zones=[vent_migration_zone],\
      distribution=starting_distribution,migration_rate=10**4)
    microbes[name] = deep_microbe
    
    #Set up the simulation
    simulation = Simulation(geometry,\
      base_carrying_capacity=base_carrying_capacity,\
      edge_carrying_capacity=edge_carrying_capacity,\
      upwards_flow_map=flow_map,\
      edge_map=edge_map,\
      temperature_map=temperature_map,microbes=microbes)

    #Generate an update function that will update the plot on each timestep of the simulation    
    microbe_color_map_names = {"Deep microbe":'Oranges',\
      'Deep microbe (lower temp)':'Greens','Planktonic microbe':'Purples'}
    microbe_color_maps ={}
    for microbe,cmap_name in microbe_color_map_names.items():
        cmap = make_transparent_cmap(get_cmap(cmap_name))
        microbe_color_maps[microbe] = cmap
    #Make the figure (but don't animate it yet)
    fig,vent_plots = initiate_simulation_plot(geometry,visual_overlay,\
      microbe_color_maps=microbe_color_maps,\
      min_heatmap_value=min_heatmap_value,\
      max_heatmap_value=max_heatmap_value)
    
    update_fn = make_plot_updater(vent_plots,simulation,\
    width=width,height=height)
    
    #The code for writing animation files is based on the
    #matplotlib tutorial here:
    # http://matplotlib.org/examples/animation/basic_example_writer.html
    
    #Run the simulation and make the movie
    ani = animation.FuncAnimation(fig, update_fn,\
      frames=simulation_length,interval=50, blit=True)
    
    #Set up a Writer object, which writes the movie to a file (requires ffmpeg)
    ffmpeg_writer = animation.writers['ffmpeg']
    writer = ffmpeg_writer(fps=15, metadata=dict(artist="None"), bitrate=1800)
    #Write the movie to an output file
    ani.save(output_movie_file, writer=writer)

if __name__ == "__main__":
    main()
