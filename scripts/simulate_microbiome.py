"""
=================
An animated image
=================

This example demonstrates how to animate an image.
"""
import numpy as np
from matplotlib.pyplot import imshow,show,figure,colorbar
from matplotlib.image import imread

import matplotlib.animation as animation
from random import randint,random
from copy import copy

class Simulation():
    def __init__(self,geometry=None,width=500,height=500,base_carrying_capacity=10**12,edge_carrying_capacity=0,
      upwards_flow_map = None,edge_map=None,temperature_map = None):
        if geometry is None: #make an empty map
            self.Geometry = np.zeros((height,width),dtype=float)
        else:
            self.Geometry = geometry
            height,width = geometry.shape
        self.TimeStep = 0
        self.Width = width
        self.Height = height
        
        self.Microbes = np.zeros((height,width),dtype=float)
        
        if temperature_map is None:
            self.Temperatures = np.zeros((height,width),dtype=float)
        else:
            self.Temperatures = temperature_map
        if upwards_flow_map is None:
            self.UpwardsFlowMap = np.zeros((height,width),dtype=float)
        else:
            self.UpwardsFlowMap = upwards_flow_map

        if edge_map is None:
            self.EdgeMap = np.zeros((height,width),dtype=float)
        else:
            self.EdgeMap = edge_map    
        
        self.CarryingCapacityMap = base_carrying_capacity + (self.EdgeMap * edge_carrying_capacity) 
        self.SpreadRates = {"up":0.025,"down":0.15,"right":0.05,"left":0.05}
        self.MigrationRate = 10 #number of timesteps on which to introduce a migrant
        
    def update(self,verbose=True):
        """Run one round of the simulation"""
        self.simulate_vent_flow()
        self.simulate_spread() #along surfaces

        planktonic_migration_zone = ((0.0,0.20),(0.0,1.0))
        deep_vent_migration_zone = ((0.95,1.0),(0.4,0.6))

        self.simulate_migration(self.MigrationRate,migration_probability=0.50,\
          migration_zone=planktonic_migration_zone)
        
        self.simulate_migration(self.MigrationRate//8,\
          migration_probability=0.20,migration_zone = deep_vent_migration_zone)       
        self.simulate_growth()
        self.simulate_mortality() #from temperature for now
        self.apply_carrying_capacity() 
        self.TimeStep += 1

        if self.TimeStep % 100 == 0 and verbose:
            print("Simulating timestep {}".format(self.TimeStep))

    def simulate_migration(self,migration_rate=1,\
      migration_zone=((0.0,1.0),(0.0,1.0)),migration_probability=1.0,migration_number=10**2):
        """Add a microbe to a random cell"""
        (min_y,max_y),(min_x,max_x) = migration_zone
        if random() < migration_probability:
            for i in range(int(migration_rate)):
                new_microbe_x = randint(int(min_x*self.Width-1),int(max_x*self.Width-1))
                new_microbe_y = randint(int(min_y*self.Height-1),int(max_y*self.Height-1))
                self.Microbes[new_microbe_y,new_microbe_x] += migration_number
    

    def apply_carrying_capacity(self):
        """Remove micorbes above the per-square carrying capacity"""
        self.Microbes = np.minimum(self.Microbes,self.CarryingCapacityMap)
        self.Microbes[self.Microbes < 0] = 0
        self.Microbes *= self.Geometry 
        #this should be 0 or 1, so any microbes in solid rock are removed

    def simulate_growth(self,growth_rate=1.15):
        """All micorbes divide"""
        #For now assume 'good' environments have both higher
        #growth rate and carrying capacity.


        random_array = self.get_uniform_random_array()
        growth_map = growth_rate
        growth_map = growth_map * random_array + 0.5 #add 50% stochasticity to growth rates
        growth_map = np.maximum(1.0,growth_map) #no mortalit in this step
        growth_map = np.minimum(2.0,growth_map) #no growth faster than 2.0
        self.Microbes = self.Microbes*growth_map
    
    def simulate_mortality(self,intrinsic_mortality=0.001,max_temp_mortality=0.12):
        death_map = 1.0 - (intrinsic_mortality + (self.Temperatures) * max_temp_mortality)
        #Not caching this because once we add evolution all microbes will have different thermal responses
        self.Microbes = death_map*self.Microbes
 
    def simulate_vent_flow(self,vent_velocity=1.00,edge_stickiness=0.70):
        offset_array = generate_offset_array(copy(self.Microbes),direction="up")
        #Flow occurs at different rates, and microbes on surfaces 'stick'
        #we subtract the edgemap to simulate this
        random_array = self.get_uniform_random_array()+0.50
        spread_rate = (self.UpwardsFlowMap - self.EdgeMap*edge_stickiness) * vent_velocity*random_array 
        offset_array = offset_array * self.Geometry * self.UpwardsFlowMap
        #Microbes can't spread into solid rock
        #apply geometry to the offset array (handled by multiplying by the geometry map)
        #any value times 0 = 0, so black areas should not be included

        new_data = self.Microbes + offset_array * spread_rate
        new_data = new_data - self.Microbes * spread_rate #conservation of microbes
        self.Microbes = new_data
    
    def get_uniform_random_array(self):
        """Return a uniform random array between 0 and 1 """
        random_array = np.random.random((self.Height,self.Width))
        return random_array

    def simulate_spread(self,edge_bonus=0.10):
        """Spread microbes to neighboring cells
        edge_bonus -- if >0.0, allow microbes on edges to spread faster
          (in any direction)
        """
        #Check that the spread does not increase cell numbers
        starting_total = np.sum(self.Microbes)
        for direction in ["up","down","left","right"]:
            spread_rate = self.SpreadRates[direction]
            offset_array = generate_offset_array(copy(self.Microbes),direction=direction)
            #Microbes can't spread into solid rock
            #apply geometry to the offset array
            #any value times 0 = 0, so black areas should not be included
            offset_array = offset_array * self.Geometry
            
            #Optionally adjust the spread rate to be higher/lower on edges
            adj_spread_rate = spread_rate + edge_bonus * self.EdgeMap

            #Add microbes to new cells and subtract them from the old cell
            new_data = self.Microbes + offset_array * adj_spread_rate  
            new_data = new_data - self.Microbes * adj_spread_rate#conservation of microbes

            self.Microbes = new_data


    def get_data(self,scaling ="log2"):
        """Return simulation data as a numpy array
        scaling -- if set, rescale the data for visualization
        """
     
        if scaling == 'linear':
            return self.Microbes
        elif scaling == 'log2':
            eps = 0.1/self.CarryingCapacityMap
            #add a tiny epsilon value to avoid divide by zero errors
            return np.log2(self.Microbes + eps)

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

def make_plot_updater(plot,simulation,width=500,height=500):
    """Make an updatefig function for the specified plot
    plot -- a matplotlib plot
    simulation -- a class object with an get_data method
      that returns a width * height numpy array of numbers
    width -- width in pixels
    height -- height in pixels

    """
    def update_fig(frame,*args):
        simulation.update()
        new_data = simulation.get_data()
        #Use plot.set_data to update the data in the plot
        #without making a new plot
        plot.set_data(new_data)
        #The comma is needed to make this a tuple and
        #therefore iterable
        return (plot,)

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
  color_scheme='viridis',min_heatmap_value=0.0,max_heatmap_value=10**12):
    """Return a figure and heatmap subplot for animation
    
    geometry - a numpy array reflecting the geometry of the map
    
    visual_overlay - a visual overlay that is plotted on top of the data
    
    color_scheme - what color set to use to convert raw numbers (e.g. how many microbes)
      into a color to plot. Examples include 'jet' and 'viridis'. 
      For valid values, see https://matplotlib.org/examples/color/colormaps_reference.html

    min_heatmap_value - what number should map to the first color on the color scheme?
    max_heatmap_value - what number should map to the last color on the color scheme?
    """
    #Initiate the base plot
    fig = figure()
    vent_plot = imshow(geometry,cmap=color_scheme,\
      vmin=min_heatmap_value,vmax=max_heatmap_value,animated = True)
    colorbar()
    imshow(visual_overlay)
    
    #fig is the overall figure
    #vent_plot is the specific heatmap part (we need both)
    return fig,vent_plot
 

def main():
    """Run the simulation and save an output movie"""
    #temporarily hard-coded user files
    #TODO: make this user specified through a commandline interface    
    
    ###Simulation parameteris
    base_carrying_capacity = 10**4 #how many microbes fit in a square?
    edge_carrying_capacity = 10**9 #how many more if the square is an edge?
    simulation_length = 1200  #how many timesteps shall we run the simulation?
    
    input_image = "../data/geometry.png"  #defines the geometry of the map
    temperature_image = "../data/temperature_soft.png" #defines which parts are hot
    flow_image = "../data/vent_flow_soft.png" #defines regions of upwards vent flow
    edge_image = "../data/vent_edges.png" #defines the edges of the map (should match geometry)
    
    ### Visual and output parameters
    visual_overlay_image = '../data/visual_overlay.png' #purely graphical overlay to make things pretty (no effect on simulation)
    color_scheme = 'viridis' 
    output_movie_file = './simulation_video.mp4'
    ###


    ### Load all the maps used in the simulation and check that they match

    geometry=load_map_image(input_image)  #read in the .png image file as a numpy array   
    height,width = geometry.shape #All other maps must have this shape!

    #Load the other maps, checking that they *EXACTLY* match in height,width
    temperature_map = load_map_image(temperature_image,height,width)
    flow_map = load_map_image(flow_image,height,width)
    edge_map = load_map_image(edge_image,height,width)
    
    #I don't use the load_map_image function here 
    #because I want to keep this RGB, not grayscale
    visual_overlay = imread(visual_overlay_image)

    #Set up the range of values to show in the heatmap 
    #If we change how data is output (e.g. log2 vs. log10) we have to change this
    min_heatmap_value = 0.0
    max_heatmap_value =np.log2(edge_carrying_capacity+base_carrying_capacity)
    
    #Make the figure (but don't animate it yet)
    fig,vent_plot = initiate_simulation_plot(geometry,visual_overlay,\
      color_scheme=color_scheme,min_heatmap_value=min_heatmap_value,max_heatmap_value=max_heatmap_value)

    #Set up the simulation
    simulation = Simulation(geometry,\
      base_carrying_capacity=base_carrying_capacity,\
      edge_carrying_capacity=edge_carrying_capacity,\
      upwards_flow_map=flow_map,\
      edge_map=edge_map,\
      temperature_map=temperature_map)

    #Generate an update function that will update the plot on each timestep of the simulation    
    update_fn = make_plot_updater(vent_plot,simulation,\
    width=width,height=height)
    
    
    #The code for writing animation files is based on the
    #matplotlib tutorial here:
    # http://matplotlib.org/examples/animation/basic_example_writer.html
    
    #Run the simulation and make the movie
    ani = animation.FuncAnimation(fig, update_fn,frames=simulation_length,interval=50, blit=True)
    
    #Set up a Writer object, which writes the movie to a file (requires ffmpeg)
    ffmpeg_writer = animation.writers['ffmpeg']
    writer = ffmpeg_writer(fps=15, metadata=dict(artist="None"), bitrate=1800)
    #Write the movie to an output file
    ani.save(output_movie_file, writer=writer)

if __name__ == "__main__":
    main()
