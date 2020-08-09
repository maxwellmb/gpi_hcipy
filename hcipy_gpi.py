from hcipy import *
from gpipsfs import *

def make_geminisouth_aperture(normalized=False, with_spiders=True):
    '''Make the Gemini South aperture.

    Parameters
    ----------
    normalized : boolean
        If this is True, the outer diameter will be scaled to 1. Otherwise, the
        diameter of the pupil will be 6.5 meters.
    with_spiders: boolean
        If this is False, the spiders will be left out.

    Returns
    -------
    Field generator
    The Gemini South aperture.
    '''
    pupil_diameter = 7.701 #m - Defined by the secondary
    spider_widths = np.array([0.014,0.01,0.01,0.01])
    support_angles = [90-43.10, 90+43.10, 270-43.10, 270+43.10]
    central_obscuration_ratio = 1.2968/pupil_diameter #From gpipsfs the secondary is 1.2968
    support_offset_y = [0.2179, -0.2179,  -0.2179,   0.2179]

    spider_offset = [0,0.34] #m

    if normalized:
        spider_widths /= pupil_diameter
        support_offset_y = [x / pupil_diameter for x in support_offset_y]
        pupil_diameter = 1.0

    obstructed_aperture = make_obstructed_circular_aperture(pupil_diameter, central_obscuration_ratio)

    if not with_spiders:
        return obstructed_aperture

    spider1 = make_spider_infinite([support_offset_y[0],0], support_angles[0], spider_widths[0])
    spider2 = make_spider_infinite([support_offset_y[1],0], support_angles[1], spider_widths[1])
    spider3 = make_spider_infinite([support_offset_y[2],0], support_angles[2], spider_widths[2])
    spider4 = make_spider_infinite([support_offset_y[3],0], support_angles[3], spider_widths[3])   

    def func(grid):
        output = obstructed_aperture(grid) * spider1(grid) * spider2(grid) * spider3(grid) * spider4(grid)
        return output.T
    return func

def make_gpi_lyot(normalized=False, with_spiders=True):
    '''Make the GPI H-band lyot mask '080m12_04'

    Parameters
    ----------
    normalized : boolean
        If this is True, the outer diameter will be scaled to 1. Otherwise, the
        diameter of the pupil will be 6.5 meters.
    with_spiders: boolean
        If this is False, the spiders will be left out.

    Returns
    -------
    Field generator
    The GPI H-band lyot mask
    '''
    
    #Many of the quantities below are defined in GPIPSFs in mm and then projected on the the primary. 
    primary_diameter = 7.7701
    magnification = primary_diameter/0.009825   # meters at primary/meters at Lyot
    
    pupil_diameter = 2*4.786*1e-3*magnification
    central_obscuration_ratio = 1.388/4.786
    spider_widths = np.array([0.404*1e-3*magnification for x in range(4)])
    support_angles = [90-43.10, 90+43.10, 270-43.10, 270+43.10]
    support_offset_y = [0.2179, -0.2179,  -0.2179,   0.2179]

    spider_offset = [0,0.34] #m

    if normalized:
        #When normalizing here, we want to normalize by the primary diameter
        spider_widths /= primary_diameter
        support_offset_y = [x / primary_diameter for x in support_offset_y]
        pupil_diameter /= primary_diameter

    obstructed_aperture = make_obstructed_circular_aperture(pupil_diameter, central_obscuration_ratio)

    if not with_spiders:
        return obstructed_aperture

    spider1 = make_spider_infinite([support_offset_y[0],0], support_angles[0], spider_widths[0])
    spider2 = make_spider_infinite([support_offset_y[1],0], support_angles[1], spider_widths[1])
    spider3 = make_spider_infinite([support_offset_y[2],0], support_angles[2], spider_widths[2])
    spider4 = make_spider_infinite([support_offset_y[3],0], support_angles[3], spider_widths[3])   

    def func(grid):
        return obstructed_aperture(grid) * spider1(grid) * spider2(grid) * spider3(grid) * spider4(grid)
    return func

def make_gpi_apodizer_no_sat_spots(normalized=True):
    '''Make the GPI Apodizer (w/o the diffraction grid)

    Parameters
    ----------
    normalized : boolean
        If this is True, the outer diameter will be scaled to 1. Otherwise, the
        diameter of the pupil will be 6.5 meters.

    Returns
    -------
    Field generator
    The gpi apodizer
    '''
    gpi_apodizer = GPI_Coronagraphic_Apodizer()
    
    pupil_diameter = 7.701 #m
    
    def func(grid):
        if grid.is_('cartesian'):
            if grid.is_separated:
                x,y = grid.separated_coords
                
                if normalized:
                    r = np.sqrt((x[np.newaxis,:]*pupil_diameter)**2+(y[:,np.newaxis]*pupil_diameter)**2)
                else:
                    r = np.sqrt((x[np.newaxis,:])**2+(y[:,np.newaxis])**2)
                print(np.max(x))
                transmission=gpi_apodizer._apod_interpolator(r).ravel()
            else:
                x,y = grid.coords
                r = np.sqrt(x**2+y**2)
                if normalized:
                    r *= (pupil_diameter)
                transmission=gpi_apodizer._apod_interpolator(r)
                
        else:
            if normalized:
                transmission = gpi_apodizer._apod_interpolator(grid.r*(pupil_diameter))
            else:
                transmission = gpi_apodizer._apod_interpolator(grid.r)
        
        return Field(transmission.astype('float'),grid)
            
    return func

def make_gpi_apodizer(normalized=True):
    '''Make the GPI Apodizer (w/o the diffraction grid)

    Parameters
    ----------
    normalized : boolean
        If this is True, the outer diameter will be scaled to 1. Otherwise, the
        diameter of the pupil will be 6.5 meters.

    Returns
    -------
    Field generator
    The gpi apodizer
    '''
    

    apod_parameters = (7.5, 585) #H-band apodizer parameters
    
    primary_diameter = 7.701
    magnification = primary_diameter/.009825   # meters at primary/meters at Lyot
    
    width = apod_parameters[0]*1e-6*magnification
    spacing = apod_parameters[1]*1e-6*magnification
    diag_spacing = spacing * np.sqrt(2)
    start = -primary_diameter/2
    
    if normalized:
        width /= primary_diameter
        spacing /= primary_diameter
        diag_spacing /= primary_diameter
        start /= primary_diameter
        pupil_diameter = 1.
    else:
        pupil_diameter = primary_diameter
        
    
    def func(grid):
        to_return = make_gpi_apodizer_no_sat_spots(normalized=normalized)(grid)
        nlines = 11
        
        for n in np.arange(nlines*2+1)-nlines:
            to_return *= make_spider_infinite([pupil_diameter/2+n*diag_spacing,pupil_diameter/2], 45, width)(grid)
            to_return *= make_spider_infinite([pupil_diameter/2,-pupil_diameter/2+n*diag_spacing], -45, width)(grid)
            
        return to_return
    return func

