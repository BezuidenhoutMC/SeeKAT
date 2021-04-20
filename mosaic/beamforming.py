import numpy as np
import datetime, logging
import katpoint
import coordinate as coord
from interferometer import InterferometryObservation
from tile import ellipseCompact, ellipseGrid
from plot import plotPackedBeam, plotBeamContour, plotBeamWithFit, plot_interferometry, plot_overlap
from beamshape import calculateBeamOverlaps
from utilities import normInverse

logger = logging.getLogger(__name__)

class PsfSim(object):
    """
    Class for simulation of beam shape.

    arguments:
    antenna -- a list of antenna coordinates whose element is
        in the order of latitude(deg), longitude(deg), altitude(meter)
    frequencies -- the central frequency of the observation in Hz

    """
    reference_antenna = coord.Antenna('ref', (-30.71106, 21.44389, 1035))
    '''speed Of Light'''
    sol = 299792458

    def __init__(self, antennas, frequencies):
        """
        constructor of the PsfSim class.

        """
        waveLengths = float(self.sol)/np.array(frequencies)
        if waveLengths.shape == (1,): waveLengths = waveLengths[0]
        self.observation = InterferometryObservation(self.reference_antenna,
                            waveLengths)
        self.antennas = PsfSim.check_antennas(antennas)

    @staticmethod
    def check_antennas(antennas):
        """
        check the type of the inputs. if they are katpoint objects,
        then extract the latitude, longitude, elevation information

        arguments:
        antennas -- either can be a list of katpoint antenna objects or
                    a list of [latitude, longitude, elevation]

        return:
        a list of antenna geographic coordinates in the order of
        [latitude, longitude, elevation]

        """
        def from_katpoint_list(antennas):
            antenna_list = []
            for antenna in antennas:
                antenna_list.append([np.rad2deg(antenna.observer.lat),
                                    np.rad2deg(antenna.observer.lon),
                                    antenna.observer.elev])
            return np.array(antenna_list)
        antennas = np.array(antennas)
        if isinstance(antennas[0], np.ndarray):
            antenna_coordinates = antennas
            names = ["%03d" % i for i in range(len(antennas))]
        elif isinstance(antennas[0], katpoint.Antenna):
            antenna_coordinates = from_katpoint_list(antennas)
            names = [ant.name for ant in antennas]
        elif isinstance(antennas[0], str):
            katpoint_antennas = [katpoint.Antenna(i) for i in antennas]
            antenna_coordinates = from_katpoint_list(katpoint_antennas)
            names = [ant.name for ant in katpoint_antennas]
        else:
            raise Exception("Antennas are passed in unknown format")

        antenna_objects = [coord.Antenna(names[idx], coordinate)
                for idx, coordinate in enumerate(antenna_coordinates)]
        return antenna_objects

    @staticmethod
    def check_source(source):
        """
        check the type of the inputs. if it is a katpoint object,
        then extract the RA, DEC information

        arguments:
        source -- either can be a katpoint target object or [RA, DEC]


        return:
        a coordinate as [RA, DEC]

        """
        if (isinstance(source, np.ndarray) or isinstance(source, list) or
                isinstance(source, tuple)):
            return source
        elif isinstance(source, katpoint.Target):
            ra = source.body._ra
            dec = source.body._dec
            return (str(ra), str(dec))
        else:
            raise Exception("source are passed in unknown format")

    def get_beam_shape(self, source, time, beam_number = 400, beam_size = None):
        """
        return the beamshape of current oservation parameters
        assuming the beam is roughly a ellipse.

        arguments:
        source -- the boresight location or telescope pointing
            in the order of RA(deg), DEC(deg)
        time -- the observation time in datatime object or epoch seconds


        return:
        a beamshape object contain properties of semi-major axis,
            semi-mino axis and orientation all in degree
        """
        if len(self.antennas) < 3:
            raise "the number of antennas should be not less than 3"
        bore_sight = PsfSim.check_source(source)
        self.observation.setBoreSight(bore_sight)
        self.observation.setObserveTime(time)
        self.observation.setBeamNumber(beam_number)
        if beam_size != None:
		self.observation.setAutoZoom(False)
		self.observation.setBeamSizeFactor(beam_size)
        self.observation.createContour(self.antennas)
        axisH, axisV, angle, image_range = self.observation.getBeamAxis()
        horizon = np.rad2deg(self.observation.getBoreSight().horizontal)
        psf = self.observation.getPointSpreadFunction()
        resolution = self.observation.getResolution()
        bore_sight_object = self.observation.getBoreSight()
        return BeamShape(axisH, axisV, angle, psf, self.antennas,
                bore_sight_object, self.reference_antenna, horizon, resolution)


class BeamShape(object):
    """
    Class of  the BeamShape object contain properties of the beamshape

    arguments:
    axisH -- length of the semi-major axis in degree
    axisV -- length of the semi-minor axis in degree
    angle -- orientation of the angle in degree

    """
    def __init__(self, axisH, axisV, angle, psf, antennas, bore_sight,
            reference_antenna, horizon, resolution):
        """
        constructor of the BeamShape class.

        """
        self.axisH = axisH
        self.axisV = axisV
        self.angle = angle
        self.psf = psf
        self.antennas = antennas
        self.bore_sight = bore_sight
        self.reference_antenna = reference_antenna
        self.horizon = horizon
        self.resolution = resolution

    def width_at_overlap(self, overlap):
        """
        return the half widths of the ellipse in major axis and
        minor axis direction given a overlap level.

        the relationship between one sigma and the full width maximal can be
        found in this link

        https://en.wikipedia.org/wiki/Full_width_at_half_maximum
        2.*np.sqrt(2.*np.log(2)) = 2.3548200450309493
        """
        sigmaH = self.axisH * (2./2.3548200450309493)
        sigmaV = self.axisV * (2./2.3548200450309493)

        widthH = normInverse(overlap, 0, sigmaH)
        widthV = normInverse(overlap, 0, sigmaV)

        return widthH, widthV

    def plot_psf(self, filename, shape_overlay = False):
        """
        plot the point spread function

        arguments:
        filename --  name and directory of the plot
        shape_overlay -- whether to add the shape overlay on the psf

        """
        if not shape_overlay:
            plotBeamContour(self.psf.image, self.psf.bore_sight.equatorial,
                    self.psf.image_range, filename, interpolation = True)
        else:
            widthH = self.axisH/self.resolution
            widthV = self.axisV/self.resolution
            plotBeamWithFit(self.psf.image, self.psf.bore_sight.equatorial,
                    self.psf.image_range, widthH, widthV,
                    self.angle, filename, interpolation = True)

    def plot_interferometry(self, filename):
        """
        plot the interferometry overview, including the antennas, the source

        arguments:
        filename --  name and directory of the plot

        """
        antennas = np.array([antenna.geo for antenna in self.antennas])
        plot_interferometry(antennas, self.reference_antenna.geo, self.horizon, filename)


class Overlap(object):
    """
    Class of overlap object contain a overlap calculation result

    arguments:
    metrics -- a measurement of overlap in a gridded area
    mode -- the mode used in the overlap calculation

    """
    def __init__(self, metrics, mode):
        """
        constructor of the Tiling class.

        """

        self.metrics  = metrics
        self.mode = mode

    def plot(self, filename):
        """
        plot the overlap result in specified filename

        """

        plot_overlap(self.metrics, self.mode, filename)

    def calculate_fractions(self):
        """
        calculation the occupancy of different overlap situation in a region.
        This method only works if the overlap is calculated in "counter" mode.
        overlap situation include: overlap, non-overlap, empty

        return
        overlapped: faction of the region where beams are overlapped
        non_overlapped: faction of the region inside a single beam
        empty: fraction of region where is not covered by any beam

        """

        if self.mode != "counter":
            raise "the fraction calculation is only supportted in counter mode"
        overlap_counter = self.metrics
        overlap_grid = np.count_nonzero(overlap_counter > 1)
        non_overlap_grid = np.count_nonzero(overlap_counter == 1)
        empty_grid = np.count_nonzero(overlap_counter == 0)
        point_num = overlap_grid+non_overlap_grid+empty_grid
        overlapped, non_overlapped, empty = np.array([overlap_grid, non_overlap_grid,
                empty_grid])/float(point_num)
        return overlapped, non_overlapped, empty


class Tiling(object):
    """
    Class of Tiling object contain a tiling result

    arguments:
    coordinates -- tiling coordinates as a list of [RA, DEC] in degree
    beamShape -- beamShape object
    raidus -- the raidus of the entire tiling
    overlap -- how much overlap between two beams, range in (0, 1)
    """
    def __init__(self, coordinates, beam_shape, radius, overlap):
        """
        constructor of the Tiling class.

        """

        self.coordinates = coordinates
        self.beam_shape = beam_shape
        self.tiling_radius = radius
        self.beam_num = len(coordinates)
        self.overlap = overlap

    def plot_tiling(self, filename, overlap = None, index = False):
        """
        plot the tiling pattern with specified file name.

        arguments:
        filename --  filename of the plot, the format and directory can be
        specified in the file name such as  "plot/pattern.png" or "pattern.pdf"
        overlap -- the overlap ration between beams,
                   default is None(using the tilling overlap)
        index -- wather to show the index of the beam, default is False
        """
        if overlap is None:
            widthH, widthV = self.beam_shape.width_at_overlap(self.overlap)
        elif overlap == 0.5:
            widthH, widthV = self.beam_shape.axisH, self.beam_shape.axisV
        else:
            widthH, widthV = self.beam_shape.width_at_overlap(overlap)
        maxDistance = np.max(
                np.sqrt(np.sum(np.square(self.coordinates), axis = 1)))
        upper_left_pixel = [-maxDistance, maxDistance] # x,y
        bottom_right_pixel = [maxDistance, -maxDistance] # x,y

	#print("ul: "+str(upper_left_pixel))
	#print("br: "+str(bottom_right_pixel))
	#print("center: "+str(self.beam_shape.bore_sight.equatorial))

        coordinates_equatorial = coord.convert_pixel_coordinate_to_equatorial(
            [upper_left_pixel, bottom_right_pixel], self.beam_shape.bore_sight.equatorial)
        equatorial_range = [
            coordinates_equatorial[0][0], coordinates_equatorial[1][0], # left, right
            coordinates_equatorial[0][1], coordinates_equatorial[1][1]] # up, bottom

	#Tiaan mod
	width = equatorial_range[1] - equatorial_range[0]
	height = equatorial_range[3] - equatorial_range[2]
	print("Size of primary beam: "+str(width) + " deg wide, " + str(height) + " deg heigh.")	

        pixel_range = [upper_left_pixel[0], bottom_right_pixel[0], # left, right
                       upper_left_pixel[1], bottom_right_pixel[1]] # up, bottom
        plotPackedBeam(self.coordinates, self.beam_shape.angle, widthH, widthV,
            self.beam_shape.bore_sight.equatorial, equatorial_range, pixel_range,
            self.tiling_radius, fileName=filename, index = index)

    def get_equatorial_coordinates(self):
        """
        convert pixel coordinates to equatorial coordinates

        return:
        coordinates_equatorial --  tiling coordinates in equatorial frame
        """
        coordinates_equatorial = coord.convert_pixel_coordinate_to_equatorial(
               self.coordinates, self.beam_shape.bore_sight.equatorial)
        return coordinates_equatorial

    def get_beam_size(self):
        """
        get the size of the beam in equatorial metric

        return:
        width1, width2 --  semi-majors of the beam in degree in equatorial plane
        """

        axis1, axis2 = self.beam_shape.width_at_overlap(self.overlap)
        width1, width2 = coord.convert_pixel_length_to_equatorial(axis1, axis2,
                self.beam_shape.angle, self.beam_shape.bore_sight.equatorial)
        return width1, width2

    def plot_sky_pattern(self, filename):
        """
        plot the ksy pattern with specified filename

        """

        heats = self.calculate_overlap("heater", new_beam_shape = None)
        heats.plot(filename)

    def calculate_overlap(self, mode, new_beam_shape = None):
        """
        calculate overlap of the tiling pattern.

        arguments:
        mode -- mode of the calculation,
                "heater" will calculate the tiling pattern as sky temperature,
                "counter" will calculate the counts of the overlap regions,
                          non-overlap regions and empty regions.
        filename --  filename of the plot, the format and directory can be
                specified in the file name such as  "plot/pattern.png"

        return:
        overlap_counter -- counter when choose "counter" mode
        overlapHeater -- heater when choose "heater" mode
        """
        if new_beam_shape is None:
            beam_shape = self.beam_shape
        else:
            beam_shape = new_beam_shape
        overlap_metrics = calculateBeamOverlaps(
                self.coordinates, self.tiling_radius, beam_shape.axisH,
                beam_shape.axisV, beam_shape.angle, self.overlap, mode)
        overlap = Overlap(overlap_metrics, mode)
        return overlap


def generate_nbeams_tiling(beam_shape, beam_num, overlap = 0.5, margin=None):
    """
    generate and return the tiling.
    arguments:
    beam_shape -- beam_shape object
    beam_num -- number of beams to tile
    overlap -- how much overlap between two beams, range in (0, 1)
    margin -- the maximum difference in beam number between required and generated
             increase this number will make the tiling process faster, while
             decrease this number will make the tiling process slower or even fail
             the number of actual generated beams will always larger then required
             default is 5% of the required beam number
    return:
    tiling -- tiling coordinates in a list of pixel coordinates pairs in degree
    """
    if margin is None:
        margin = int(round(beam_num * 0.05))
    widthH, widthV = beam_shape.width_at_overlap(overlap)
    tiling_coordinates, tiling_radius = ellipseCompact(
            beam_num, widthH, widthV, beam_shape.angle, margin,
            seed = beam_shape.bore_sight.equatorial[0])

    tiling_obj = Tiling(tiling_coordinates, beam_shape, tiling_radius, overlap)
    width1, width2 = tiling_obj.get_beam_size()

    logger.info("tiling: overlap: {},radius: {:.3g} in pixel plane".format(
                overlap, tiling_radius))
    logger.info("tiling: width1: {:.3g} arcsec, width2: {:.3g} arcsec "
                "in equatorial plane".format(width1.arcsec, width2.arcsec))

    return tiling_obj

def generate_radius_tiling(beam_shape, tiling_radius, overlap = 0.5):
    """
    return the tiling inside a specified region
    arguments:
    beam_shape -- beam_shape object
    tilingRadius -- the radius of the region to tile
    overlap -- how much overlap between two beams, range in (0, 1)
    return:
    tiling -- tiling coordinates in a list of pixel coordinates pairs in degree
    """
    widthH, widthV = beam_shape.width_at_overlap(overlap)
    tiling_coordinates = ellipseGrid(
            tiling_radius, widthH, widthV, beam_shape.angle)

    tiling_obj = Tiling(tiling_coordinates.T, beam_shape, tiling_radius, overlap)

    return tiling_obj

def dict_to_ordered_list(dict_obj):
    ordered_list = []
    for key in sorted(dict_obj.iterkeys()):
        ordered_list.append(dict_obj[key])
    return ordered_list

def dict_to_antenna_ordered_list(dict_obj, antennas, pol='h'):
    ordered_list = []
    for antenna in antennas:
        antenna_key = "{}{}".format(antenna.name, pol)
        ordered_list.append(dict_obj[antenna_key])

    return ordered_list

class DelayPolynomial(object):
    """
    Class for generation of  delay polynomial

    arguments:
    antennas -- a list of antenna objects or coordinate in csv format
    targets -- a list of beam location in equatorial coordinates
    frequencies -- a list of frequencies on which the polynomail is calculated in Hz
    reference -- the reference antenna for delay calculation
    """
    def __init__(self, antennas, bore_sight, targets, reference):
        """
        constructor of the Delay Polynomial class.

        """
        self.antennas = antennas
        self.targets = DelayPolynomial.check_targets(targets)
        self.frequency = 1.4e9
        self.reference = reference
        self.bore_sight = DelayPolynomial.check_targets([bore_sight,])[0]

    @staticmethod
    def check_targets(targets):
        """
        check the target data type, the arguments will be converted to katpoint
            object if they are not.

        arguments:
        targets -- a list of target objets in the format of
            katpoint target object or set of [longitude, latitude, altitude]

        return:
            targets in katpoint object

        """
        if isinstance(targets[0], katpoint.Target):
            return targets
        else:
            return DelayPolynomial.make_katpoint_target(targets)

    @staticmethod
    def check_time(time):
        """
        check the the data type of the time value. If the values are datetime
            objects, they will be converted to seconds

        arguments:
        time -- epoch seconds or datetime objects

        return:
        time in epoch seconds
        """
        if type(time) != int and type(time) != float:
            return coord.datetimeToEpoch(time)
        else:
            return time

    @staticmethod
    def make_katpoint_target(sources):
        """
        check the the data type of the source. If the values are in (RA, DEC) pair,
        they will be converted to katpoint target object

        arguments:
        source -- source in either (RA, DEC) pairs or katpoint target objects

        return:
        katpoint target objects
        """

        targets = []
        for source in sources:
            target_string = ",".join(['radec',
                            coord.angleToHour(source[0]),
                            coord.angleToDEC(source[1])])
            targets.append(katpoint.Target(target_string))
        return targets

    def get_delay_polynomials(self, epoch, duration=10.0):
        """
        calculate and return the polynomials

        Arguments:
        timestamp -- the observation time in datatime object or epoch seconds
        duration -- the duration in which the polynomial is calcuated

        return:
        polynomials in the order of beam, antenna, (delay, rate)

        """
        timestamp = DelayPolynomial.check_time(epoch)
        antennaObjectList = self.antennas
        timestamp = (timestamp, timestamp + duration)

        target_array = []
        for target in self.targets:
            dc = katpoint.DelayCorrection(self.antennas,
                    self.reference, self.frequency)
            delay, phase = dc.corrections(target, timestamp)
            delayArray = np.array(dict_to_antenna_ordered_list(
                        delay, self.antennas))
            """
            [:, 0, :]: only take first rate output
            """
            target_array.append(delayArray[:,0,:])
        target_array = np.array(target_array)
        """
        subtract the boresight beam form the offset beams
        """
        dc = katpoint.DelayCorrection(self.antennas,
            self.reference, self.frequency)
        delay, phase = dc.corrections(self.bore_sight, timestamp)
        bore_sight_delay = np.array(dict_to_antenna_ordered_list(
                    delay, self.antennas))[:,0,:]

        target_array = target_array - bore_sight_delay
        return target_array
