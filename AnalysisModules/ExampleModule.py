import numpy as np
import scipy as sp


def AwesomeAnalysisModule(event, some_optional_number=7,
                          another_option_to_consider=True,
                          option3=np.array([0, 3.14159], dtype=np.float64)
                          ):
    # First some housekeeping:
    #
    # Define some default values to output -- if things go wrong down the
    # road you might need these, but make sure they won't be confused with
    # real answers -- -1 or -99 are common possibilities
    # Note - all the outputs should be np datatypes of fixed length - ndarrays
    # of any dimension are ok, as are scalars, as long as the shape is fixed
    # (i.e. shape of array is same for every event)
    #
    # Exception -- some modules will have one variable length dimension that's
    # common to all the keys -- e.g. # of PMT hits in the PMTtraces analysis.
    # In this case, that dimension should be the first dimension in each key.

    default_output = dict(t0=np.float64(-1),
                          bubblebigness=np.array([-99, -99], dtype=np.float64),
                          bubblebrightness=np.float64(-1))

    # Check that event loaded what you needed -- in this case we'll
    # say we need fastDAQ and PMTtraces
    if not (event['fastDAQ']['loaded'] and event['PMTtraces']['loaded']):
        return default_output

    # and now that you have what you need, make your real output
    output_dict = dict()
    output_dict['t0'] = GreatT0Finder(event['fastDAQ'], option3)

    if another_option_to_consider:
        (output_dict['bubblebigness'], output_dict['bubblebrightness']) =\
            SomeOtherFunction(event, some_optional_number, output_dict['t0'])
    else:
        output_dict['bubblebigness'] = np.float64(7.7)
        output_dict['bubblebrightness'] = np.random.randn(1)

    return output_dict


def GreatT0Finder(event, someoption):
    return np.float64(3)


def SomeOtherFunction(event, anotheroption, t0):
    return (np.float64(1.1), np.float64(88))
