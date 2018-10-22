"""
    .. module:: time
      :synopsis: This module implements a function for testing / benchmarking
                 the functions. More function smay be implemented.
"""

# Third-party modules
from functools import wraps
import time


def fn_timer(function):
    """This function is a wrapper to benchmarck a function
    when trying to optimize it. It returns the total running time of the function

        Args:
            function: The function to benchmarck

        Returns:
            float: The time the function took to run

    """
    @wraps(function)
    def function_timer(*args, **kwargs):
        """It calculates the time during the run of a function"""
        t_0 = time.time()
        result = function(*args, **kwargs)
        t_1 = time.time()
        print("Total time running {:s}: {:s} seconds"\
            .format(function.__name__, str(t_1-t_0)))
        return result
    return function_timer
