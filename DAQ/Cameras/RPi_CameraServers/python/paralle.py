#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 11:59:50 2020

@author: pi
"""

import os                                                                       
from multiprocessing import Pool                                                
                                                                                

class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""

class Timer:
    def __init__(self):
        self._start_time = None

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError( f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None
        print(f"Elapsed time: {elapsed_time:0.4f} seconds")



t = Timer()
t.start()
processes = ('process1.py')                                    
other = ('process2.py','process3.py')
final = ('process4.py')
                                                                                
def run_process(process):                                                             
    os.system('python {}'.format(process))                                       
                                                                                
                                                                                
pool = Pool(processes=4)                                                        
pool.map(run_process, processes) 
pool.map(run_process, other)
pool.map(run_process, final)
t.stop()

