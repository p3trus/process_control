#! /usr/bin/env python
#  -*- coding: utf-8 -*-
#
# process_control, (c) 2014, see AUTHORS.  Licensed under the GNU GPL.
"""A simple python library to estimate first order plus deadtime models and
calculate pid values.
"""
__version__ = '0.1.0'

import collections
import numpy as np
import numpy.linalg


def _estimate_fopdt(u, y, deadtime):
    
    b = y[1 + deadtime:]
    
    A = np.empty((len(b), 3))
    A[:, 0] = y[deadtime:-1]
    A[:, 1] = u[:-(1 + deadtime)]
    A[:, 2] = 1.
    x, residuals, rank, s = np.linalg.lstsq(A, b)
    return np.append(x, deadtime), residuals, rank, s


class FOPDT(object):
    """Implements a discrete first order plus dead time model.
    
    The corresponding difference equation is::
    
        y[n + 1] = a * y[n] + b * u[n - deadtime] + c

    with::

        a = exp(dt / time_constant)
        b = (1 - a) * gain
        c = (1 - a) * steady_state
    
    ::

        u, y = np.loadtxt('testdata.dat', unpack=True)
        #Fit model to data
        model = FOPDT.fit(u, y, deadtime=range(0, 20), dt=.5)
    
        # Simulate Output
        uo, yo = model.output(u, y0=y[0:model.deadtime + 1])

    """
    def __init__(self, a, b, c, deadtime_shift, dt=1):
        self.a = a
        self.b = b
        self.c = c
        self.deadtime_shift = int(deadtime_shift)
        self.dt = dt

    def output(self, u, y0=None):
        """Calculates the response of the FOPDT model
        to the given input vector.
        
        :param u: The input vector.
        :param y0: The initial condition.
        
        """
        deadtime = self.deadtime_shift
        y = np.zeros_like(u)
        if y0 is not None:
            y[:deadtime + 1] = y0
        for i in range(len(y) - (1 + deadtime)):
            y[i + 1 + deadtime] = y[i + deadtime] * self.a + u[i] * self.b + self.c
            
        return u, y
  
    @staticmethod
    def fit(u, y, deadtime_shift, dt=1):
        """Fits a first order plus deadtime model to the data.
        
        :param u: The input vector.
        :param y: The output vector.
        :param dt: The sampling time.
        :param deadtime_shift: The dead time in integer multiples of the
            sampling time. If the deadtime is known beforehand it can
            be used directly. Otherwise a range of deadtimes
            should be specified. The deadtime giving the best fit
            is then used.

        :returns: A :class:`~ .FOPDT` instance.
        
        """
        assert len(u) == len(y)
        if isinstance(deadtime_shift, collections.Iterable):
            fits = [_estimate_fopdt(u, y, dtime) for dtime in deadtime_shift]
            # get fit with smallest residual error
            x, residuals, rank, s = min(fits, key=lambda f: f[1])
        else:
            x, residuals, rank, s = _estimate_fopdt(u, y, deadtime_shift)
            
        return FOPDT(*x, dt=dt)
    
    @property
    def time_constant(self):
        """The time constant of the process model."""
        return - self.dt / np.log(self.a)
    
    @property
    def gain(self):
        """The gain of the process model."""
        return self.b / (1. - self.a)
    
    @property
    def steady_state(self):
        """The steady state value of the process model."""
        return self.c / (1. - self.a)
        
    @property
    def deadtime(self):
        """The process deadtime."""
        return self.deadtime_shift * self.dt
        
        
def cohen_coon(gain, time_constant, deadtime, stability_margin=1., mode='pi'):
    """calculates the cohen coon pid tuning rules.
    
    .. note::

        1. The cohen coon tuning rules are designed for non
           interactive pid controllers.
        2. Be careful to use the correct units for `gain`,
           `time_constant` and `deadtime`, e.g. if controller
           uses minutes, `deadtime` and `time_constant` should
           be in minutes to get the correct results.
           
        
    :param gain: The process gain.
    :param time_constant: The process time constant.
    :param dead_time: The process deadtime.
    :param stability_margin: Cohen coon aims for quarter amplitude damping.
        To detune the algorithm use a stability margin greater than 1.
    :param mode: Decides wether tuning rules are calculated for
        'P', 'PI', 'PD or 'PID' controllers.
    """
    def controller_gain(a, b):
        return a / gain * (time_constant / deadtime + b) / stability_margin
    def time_const(a, b, c):
        return a * deadtime * (time_constant + b * deadtime) / (time_constant + c * deadtime)
    
    if mode == 'P':
        return controller_gain(1.03, 0.34), 0., 0.
    if mode == 'PI':
        return controller_gain(0.9, 0.092), time_const(3.33, 0.092, 2.22), 0.
    if mode == 'PD':
        return controller_gain(1.24, 0.129), 0., time_const(0.27, -0.324, 0.129)
    if mode == 'PID':
        return controller_gain(1.35, 0.185), time_const(2.5, 0.185, 0.611), time_const(0.37, 0., 0.185)


def lambda_tuning(gain, time_constant, deadtime, closed_loop=3.):
    """calculates the lambda pid tuning rules.
    
    .. note::

        1. The lambda tuning rules are designed for non
           interactive  and interactive pid controllers.
        2. Be careful to use the correct units for `gain`,
           `time_constant` and `deadtime`, e.g. if controller
           uses minutes, `deadtime` and `time_constant` should
           be in minutes to get the correct results.
           
        
    :param gain: The process gain.
    :param time_constant: The process time constant.
    :param dead_time: The process deadtime.
    :param closed_loop: The closed loop time constant mulitplier.

    """
    return time_constant/(gain * (closed_loop * time_constant + deadtime)), time_constant, 0.