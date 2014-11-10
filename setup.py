#! /usr/bin/env python
#  -*- coding: utf-8 -*-
#
# process_control, (c) 2014, see AUTHORS.  Licensed under the GNU GPL.
from setuptools import setup, find_packages


desc = ('A simple python library to estimate first order plus deadtime models'
        'and calculate pid values.')

setup(
    name='process_control',
    version=__import__('process_control').__version__,
    author='Marco Halder',
    author_email='marco.halder@frm2.tum.de',
    license = 'GNU General Public License (GPL), Version 3',
    classifiers=[
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Operating System :: MacOS',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Software Development',
        'Topic :: Scientific/Engineering',
    ],
    url='https://github.com/p3trus/process_control',
    description=desc,
    long_description=open('README.md').read(),
    install_requires=['numpy'],
    
    packages=find_packages(),
    include_package_data=True,
)
