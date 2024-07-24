#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup

requirements = ['Click>=7.0', 'matplotlib >= 3.1.2', 'seaborn>=0.9.0','pyyaml>=5.1.2', 'plotly>=4.3.0',
                'numpy>=1.12.1', 'pandas>=0.24.2', 'biopython==1.79', 'scipy>=1.2.1', 'statsmodels',
                'jinja2', 'tqdm>=4.66.0', 'mpmath>=1.3.0']

setup(
    author="Guy Assa",
    author_email='guyassa@gmail.com',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: Other/Proprietary License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.8',

    ],
    description="CRISPECTOR2 - Genome Editing Analysis Tool, with allele extension",
    entry_points={
        'console_scripts': [
            'crispector = crispector2.cli:main',
        ],
    },
    python_requires='>=3.8',
    install_requires=requirements,
    include_package_data=True,
    keywords='crispector2',
    name='crispector2',
    package_dir={'crispector2': 'crispector2'},
    packages=['crispector2', 'crispector2.algorithm', 'crispector2.allele', 'crispector2.config',
              'crispector2.input_processing', 'crispector2.modifications', 'crispector2.report', 'crispector2.utils'],
    url='https://github.com/theAguy/crispector2',
    version='2.1.2',
    zip_safe=False,
)