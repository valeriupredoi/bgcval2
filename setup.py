#!/usr/bin/env python
"""BGCVal2 installation script."""
# This script only installs dependencies available on PyPI

import json
import os
import re
import sys
from pathlib import Path

from setuptools import Command, setup

sys.path.insert(0, os.path.dirname(__file__))

PACKAGES = [
    'bgcval2',
]

REQUIREMENTS = {
    # Installation script (this file) dependencies
    'setup': [
        'setuptools_scm',
    ],
    # Installation dependencies
    # Use with pip install . to install from source
    'install': [
        'basemap>=1.3.6',
        'cartopy',
        'matplotlib',
        'nctoolkit>=0.8.7',  # use linux64 build
        'netcdf4',
        'numpy!=1.24.3',
        'pip!=21.3',
        'pyyaml',
        'scikit-learn',
        'scipy',
    ],
    # Test dependencies
    # Execute 'python setup.py test' to run tests
    'test': [
        'flake8',
        'pytest>=3.9,!=6.0.0rc1,!=6.0.0',
        'pytest-cov>=2.10.1',
        'pytest-env',
        'pytest-flake8>=1.0.6',
        'pytest-html!=2.1.0',
        'pytest-metadata>=1.5.1',
        'pytest-mypy',
        'pytest-mock',
        'pytest-xdist',
    ],
    # Development dependencies
    # Use pip install -e .[develop] to install in development mode
    'develop': [
        'autodocsumm',
        'codespell',
        'docformatter',
        'isort',
        'pre-commit',
        'prospector[with_pyroma,with_mypy]>=1.9.0',
        'sphinx>2',
        'sphinx_rtd_theme',
        'vprof',
        'yamllint',
        'yapf',
    ],
}


def discover_python_files(paths, ignore):
    """Discover Python files."""
    def _ignore(path):
        """Return True if `path` should be ignored, False otherwise."""
        return any(re.match(pattern, path) for pattern in ignore)

    for path in sorted(set(paths)):
        for root, _, files in os.walk(path):
            if _ignore(path):
                continue
            for filename in files:
                filename = os.path.join(root, filename)
                if (filename.lower().endswith('.py')
                        and not _ignore(filename)):
                    yield filename


class CustomCommand(Command):
    """Custom Command class."""
    def install_deps_temp(self):
        """Try to temporarily install packages needed to run the command."""
        if self.distribution.install_requires:
            self.distribution.fetch_build_eggs(
                self.distribution.install_requires)
        if self.distribution.tests_require:
            self.distribution.fetch_build_eggs(self.distribution.tests_require)


class RunLinter(CustomCommand):
    """Class to run a linter and generate reports."""

    user_options: list = []

    def initialize_options(self):
        """Do nothing."""

    def finalize_options(self):
        """Do nothing."""

    def run(self):
        """Run prospector and generate a report."""
        check_paths = PACKAGES + [
            'setup.py',
            'tests',
        ]
        ignore = [
            'doc/',
        ]

        # try to install missing dependencies and import prospector
        try:
            from prospector.run import main
        except ImportError:
            # try to install and then import
            self.distribution.fetch_build_eggs(['prospector[with_pyroma]'])
            from prospector.run import main

        self.install_deps_temp()

        # run linter

        # change working directory to package root
        package_root = os.path.abspath(os.path.dirname(__file__))
        os.chdir(package_root)

        # write command line
        files = discover_python_files(check_paths, ignore)
        sys.argv = ['prospector']
        sys.argv.extend(files)

        # run prospector
        errno = main()

        sys.exit(errno)


def read_authors(filename):
    """Read the list of authors from .zenodo.json file."""
    with Path(filename).open() as file:
        info = json.load(file)
        authors = []
        for author in info['creators']:
            name = ' '.join(author['name'].split(',')[::-1]).strip()
            authors.append(name)
        return ', '.join(authors)


setup(
    name='BGCVal2',
    version="1.0",
    author=read_authors(".zenodo.json"),
    description="BGCVal2: modernized version of BGCVal",
    long_description=Path('README.md').read_text(),
    long_description_content_type='text/markdown',
    url='https://github.com/valeriupredoi/bgcval2',
    download_url='https://github.com/valeriupredoi/bgcval2',
    license='probably GPL',
    classifiers=[
        'Development Status :: 0 - initial',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: probably GPL',
        'Natural Language :: English',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Topic :: Scientific/Engineering :: GIS',
        'Topic :: Scientific/Engineering :: Hydrology',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    packages=PACKAGES,
    # Include all version controlled files
    include_package_data=True,
    setup_requires=REQUIREMENTS['setup'],
    install_requires=REQUIREMENTS['install'],
    tests_require=REQUIREMENTS['test'],
    extras_require={
        'develop': REQUIREMENTS['develop'] + REQUIREMENTS['test'],
        'test': REQUIREMENTS['test'],
    },
    entry_points={
        'console_scripts': [
            'analysis_compare = bgcval2.analysis_compare:main',
            'analysis_level3_amoc = bgcval2.analysis_level3_amoc:main',
            'analysis_level3_dms = bgcval2.analysis_level3_dms:main',
            'analysis_level3_omz = bgcval2.analysis_level3_omz:main',
            'analysis_level3_sameYear = bgcval2.analysis_level3_sameYear:main',
            'analysis_p2p = bgcval2.analysis_p2p:run',
            'analysis_timeseries = bgcval2.analysis_timeseries:main',
            'download_from_mass = bgcval2.download_from_mass:main',   
            'bgcval2_make_report = bgcval2.bgcval2_make_report:main',
        ],
    },
    cmdclass={
        #         'test': RunTests,
        'lint': RunLinter,
    },
    zip_safe=False,
)
