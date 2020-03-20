#!/usr/bin/env python
from setuptools import setup, find_packages

install_requires = [
    'pyfaidx>=0.5.3.1',
    'pysam>=0.11.2.2y7',
    'pysamstats>=1.0.1'
]

tests_require = [
    'nose',
    'mock'
]

extras_require = {
    'docs': [
        'Sphinx>=1.1',
        'numpydoc>=0.5'
    ]
}


setup(name='mutbox',
      version='0.0.1',
      url='https://github.com/danielmsk/mutbox',
      license='MIT',
      author='Daniel Minseok Kwon',
      author_email='daniel.minseok.kwon@gmail.com',
      description='Genome portal site generator',
      download_url='https://github.com/danielmsk/mutbox/archive/0.1.tar.gz',
      keywords=['genomics', 'bioinformatics'],
      classifiers=[
          'Operating System :: OS Independent',
          'Topic :: Software Development :: Libraries',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
      ],
      # packages=find_packages(exclude=['tests']),
      packages=find_packages('src'),
      package_dir={'': 'src'},
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      zip_safe=False,
      install_requires=install_requires,
      setup_requires=['nose>=1.0'],
      test_suite='nose.collector',
      entry_points={
          'console_scripts': [
              'mutbox=mutbox:cli',
          ]
      })
