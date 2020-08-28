#!/usr/bin/env python
from setuptools import setup, find_packages

install_requires = [
    'pyfaidx>=0.5.3.1',
    'pysam>=0.11.2.2y7',
    'pytabix>=0.0.2',
    'pyliftover>=0.4',
    'xmltodict>=0.12.0',
    'pyvcf>=0.6.8',
    'tqdm>=4.48.2'
]

tests_require = [
    'nose',
    'mock'
]

extras_require = {
    'docs': [
        'Sphinx>=1.1'
    ]
}

setup(name='mutanno',
      version='0.4.5',
      url='https://github.com/dbmi-bgm/mutanno',
      license='MIT',
      author='CGAP team at Harvard Medical School',
      author_email='daniel.minseok.kwon@gmail.com',
      description='Mutation Annotation Tool',
      download_url='https://github.com/dbmi-bgm/mutanno/archive/0.4.0.tar.gz',
      keywords=['genomics', 'bioinformatics'],
      classifiers=[
          'Operating System :: OS Independent',
          'Topic :: Software Development :: Libraries',
          'License :: OSI Approved :: MIT License',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
      ],
      # packages=find_packages(exclude=['tests']),
      python_requires=">=3.4",
      packages=find_packages('src'),
      package_dir={'': 'src'},
      package_data={'mutanno.data': ['*.json', '*.gz'],
                    'mutanno.web.templates': ['*.html'],
                    'mutanno.web.static.css': ['*.css'],
                    'mutanno.web.static.fonts': ['*.woff', '*.woff2'],
                    'mutanno.web.static.img': ['*.png', '*.gif'],
                    'mutanno.web.static.js': ['*.js', '*.png', '*.map'],
                    'mutanno.web.static.semantic': ['*.js', '*.json', '*.css', '*.md', 'LICENSE'],
                    'mutanno.web.static.semantic.components': ['*.js', '*.css'],
                    'mutanno.web.static.semantic.themes.default.assets.images': ['*.png'],
                    'mutanno.web.static.semantic.themes.default.assets.fonts': ['*.eot', '*.svg', '*.ttf', '*.woff', '*.woff2', '*.otf']
                    },
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      zip_safe=False,
      install_requires=install_requires,
      setup_requires=['nose>=1.0'],
      test_suite='nose.collector',
      entry_points={
          'console_scripts': [
              'mutanno=mutanno:cli',
          ]
      })
