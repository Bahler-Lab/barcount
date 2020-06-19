from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()
    
    
setup(name='barcount',
      version='0.9',
      description='Python script for counting DNA barcodes in next generation sequencing data of pooled library experiments',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/Bahler-Lab/barcount',
      author='Stephan Kamrad',
      author_email='stephan.kamrad@crick.ac.uk',
      license='MIT',
      packages=['barcount'],
      install_requires=[
          'pandas',
          'numpy',
          'biopython',
         ],
      entry_points = {
        'console_scripts': [
            'barcount = barcount:cli',
                ],
            },
      classifiers=[
        "Development Status :: 4 - Beta", 
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        ],
      zip_safe = False,
      python_requires='>=3.7')
