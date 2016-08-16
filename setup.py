from setuptools import setup

setup(name='rapclust',
      version='0.1.1',
      scripts=['bin/RapClust'],
      description='Accurate, Fast and Lightweight Clustering of de novo Transcriptomes using Fragment Equivalence Classes',
      url='https://github.com/COMBINE-lab/RapClust',
      author='Avi Srivastava, Hirak Sarkar, Laraib Malik, Rob Patro',
      author_email='rob.patro@cs.stonybrook.edu',
      license='BSD with attribution',
      packages=['rapclust'],
      install_requires=[
          'PyYAML',
          'click',
          'networkx',
          'numpy',
          'pandas'
      ],
      zip_safe=False)
