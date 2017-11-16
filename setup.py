from distutils.core import setup
import os

setup(name='UQuant',
      version='1.0',
      description='Uncertainty Quantification package',
      author='Soheil Hajian and Nikolai Strogies',
      url='https://',
      packages=['UQuant'],
      package_data={'UQuant': ['lib/CWrapper.so']}
     )
