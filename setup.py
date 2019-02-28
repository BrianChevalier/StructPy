#from setuptools import setup, find_packages
from distutils.core import setup

path = '/private/var/mobile/Containers/Shared/AppGroup/A79EC7C7-D188-4D79-9F2F-9CE43899AC1C/Pythonista3/Documents/site-packages-3'		

packages = ['Resources','StructPy']

setup(name='StructPy',
      version='0.1',
      description='Structural Analysis in Python',
      author='Brian Chevalier',
      author_email='Brian.Chevalier@gmail.com',
      url='https://github.com/BrianChevalier/StructPy',
      download_url='https://github.com/BrianChevalier/StructPy.git',
      license='MIT',
      install_requires=[
          'numpy',
          'matplotlib'],
      zip_safe=False,
      packages=packages)


