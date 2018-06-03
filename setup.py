from setuptools import setup, find_packages

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
      packages=['StructPy.StructPy', 'StructPy.Resources'])
 
