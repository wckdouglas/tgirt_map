from distutils.core import setup, Extension
import glob

setup(
    name='tgirt_map',
    version='0.1',
    description='TGIRT RNA-seq pipeline',
    url='',
    author='Douglas Wu',
    author_email='wckdouglas@gmail.com',
    license='MIT',
    packages=['tgirt_map'],
    package_data={'tgirt_map': ['snakemake/*.smk']},
    scripts = glob.glob('bin/*.py')
)
