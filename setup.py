from distutils.core import setup, Extension

setup(
    name='tgirt_map',
    version='0.1',
    description='TGIRT RNA-seq pipeline',
    url='',
    author='Douglas Wu',
    author_email='wckdouglas@gmail.com',
    license='MIT',
    packages=['tgirt_map'],
    scripts = ['bin/tgirt_count.py']
)
