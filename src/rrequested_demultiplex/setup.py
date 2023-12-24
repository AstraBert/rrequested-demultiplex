from setuptools import setup, find_packages

setup(
    name='rrequested_pkg',
    version='1.2.0',
    description='rrequested-pkg: A versatile package to preprocess raw basecalled reads with quality and size filtering and ex-novo demultiplexing',
    author='Astra Bertelli',
    author_email='astra.bertelli01@universitadipavia.it',
    url='https://github.com/AstraBert/rrequested-pkg',
    packages=find_packages(),
    install_requires=[
        'biopython==1.81',
        'edlib==1.3.9',
    ],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.10',    
    ],
    long_description=open('readme.md').read(),
    long_description_content_type='text/markdown',
)