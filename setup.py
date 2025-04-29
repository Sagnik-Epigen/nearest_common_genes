from setuptools import setup, find_packages

setup(
    name='nearest_common_genes',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'gget',
        'pybedtools',
        'pandas',
    ],
    entry_points={
        'console_scripts': [
            'nearest-genes=nearest_common_genes.core:main',
        ],
    },
    author='Your Name',
    description='Fetch nearest genes from two BED files and find common ones using gget.',
    python_requires='>=3.7',
)
