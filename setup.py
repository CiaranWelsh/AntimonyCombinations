from setuptools import setup

MAJOR = 0
MINOR = 0
MICRO = 2

version = f'{MAJOR}.{MINOR}.{MICRO}'

setup(
    name='AntimonyCombinations',
    version=version,
    packages=['antimony_combinations'],
    license='MIT',
    long_description=open('README.md').read(),
    author='Ciaran Welsh',
    author_email='ciaran.welsh@newcastle.ac.uk',
    url='https://github.com/CiaranWelsh/AntimonyCombinations',
    keywords=['SBML', 'antimony', 'tellurium', 'model selection',
              'combinations', 'model merging'],
    install_requires=[
        'tellurium', 'pycotools3', 'pandas', 'numpy'
    ]
)
