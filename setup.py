from setuptools import setup

setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='XAnalysis',
    url='https://github.com/Bribak/xanalysis',
    author='Daniel Bojar',
    author_email='daniel@bojar.net',
    # Needed to actually package something
    packages=['xanalysis'],
    # Needed for dependencies
    install_requires=['biopython','scikit-learn','pandas','numpy'],
    # *strongly* suggested for sharing
    version='0.2',
    # The license can be anything you like
    license='CC0',
    description='Analyzing DNA sequences',
    long_description=open('README.md').read(),
)
